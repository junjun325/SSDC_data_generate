%% 产生数据
clc
clear variables
close all
M=8;
f0=1e6;
fc=1e6;
fs=4*f0;
times=500;
MaxItr = 800;
ErrorThr = 1e-3;     % RVM 终止误差

D_start=-30;
D_stop=29;
K=2;
theta=D_start:1:D_stop; %-30到29共60个角度间隔
L=length(theta);
A=exp(1i*pi*fc*(0:M-1)'*sind(theta)/f0);
H=zeros(M*M,L);
for i=1:M
    fhi=A*diag(exp(-1i*pi*(i-1)*sind(theta)));
    H((i-1)*M+1:i*M,:)=fhi;
end
snapshot=50:50:400;
SNR=15;
S_label=zeros(times,L,length(snapshot));
R_est=zeros(times,M*(M-1),length(snapshot));
RS_est=zeros(times,M*(M-1),length(snapshot));
gamma=zeros(times,L,length(snapshot));
gamma_R=zeros(times,L,length(snapshot));
T_SBC=zeros(times,length(snapshot));
T_SBC_R=zeros(times,length(snapshot));
S_est=zeros(times,L,2,length(snapshot));
SS_est=zeros(times,L,2,length(snapshot));
DOA_train=zeros(2,times,length(snapshot));
XX=zeros(times,M,M,length(snapshot));
S=zeros(times,M,M,length(snapshot));

for k=1:length(snapshot)
    for i=1:times
%         DOA_train(1,i,k)=-12+1*rand; %cita1，cita2 左右
%         DOA_train(2,i,k)=DOA_train(1,i,k)+15;
        DOA_train(1,i,k)=-5.5; %cita1，cita2 左右
        DOA_train(2,i,k)=8.5;
        
%         [X1,~]=signal_generate(M,snapshot,DOA_train(1,i,k),f0,fc,fs,1);
%         [X2,~]=signal_generate(M,snapshot,DOA_train(2,i,k),f0,fc,fs,1);
%         temp1=awgn(X1,SNR(k),'measured');
%         temp2=awgn(X2,SNR(k),'measured');
%         X= temp1+ temp2;
       X=signal_generate_all_snapshot(M,snapshot,DOA_train(:,i,k),f0,fc,fs,1,k);
       X=awgn(X,SNR,'measured');
        [R_est(i,:,k),Rx,RS,RS_est(i,:,k)]=feature_extract_R(X) ;%R_est严格上三角的实部和虚部，Rx为协方差矩阵
        %RS_est严格上三角的实部和虚部，RS为信号协方差矩阵
        temp=H'*vec(Rx);
        temp=temp/norm(temp);
        S_est(i,:,1,k)=real(temp);
        S_est(i,:,2,k)=imag(temp);%计算出每个样本为60*2的输入，与RX对应
        
        XX(i,:,:,k)=Rx;%Rx为协方差矩阵保存
        S(i,:,:,k)=RS;%RS为信号协方差矩阵保存

        temp2=H'*vec(RS);
        temp2=temp2/norm(temp2);
        SS_est(i,:,1,k)=real(temp2);
        SS_est(i,:,2,k)=imag(temp2);%计算出每个样本为60*2的输入，与RS对应

        S_label(i,round(DOA_train(1,i,k))+31,k)=1;
        S_label(i,round(DOA_train(2,i,k))+31,k)=1;
        tic
        [gamma_R(i,:,k)] = SBAC_R(X,H,MaxItr,ErrorThr, S_label(i,:),2);
        T_SBC_R(i,k)=toc;
        tic
        [gamma(i,:,k)] = SBAC(X,A,MaxItr,ErrorThr);
        T_SBC(i,k)=toc;
        
        i
        k
    end
   
end

close all
plot(theta,S_est(i,:,1,k))
plot(theta,SS_est(i,:,1,k))
xlim([-30,30])
hold on
plot(theta,(S_label(i,:,k)'))
grid on
plot(theta,(gamma(i,:,k)'))
plot(theta,(gamma_R(i,:,k)'))
legend('A','true','S','R')
save('data_snapshot_500mc_snr15.mat','R_est','DOA_train','RS_est','XX','S',...
    'S_label','S_est','SS_est','theta','SNR','gamma','gamma_R','T_SBC','T_SBC_R')


