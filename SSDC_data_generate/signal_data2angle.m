%% 产生数据
clc
clear variables
close all
M=8;snapshot=256;
f0=1e6;
fc=1e6;
fs=4*f0;
C=M*(M-1);
%%
MaxItr = 800;
ErrorThr = 1e-3;     % RVM 终止误差
D_start=-30;
D_stop=29;
LM=500;
theta=D_start:1:D_stop;
L=length(theta);
A=exp(1i*pi*fc*(0:M-1)'*sind(theta)/f0);
H=zeros(M*M,L);
for i=1:M
    fhi=A*diag(exp(-1i*pi*(i-1)*sind(theta)));
    H((i-1)*M+1:i*M,:)=fhi;
end
Angle=0:2:40;
S_label=zeros(LM,L,length(Angle));
gamma=zeros(LM,L,length(Angle));
gamma_R=zeros(LM,L,length(Angle));
R_est=zeros(LM,C,length(Angle));
RS_est=zeros(LM,C,length(Angle));
S_est=zeros(LM,L,2,length(Angle));
SS_est=zeros(LM,L,2,length(Angle));
DOA_train=zeros(2,LM,length(Angle));
S_abs=zeros(LM,2*L,length(Angle));
SS_abs=zeros(LM,2*L,length(Angle));
T_SBC=zeros(LM,length(Angle));
T_SBC_R=zeros(LM,length(Angle));
XX=zeros(LM,M,M,length(Angle));
S=zeros(LM,M,M,length(Angle));

for k=1:length(Angle)
    for i=1:LM
%         DOA_train(1,i,k)=-Angle(k)/2+1*rand-0.5;%在固定的范围内随机生成的一个角度
%         DOA_train(2,i,k)=DOA_train(1,i,k)+Angle(k);%保持两个角度是固定的间隔
        DOA_train(1,i,k)=-15.5;%
        DOA_train(2,i,k)=DOA_train(1,i,k)+Angle(k);%保持两个角度是固定的间隔
        X=signal_generate_all(M,snapshot,DOA_train(:,i,k),f0,fc,fs,1);
        X=awgn(X,0,'measured');
        [R_est(i,:,k),Rx,RS,RS_est(i,:,k)]=feature_extract_R(X) ;%R_est严格上三角的实部和虚部，Rx为协方差矩阵
        %RS_est严格上三角的实部和虚部，RS为信号协方差矩阵
%    
%         [X1,~]=signal_generate(M,snapshot,DOA_train(1,i,k),f0,fc,fs,1);
%         [X2,~]=signal_generate(M,snapshot,DOA_train(2,i,k),f0,fc,fs,1);
%         temp1=awgn(X1,0,'measured');
%         temp2=awgn(X2,0,'measured');
%         X= temp1+ temp2;
%         [R_est(i,:,k),Rx]=feature_extract_R(X) ;
        temp=H'*vec(Rx);
        temp=temp/norm(temp);
        S_est(i,:,1,k)=real(temp);
        S_est(i,:,2,k)=imag(temp);
 
        XX(i,:,:,k)=Rx;%Rx为协方差矩阵保存
        S(i,:,:,k)=RS;%RS为信号协方差矩阵保存


        temp2=H'*vec(RS);
        temp2=temp2/norm(temp2);
        SS_est(i,:,1,k)=real(temp2);
        SS_est(i,:,2,k)=imag(temp2);%计算出每个样本为60*2的输入，与RS对应

        S_abs(i,:,k)=[real(temp);imag(temp)];
        SS_abs(i,:,k)=[real(temp2);imag(temp2)];
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
xlim([-60,60])
hold on
plot(theta,(S_label(i,:,k)'))
grid on
plot(theta,(gamma(i,:,k)'))
plot(theta,(gamma_R(i,:,k)'))
legend('A','true','S','R')


save('data_angle_20_2.mat','R_est','RS_est','XX','S','DOA_train','gamma_R',...
    'S_label','S_est','SS_est','S_abs','SS_abs','theta','gamma','T_SBC','T_SBC_R','Angle')
%data2_angle1  一个角度固定为0
%data2_angle  法线对称
figure
 plot(mean(T_SBC))
 hold on
plot(mean(T_SBC_R))
