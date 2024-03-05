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

MaxItr = 1000;
ErrorThr = 1e-3;     % RVM 终止误差
D_start=-30;
D_stop=29;
DOA11=-2:-2:-10;
DOA22=2:2:10;
DOA_train=[DOA11;DOA22]+0.1*rand(2,length(DOA11));
LM=size(DOA_train,2);
XX=zeros(LM,M,M);
S=zeros(LM,M,M);
theta=D_start:1:D_stop;
L=length(theta);
A=exp(1i*pi*fc*(0:M-1)'*sind(theta)/f0);
H=zeros(M*M,L);
for i=1:M
    fhi=A*diag(exp(-1i*pi*(i-1)*sind(theta)));
    H((i-1)*M+1:i*M,:)=fhi;
end

S_label=zeros(LM,L);
gamma=zeros(LM,L);
gamma_R=zeros(LM,L);

R_est=zeros(LM,C);
RS_est=zeros(LM,C);
S_est=zeros(LM,L,2);
SS_est=zeros(LM,L,2);
for i=1:LM
%     [X1,~]=signal_generate(M,snapshot,DOA_train(1,i),f0,fc,fs,1);
%     [X2,~]=signal_generate(M,snapshot,DOA_train(2,i),f0,fc,fs,1);
%     temp1=awgn(X1,0,'measured');
%     temp2=awgn(X2,0,'measured');
%     X= temp1+ temp2;
%     [R_est(i,:),Rx]=feature_extract_R(X) ;
      X=signal_generate_all(M,snapshot,DOA_train(:,i),f0,fc,fs,1);
      X=awgn(X,0,'measured');
      [R_est(i,:),Rx,RS,RS_est(i,:)]=feature_extract_R(X) ;%R_est严格上三角的实部和虚部，Rx为协方差矩阵
        %RS_est严格上三角的实部和虚部，RS为信号协方差矩阵
   
    temp=H'*vec(Rx);
    temp=temp/norm(temp);
    S_est(i,:,1)=real(temp);
    S_est(i,:,2)=imag(temp);

    XX(i,:,:)=Rx;%Rx为协方差矩阵保存
    S(i,:,:)=RS;%RS为信号协方差矩阵保存

    temp2=H'*vec(RS);
    temp2=temp2/norm(temp2);
    SS_est(i,:,1)=real(temp2);
    SS_est(i,:,2)=imag(temp2);
    
    
    S_label(i,round(DOA_train(1,i))+31)=1;
    S_label(i,round(DOA_train(2,i))+31)=1;
    
    [gamma_R(i,:)] = SBAC_R(X,H,MaxItr,ErrorThr, S_label(i,:),2);
    
    [gamma(i,:)] = SBAC(X,A,MaxItr,ErrorThr);
    
    
    i
end
% % %

close all
i=1;
plot(theta,S_est(i,:,1))
xlim([-30,30])
hold on
plot(theta,(S_label(i,:)'))
grid on
plot(theta,(gamma(i,:)'))
plot(theta,(gamma_R(i,:)'))
legend('A','true','S','R')


save('data_test_2.mat','R_est','DOA_train',...
    'S_label','S_est','SS_est','RS_est','XX','S','theta','gamma','gamma_R')

