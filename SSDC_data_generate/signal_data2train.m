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
DOA11=[];
DOA22=[];
k1=[1:1:40];
k=repmat(k1,1,10);
D_start=-30;
D_stop=29;
for l=1:length(k)
    DOA1=D_start:1:D_stop-k(l);
    DOA2=D_start+k(l):1:D_stop;
    DOA11=[DOA11,DOA1];
    DOA22=[DOA22,DOA2];
end
DOA_train=[DOA11;DOA22];

theta=D_start:1:D_stop;
L=length(theta);
A=exp(1i*pi*fc*(0:M-1)'*sind(theta)/f0);
H=zeros(M*M,L);
for i=1:M
    fhi=A*diag(exp(-1i*pi*(i-1)*sind(theta)));
    H((i-1)*M+1:i*M,:)=fhi;
end
for ii=1:L
    aa = A(:,ii);
    conjaa = conj(aa);
    H2(:,ii) = kron(conjaa,aa);
end


S_label=zeros(length(DOA_train),L);
RS_est=zeros(length(DOA_train),C);
RRS_est=zeros(length(DOA_train),C);
R_est=zeros(length(DOA_train),C);
XX=zeros(length(DOA_train),M,M);
S=zeros(length(DOA_train),M,M);
RRRS=zeros(length(DOA_train),M,M);
S_est=zeros(length(DOA_train),L,2);
SS_est=zeros(length(DOA_train),L,2);
RRRS_est=zeros(length(DOA_train),L,2);
S_abs=zeros(length(DOA_train),2*L);
SS_abs=zeros(length(DOA_train),2*L);
RRRS_abs=zeros(length(DOA_train),2*L);
for i=1:length(DOA_train)
    i
    X=signal_generate_all(M,snapshot,DOA_train(:,i),f0,fc,fs,1);
    X=awgn(X,-20*rand,'measured');
    [R_est(i,:),Rx,RS,RS_est(i,:),RRS,RRS_est(i,:)]=feature_extract_R(X) ;%R_est严格上三角的实部和虚部，Rx为协方差矩阵
        %RS_est严格上三角的实部和虚部，RS为信号协方差矩阵,RRS为前后项信号协方差矩阵的和求平均,RRS_est为RRS的严格上三角的实部和虚部
%     [X1,~]=signal_generate(M,snapshot,DOA_train(1,i),f0,fc,fs,1);
%     [X2,~]=signal_generate(M,snapshot,DOA_train(2,i),f0,fc,fs,1);
%     temp1=awgn(X1,-10*rand,'measured');
%     temp2=awgn(X2,-10*rand,'measured');
%     X= temp1+ temp2;
%     [R_est(i,:),Rx]=feature_extract_R(X) ;
     XX(i,:,:)=Rx;%Rx为协方差矩阵保存
     S(i,:,:)=RS;%RS为信号协方差矩阵保存
     RRRS(i,:,:)=RRS;%RRRS为信号协方差RRS矩阵保存
    temp=H'*vec(Rx);
    temp=temp/norm(temp);
    S_est(i,:,1)=real(temp);
    S_est(i,:,2)=imag(temp);

     temp2=H'*vec(RS);
     temp2=temp2/norm(temp2);
     SS_est(i,:,1)=real(temp2);
     SS_est(i,:,2)=imag(temp2);%计算出每个样本为60*2的输入，与RS对应

      temp3=H'*vec(RRS);
     temp3=temp3/norm(temp3);
     RRRS_est(i,:,1)=real(temp3);
     RRRS_est(i,:,2)=imag(temp3);%计算出每个样本为60*2的输入，与RRS对应

    S_abs(i,:)=[real(temp);imag(temp)];
    SS_abs(i,:)=[real(temp2);imag(temp2)];
    RRRS_abs(i,:)=[real(temp3);imag(temp3)];
    S_label(i,round(DOA11(i))+31)=1;
    S_label(i,round(DOA22(i))+31)=1;

end
% % %

close all
plot(theta,S_est(i,:,1))
plot(theta,SS_est(i,:,1))
xlim([-30,30])
hold on
plot(theta,(S_label(i,:)'))
grid on
legend('Rx-est','true','RS-est','S_abs','SS_abs')


save('data_train_data_lowsnr_new_snr_20.mat','R_est','DOA_train','RS_est','XX','S',...
    'S_label','S_est','SS_est','S_abs','SS_abs','RRRS','RRS_est','RRRS_est','RRRS_abs')

