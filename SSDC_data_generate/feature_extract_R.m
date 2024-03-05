function [r_doa1,Rx,RS,rs_doa1,RRS,rrs_doa1]=feature_extract_R(X)
%clc; clear all;X=randn(5,5);
[M,N]=size(X);
Rx=X*X'/N;%计算阵列协方差矩阵 时域协方差矩阵
r_doa=zeros(1,M*(M-1)/2);
rs_doa=zeros(1,M*(M-1)/2);
rrs_doa=zeros(1,M*(M-1)/2);
%%  协方差矩阵上三角矩阵
k=1;
for i=1:M
    for j=i+1:M
        r_doa(k)=Rx(i,j);
        k=k+1;
    end
end

% r_doa=r_doa/norm(r_doa,2);
% r_doa1=1*[real(r_doa) imag(r_doa)];

%
% r_doa=r_doa-min(real(r_doa))-1j*min(imag(r_doa));
r_doa=r_doa/norm(r_doa);
r_doa1=1*[real(r_doa) imag(r_doa)];
% temp=H'*vec(Rx);
% SS=temp/norm(temp);
% SS1=1*[real(SS) imag(SS)];
% ;
% r_doa=r_doa-mean(r_doa);
% r_doa=r_doa/norm(r_doa,2);
% r_doa1=[real(r_doa) imag(r_doa)];
%r_doa1=[angle(r_doa) abs(r_doa)].';

[EV,D]=eig(Rx);%%%%
EVA=diag(D)';
[EVA,I]=sort(EVA);
EVA=fliplr(EVA);
EV=fliplr(EV(:,I));
RS=EV(:,1)*diag(EVA(1))*EV(:,1)';%信号协方差矩阵
%%  协方差矩阵上三角矩阵
k=1;
for i=1:M
    for j=i+1:M
        rs_doa(k)=RS(i,j);
        k=k+1;
    end
end

rs_doa=rs_doa/norm(rs_doa);
rs_doa1=1*[real(rs_doa) imag(rs_doa)];
%%
J=fliplr(eye(M,M));
RRS=(RS+J*conj(RS)*J')./2;%（信号协方差矩阵+反向信号协方差矩阵）/2
%%  协方差矩阵上三角矩阵
k=1;
for i=1:M
    for j=i+1:M
        rrs_doa(k)=RRS(i,j);
        k=k+1;
    end
end

rrs_doa=rrs_doa/norm(rrs_doa);
rrs_doa1=1*[real(rrs_doa) imag(rrs_doa)];
