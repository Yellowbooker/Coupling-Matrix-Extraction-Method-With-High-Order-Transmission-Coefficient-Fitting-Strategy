% =========================================================================

% Calculate the S-parameters from a N+2 order coupling matrix

% Data: 2020/03/21 By YB

% =========================================================================

function [S11,S22,S12]=analyseCM(n,F_band,M,CF,FBW,N_s)
% N,N_pole
% Fstr=Fd(1,1)
% Fsto=Fd(N_sample,1)
% M= coupling matrix
% N_s= N_sample

R=zeros(n+2,n+2);
R(1,1)=1;
R(n+2,n+2)=1;
W=eye(n+2,n+2);
W(1,1)=0;
W(n+2,n+2)=0;
for k=1:N_s
    wp(k)=(F_band(k)/CF-CF/F_band(k))/FBW*(1i);
    A=R+W.*wp(k)+(1i).*M;
    AP=inv(A);
    S12(k)=2*AP(n+2,1);
    S11(k)=-1+2*AP(1,1);
    S22(k)=-1+2*AP(n+2,n+2);
end
end