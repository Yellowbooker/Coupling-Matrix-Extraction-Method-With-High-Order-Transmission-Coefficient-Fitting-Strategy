%% Cauchy Method With Transmission Zeros Selection

% =========================================================================

% It is a code for extracting the coupled-resonator circuit of a microwave
% filter. The method used to fit the sampled data is the Cauchy Method.

% Transmission zeros selecting operation is used to improve the fitting 
% accuracy of the sampled data.

% Operating steps;
% 1. Remove the phase loading by using high-order rational function.
% 2. Fit the sampled data by Cauchy Method.
% 3. Transmission zeros selection.
% 4. Calculate the transversal coupling matrix by using the Y-parameters.
% 5. Transform the transversal coupling matrix into the matrix with the
% target topology.

% Note: only suit the filter with the real-frequency transmission zeros.

% Data: 2023/06/01 By YellowBook

% Ref: A new computer-aided tuning scheme for general lossy 
% coupled-resonator bandpass filters based on the Cauchy method

% =========================================================================

clear
close all
warning off
% The warning is about the polt issue... ignore it

format short
Size = 18; % size of figure
TZrange = 0.25; % for selecting TZs with the realpart lower than range
% ======================= Loop Parameters =======================

X = 68; % Max number of iterations
err_ak_limit = -120;    % in dB
err_cp_limit = -120;	% in dB

% ======================= Filter Parameters =======================

N = 9; % Order of the filter
Nz = 3; % Number of the transmission zeros 
% If the extracted Tz is not correct,plz increase Nz and 
% create a new parasitic coupling to explain the additional Tz
CF = sqrt(1920*1980)*1e6; % Hz center mapping frequency
BW = (1980-1920)*1e6; % Hz mapping bandwidth
FBW = BW/CF;
M = 2; % Additional order of the rational function to remove the phase loading
Line = 10; % Radius to classify the poles and zeros

% % Set fitting frequency range in the lowpass range
w1 = -3.5;     
w2 = 3.5;

% ======================= Sampling data import =======================

[F_band,Origin_S11,Origin_S12,Origin_S22]=Reads2p();

Ns=length(F_band);% number of sampling data


Origin_dBS11(:,1)=20*log10(abs(Origin_S11(:,1)));
Origin_dBS12(:,1)=20*log10(abs(Origin_S12(:,1)));
Origin_dBS22(:,1)=20*log10(abs(Origin_S22(:,1)));
w_low=(F_band./CF-CF./F_band)/FBW;
s_low=1i*w_low;


[~,i1]=min(abs(w_low-ones(Ns,1)*w1));
[~,i2]=min(abs(w_low-ones(Ns,1)*w2));

% De_Embedding_2p is a function to fit S11 and S22 independently

[S11,S12,S22]=De_Embedding_2p_onlyarg(F_band,Origin_S11,Origin_S12,Origin_S22,N,CF,BW,M,Line);

% find the fitting range
Ns = i2 - i1 + 1;
S11 = S11(i1:i2);
S12 = S12(i1:i2);
S22 = S22(i1:i2);
s_low = s_low(i1:i2);
F_band = F_band(i1:i2);
Origin_dBS11 = Origin_dBS11(i1:i2);
Origin_dBS12 = Origin_dBS12(i1:i2);
Origin_dBS22 = Origin_dBS22(i1:i2);


% % 添加高斯白噪声，幅度值为-80dB--0.0001
% Mag_noise = 10^(-1*((140+3)/20)); %dB
% noise11 = Mag_noise*(randn(Ns,1) + 1i*randn(Ns,1));
% noise22 = Mag_noise*(randn(Ns,1) + 1i*randn(Ns,1));
% noise12 = Mag_noise*(randn(Ns,1) + 1i*randn(Ns,1));
% S11 = S11 + noise11;
% S22 = S22 + noise22;
% S12 = S12 + noise12;


M_S11 = diag(S11);
M_S22 = diag(S22);
M_S12 = diag(S12);

Vn = ones(Ns,N+1);
Vnz = ones(Ns,Nz+1);
for k = 1:Ns
    for index = 1:N
        Vn(k,index+1) = s_low(k)^(index);
    end
end
left=[M_S12*Vn        zeros(Ns,N+1)   -1.*M_S11*Vn;
      zeros(Ns,N+1)   M_S12*Vn        -1.*M_S22*Vn;];
[U,S,V] = svd(left);
fgp = V(:,end);
F = fgp(N+1:-1:1);
F22 = fgp(2*N+2:-1:N+2);
Pall = fgp(2*N+3+N:-1:2*N+3);

Tzall = roots(Pall);

for i = 1:N
    if abs(real(Tzall(i))) > TZrange
        Tz_t(i)  = 10;
    else
        Tz_t(i) = Tzall(i);
    end
end
Tz_t = sort(Tz_t);
Tz = Tz_t(1:Nz);

figure('name','TZ')
plot(real(Tzall),imag(Tzall),'bo','MarkerSize', 10,'linewidth',2);
hold on
plot(real(Tz),imag(Tz),'rs','MarkerSize', 15,'linewidth',2);
hold on
legend('拟合的零点','选择的零点', 'Location', 'NorthWest')
% legend('S11 Mea.','S12 Mea.','S11 CM.','S12 CM.', 'Location', 'NorthWest')
ylabel('虚部','fontsize',Size);
% set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
xlabel('实部','fontsize',Size);
% set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
% xlim([-3,3]);
% ylim([-3,3]);
% yticks([-3:1:3]);
% xticks([-3:1:3]);
grid on


P = poly(Tz);

F = F./F(1);
F22 = F22./F22(1);

% Fit = polyval(P,s_low)./polyval(F,s_low); % 估计值
% Sim = S12./S11; % 仿真值
% Er = Sim\Fit;

%确定系数,通过S21/S11 or S21/S22来确定

Fit = abs(polyval(P,s_low)./polyval(F,s_low)); % 估计值
Sim = abs(S12./S11); % 仿真值
Er = (Sim.'*Fit)/(Sim.'*Sim);
if mod(N-Nz,2) == 0
    Er = Er*1i;
end

P = P./Er;
P = [zeros(1,N-Nz), P].';
TE = conv(F,F22) - conv(P,P);
allpole = roots(TE);

% %确定系数,通过S21/S11 or S21/S22来确定 在0频率确定
% [mins,index] = min(abs(s_low));
% Fit = polyval(P,s_low(index))./polyval(F,s_low(index)); % 估计值
% Sim = S12(index)/S11(index); % 仿真值
% Er = Fit/Sim;
% P = P./Er;
% P = [zeros(1,N-Nz), P].';
% TE = conv(F,F22) - conv(P,P);
% allpole = roots(TE);


figure('name','K')
plot(imag(s_low),db(Sim),'b-',imag(s_low),db(Fit./Er),'k--','Linewidth',2);
% legend('Pole', 'NorthWest')
legend('Simulated K','Fitted K', 'NorthWest')
ylabel('Magnitude (dB)','fontsize',Size);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',Size);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
grid on



allpole = sort(allpole,'ComparisonMethod','real');
poleL = allpole(1:N);
poleR = allpole(N+1:2*N);

E = poly(poleL).';
T = poly(poleR).';


figure('name','P')
plot(real(poleL),imag(poleL),'bx',real(poleR),imag(poleR),'b+','Linewidth',1.5,'Markersize',10);
legend('Roots of E','Roots of T', 'NorthWest')
ylabel('Imaginary part','fontsize',Size);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
xlabel('Real part','fontsize',Size);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
grid on



% E = E.*F22(1);
% T = T.*F(1);

Yd = E + F + F22 + T;
Y11n = E - F + F22 - T;
Y12n = -2*P;
Y22n = E + F - F22 - T;

[residueY11,Ypole,~] = residue(Y11n,Yd);
[residueY12,Ypole,~] = residue(Y12n,Yd);
[residueY22,Ypole,~] = residue(Y22n,Yd);

% ======================= Y to TCM =======================
M=zeros(N+2,N+2);
for index=1:N
    M(index+1,index+1)=1i*Ypole(index,1);
    if abs(residueY11(index,1))>=abs(residueY22(index,1))
        M(1,index+1)=sqrt(residueY11(index,1));
        M(N+2,index+1)=residueY12(index,1)/sqrt(residueY11(index,1));
    else if abs(residueY11(index,1))<abs(residueY22(index,1))
        M(N+2,index+1)=sqrt(residueY22(index,1));
        M(1,index+1)=residueY12(index,1)/sqrt(residueY22(index,1));
        end
    end
    M(index+1,1)=M(1,index+1);
    M(index+1,N+2)=M(N+2,index+1);
end

CM=to_foldedCM(N,M);
fprintf('The extracted folded coupling matrix (real part):\n');
display(real(CM));  
fprintf('The extracted folded coupling matrix (imaginary part):\n');
display(imag(CM)); 

% % ======================= Analyse CM  =======================

[ExtrS11,ExtrS22,ExtrS12]=analyseCM(N,F_band,CM,CF,FBW,Ns);
dBExtrS11(1,:)=20*log10(abs(ExtrS11(1,:)));
dBExtrS22(1,:)=20*log10(abs(ExtrS22(1,:)));
dBExtrS12(1,:)=20*log10(abs(ExtrS12(1,:)));
figure('name','S-Parameters');

plot(F_band/10^6,db(S11),'-','Color',[0.75, 0.75, 0.75],'Linewidth',2);
hold on
plot(F_band/10^6,db(S12),'-','Color',[0.75 0.75 0.75],'Linewidth',2);
hold on
plot(F_band/10^6,dBExtrS11,'k:',F_band/10^6,dBExtrS12,'k--','Linewidth',2);
hold on
legend('S11 仿真','S12 仿真','S11 提取','S12 提取', 'Location', 'NorthWest')
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
ylabel('S参数幅度 (dB)','fontsize',Size);
% set(gca,'FontName','Times New Roman');
set(gca,'linewidth',1.2);
xlabel('频率 (MHz)','fontsize',Size);
% set(gca,'FontName','Times New Roman');
set(gca,'linewidth',1.2);
% xlim([-5,5]);
ylim([-150,0]);
% xticks([-5:1:5]);
yticks([-150:20:0]);
% grid on
xlim([F_band(1)/10^6,F_band(Ns)/10^6]);
% grid on
