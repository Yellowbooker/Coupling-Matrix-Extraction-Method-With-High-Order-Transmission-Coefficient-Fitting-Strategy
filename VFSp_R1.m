%% Vector Fitting for S-parameters_R1

% =========================================================================

% It is a code for extracting the coupled-resonator circuit of a microwave
% filter. The method used to fit the sampled data is the vector fitting (VF).

% Operating steps;
% 1. Remove the phase loading by using high-order rational function.
% 2. Fit the sampled data by the VF.
% 3. Transmission zeros selection.
% 4. Calculate the transversal coupling matrix by using the Y-parameters.
% 5. Transform the transversal coupling matrix into the matrix with the
% target topology.

% Note: only suit the filter with the real-frequency transmission zeros.

% Data: 2023/06/01 By YellowBook

% Ref1: Circuit Model Extraction of Parallel-Connected Dual-Passband
% Coupled-Resonator Filters
% Ref2: A new computer-aided tuning scheme for general lossy 
% coupled-resonator bandpass filters based on the Cauchy method

% Update Record (2024/04/14 By YellowBook): 
% 1. The weight of sqrt (P/F) is used in the calculation of Er(the ripple factor). 
% This weight is used to increase the effect of the small S11 data within 
% the passband. The default weight is set to 0.5, it can be revised in Filter 
% Parameters.
% 2. A new iterative technique is applied to identify poles, it may can be
% used when the order is higher than 15. The technique is inspired by Ref3.
% Ref3: Improving Accuracy in Solving Feldtkeller Equation

% =========================================================================


clear
close all
warning off
% The warning is about the polt issue... ignore it


format short
Size = 18; % size of figure
TZrange = 0.25; % for selecting TZs with the realpart lower than range
% ======================= Loop Parameters =======================

X = 18; % Max number of iterations
err_ak_limit = -60;    % in dB
err_cp_limit = -60;	% in dB

% ======================= Filter Parameters =======================

%--------------------------------------------------------------------------
% Example 1 (simualted data of a 9th-order filter with 3 Tzs)：
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
weight = 0.5;
%--------------------------------------------------------------------------
% Example 2 (measured data of a 6th-order filter with 2 Tzs)：
% N = 6; % Order of the filter
% Nz = 2; % Number of the transmission zeros 
% % If the extracted Tz is not correct,plz increase Nz and 
% % create a new parasitic coupling to explain the additional Tz
% CF = sqrt(1920*1980)*1e6; % Hz center mapping frequency
% BW = (1980-1920)*1e6; % Hz mapping bandwidth
% FBW = BW/CF;
% M = 1; % Additional order of the rational function to remove the phase loading
% Line = 10; % Radius to classify the poles and zeros
% % % Set fitting frequency range in the lowpass range
% w1 = -4;     
% w2 = 4;
% weight = 0.5;
%--------------------------------------------------------------------------

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
% Origin_S11 = Origin_S11(i1:i2);
% Origin_S12 = Origin_S12(i1:i2);
% Origin_S22 = Origin_S22(i1:i2);
S11 = S11(i1:i2);
S12 = S12(i1:i2);
S22 = S22(i1:i2);
s_low = s_low(i1:i2);
F_band = F_band(i1:i2);
Origin_dBS11 = Origin_dBS11(i1:i2);
Origin_dBS12 = Origin_dBS12(i1:i2);
Origin_dBS22 = Origin_dBS22(i1:i2);

% [S11,S12,S22]=De_Embedding_2p_onlyarg(F_band,Origin_S11,Origin_S12,Origin_S22,N,CF,BW,M,Line);
% % 添加高斯白噪声，幅度值为-80dB--0.0001
% Mag_noise = 10^(-1*((100+3)/20)); %dB
% noise11 = Mag_noise*(randn(Ns,1) + 1i*randn(Ns,1));
% noise22 = Mag_noise*(randn(Ns,1) + 1i*randn(Ns,1));
% noise12 = Mag_noise*(randn(Ns,1) + 1i*randn(Ns,1));
% S11 = S11 + noise11;
% S22 = S22 + noise22;
% S12 = S12 + noise12;

% ======================= Initial pole ak =======================

for index=1:N
    ak(index,1)=-1+2.1/(N-1)*(index-1);
    ak(index,1)=-0.01*abs(ak(index,1))+ak(index,1)*1i;
end

% ======================= Iterate Max X times =======================

A=zeros(N,N);
b=ones(N,1);A1=zeros(Ns,N+1);A2=zeros(Ns,N);A3=zeros(Ns,Nz+1);
M_S11 = diag(S11);
M_S12 = diag(S12);
M_S22 = diag(S22);


% Weight of the fitting
W12 = diag(1./(abs(S12).^0.5));
W11 = eye(Ns);
W22 = eye(Ns);

tmp_ak = ak;
cp_plot = zeros(N, X);
err_ctl = zeros(1, X);

for index=1:X
    
    A_den=ones(Ns,1);
    
%     indexlimit = real(ak)>0;
%     ak(indexlimit) = -conj(ak(indexlimit));
    
    for index1 = 1:N
        A2(:,index1) = 1./(s_low - ones(Ns,1).*ak(index1));
        A_den(:,1) = A2(:,index1).*A_den(:,1);
    end
    
    for index1 = 1:N+1
        A3(:,index1) = (s_low.^(index1-1)).*A_den;
    end
 
    A1 = [A2,ones(Ns,1)];
    A = diag(ak);
%============================ Original=====================================    
    left=[W11*A1           zeros(Ns,N+1)   zeros(Ns,N+1)   -1*W11*M_S11*A2;
          zeros(Ns,N+1)    W12*A3          zeros(Ns,N+1)   -1*W12*M_S12*A2;
          zeros(Ns,N+1)    zeros(Ns,N+1)   W22*A1          -1*W22*M_S22*A2];
    right=[W11*S11;
           W12*S12;
           W22*S22];
    Call=lsqminnorm(left,right); % This function can be replaced by \
    
    C=mat2cell(Call,[N+1,N+1,N+1,N],1);
    cp=C{4,1};
    cp_plot(:,index)=10*log10(abs(cp));
    ak=eig(A-b*cp.');
	err_ctl(index) = sum(abs(ak - tmp_ak));
	if (10*log10(err_ctl(index)) < err_ak_limit) && (max(cp_plot(:,index)) < err_cp_limit)
		break
	end
	tmp_ak = ak;
%==============================end=========================================
end

%==============================	Disp loop time	======
disp('Loop time:')
disp(index)

figure('name','MVF-error')
subplot(1,2,1);
plot([1:index],max(cp_plot(:,[1:index])),'b','Linewidth',1.5);
legend('cp', 'Location', 'NorthEast')
title('Cp error in Log')
grid on

subplot(1,2,2);
plot([1:index],10*log10(err_ctl([1:index])),'r','Linewidth',1.5);
legend('ak', 'Location', 'NorthEast')
title('ak error in Log')
grid on

C11=C{1,1};
C12=C{2,1};
C22=C{3,1};

% ======================= Calculate ExtrS =======================

ES11 = A1*C11;
ES12 = A3*C12;
ES22 = A1*C22;


dBES11(:,1)=20*log10(abs(ES11(:,1))); dBES12(:,1)=20*log10(abs(ES12(:,1)));
dBES22(:,1)=20*log10(abs(ES22(:,1)));

figure('name','S')
plot(imag(s_low),Origin_dBS11,'r',imag(s_low),Origin_dBS12,'b',imag(s_low),dBES11,'r:',imag(s_low),dBES12,'b:','linewidth',2);
legend('S11 Sim.','S12 Sim.','S11 Fit.','S12 Fit.', 'Location', 'NorthWest')
% legend('S11 Mea.','S12 Mea.','S11 CM.','S12 CM.', 'Location', 'NorthWest')
ylabel('S-parameters (dB)','fontsize',Size);
% set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
xlabel('\omega (rad/s)','fontsize',Size);
% set(gca,'FontName','Times New Roman');
set(gca,'FontSize',Size);
set(gca,'linewidth',1.2);
% xlim([-5,5]);
ylim([-150,0]);
yticks([-150:20:0]);
xlim([s_low(1)/1i,s_low(Ns)/1i]);
xticks([-2:1:2]);
grid on

Tzall = roots(flipud(C12));

% Selecting TZ

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

% 通过留数与极点确定零点
A = diag(ak);
Rz11 = eig(A-b*C11(1:end-1).'./C11(end));
Rz22 = eig(A-b*C22(1:end-1).'./C22(end));



F = poly(Rz11);
F22 = poly(Rz22);
P = poly(Tz);


Fit = (abs(polyval(P,s_low)./polyval(F,s_low))).^weight; % 估计值
Sim = (abs(S12./S11)).^weight; % 仿真值
Er = ((Sim.'*Fit)/(Sim.'*Sim)).^(1/weight);
if mod(N-Nz,2) == 0
    Er = Er*1i;
end


P = P./Er;

P = [zeros(1,N-Nz), P];
TE = conv(F,F22) - conv(P,P);

%==========================================================================
% % New iterative technique to identify poles 
% (It can be used when the order is higher than 15)

% pk = roots(TE);
% % pk = [ak;conj(ak);];
% % Q = poly(pk);
% % [ck,~,~] = residue(TE,Q);
% dk = zeros(length(pk),0);
% A = eye(length(pk));
% b = ones(length(pk),1);
% for i = 1:24
% %     Q = poly(pk);
%     [ck,pk,~] = residue(TE,poly(pk));
% %     dk = ck;
%     A = diag(pk);
%     pk = eig(A - b*ck.'); % recalculate poles
% end
% 
% allpole = pk;
%==========================================================================
% Old iterative technique to identify poles
allpole = roots(TE);
%==========================================================================

figure('name','K')
plot(imag(s_low),db((Sim).^(1/weight)),'b-',imag(s_low),db((Fit).^(1/weight)./Er),'k--','Linewidth',2);
% plot(imag(s_low),db((Sim).^(1/1)),'b-',imag(s_low),db((Fit).^(1/1)./Er),'k--','Linewidth',2);
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

E = poly(poleL);
T = poly(poleR);


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
CM = to_foldedCM(N,M);

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

figure()
plot(F_band/10^6,db(S22),'-','Color',[0.75, 0.75, 0.75],'Linewidth',2);
hold on
plot(F_band/10^6,dBExtrS22,'k--','Linewidth',2);
hold on
legend('S22 仿真','S22 提取', 'Location', 'NorthWest')
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
