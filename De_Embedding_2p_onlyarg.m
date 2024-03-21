%% A function to remove the phase loading at the ports of the coupled-resonator filter

%==========================================================================

% Input:
% F_band is the sampled frequency points
% S11, S12, S22 are the sampled S-parameters data
% N is the order of the filter
% CF is the center frequency of the filter
% BW is the bandwidth of the filter
% M is the additional order of the raditon function that used to repersent
% the loaded transmission line and constant phase shift
% Line is the radius that distinguish the poles and zeros

% Output:
% New_S11, New_S12 and New_S22 are S-parameters data for which the loaded 
% phase has been removed

% Ref: Phase De-Embedding of Narrowband CoupledResonator Networks by Vector Fitting

% Data: 2020/03/27 by YB

%==========================================================================

function [New_S11,New_S12,New_S22]=De_Embedding_2p_onlyarg(F_band,S11,S12,S22,N,CF,BW,M,Line)

X = 18; % maximum of the loop times
err_ak_limit = -60;  %in dB
err_cp_limit = -60;	% in dB

FBW = BW/CF;
N_s = length(F_band); % number of sampling data

w_low = (F_band./CF - CF./F_band)/FBW;
s_low = 1i*w_low;

A = zeros(N+M,N+M);
b = ones(N+M,1);
A1 = zeros(N_s,N+M+1);
A2 = zeros(N_s,N+M);
Sxx{1} = S11;
Sxx{2} = S22;

for INDEX = 1:length(Sxx)
    
    M_Sxx{INDEX} = diag(Sxx{INDEX});
    
    % ======================= Initial Common Pole ak [-3,2.3]

    for index = 1:N+M
        ak(index,1) = -3 + 5.3/(N + M - 1)*(index - 1);
        ak(index,1) = -0.01*abs(ak(index,1)) + ak(index,1)*1i;
    end
    tmp_ak = ak;
    cp_plot = zeros(N+M, X);
    err_ctl = zeros(1, X);

    for index = 1:X % loop of the VF
        
        for index1 = 1:N+M
            A2(:,index1) = 1./(s_low - ones(N_s,1).*ak(index1));
        end
        
        A1 = [A2,ones(N_s,1)];
        A = diag(ak,0);
        left = [A1   -1*M_Sxx{INDEX}*A2];
        right = Sxx{INDEX};
    %     Call = left\right;
        C_all = lsqminnorm(left,right); % Matlab R2018b
        C = mat2cell(C_all,[N+M,1,N+M],1);
        cp = C{3,1};
        cp_plot(:,index)=10*log10(abs(cp));
        ak = eig(A-b*cp.');
        err_ctl(index) = sum(abs(ak - tmp_ak));
        if 10*log10(err_ctl(index)) < err_ak_limit && max(cp_plot(:,index)) < err_cp_limit
            break
        end
        tmp_ak = ak;
    end
    disp('Loop time:')
    disp(index)
    pole{INDEX} = ak;
    residue{INDEX} = C{1,1};
    d{INDEX} = C{2,1};
    A11 = diag(pole{INDEX},0);
    zero{INDEX} = eig(A11-b*residue{INDEX}.'/d{INDEX});
    % ======================= Caculate Phase Factor
    index2 = 0;
    index3 = 0;
    for index=1:N+M
            if abs(zero{INDEX}(index,1))>=Line
                    index2 = index2+1;
                    alfa_zk{INDEX}(index2,1) = zero{INDEX}(index,1);
            end
            if abs(pole{INDEX}(index,1))>=Line
                    index3 = index3+1;
                    alfa_ak{INDEX}(index3,1) = pole{INDEX}(index,1);
            end
    end
    if index2 == 0 && index3 == 0
        for index = 1:N_s
        Extralfa{INDEX}(index,1) =1;
        end
    else
    TestNum = poly(alfa_zk{INDEX});
    TestDen = poly(alfa_ak{INDEX});
    for index = 1:N_s
        Extralfa{INDEX}(index,1) = d{INDEX}*polyval(TestNum,s_low(index,1))/polyval(TestDen,s_low(index,1));
        ExtrSxx{INDEX}(index,1) = d{INDEX}*polyval(poly(zero{INDEX}),s_low(index,1))/polyval(poly(pole{INDEX}),s_low(index,1));
    end
    end
    unwrap_Palfa{INDEX} = unwrap(angle(Extralfa{INDEX}));
    
end

% Size = 22;
% figure('name','P-Z of S11')
% plot(real(zero{1}),imag(zero{1}),'ro',real(pole{1}),imag(pole{1}),'bx','Linewidth',1.5,'MarkerSize',5);
% legend('Zero S11', 'Pole S11', 'NorthWest')
% set(gca,'FontName','Times New Roman');
% set(gca,'FontSize',24);
% ylabel('Imaginary part','fontsize',Size);
% set(gca,'FontName','Times New Roman');
% set(gca,'FontSize',Size);
% set(gca,'linewidth',1.2);
% xlabel('Real part','fontsize',Size);
% set(gca,'FontName','Times New Roman');
% set(gca,'FontSize',Size);
% set(gca,'linewidth',1.2);
% grid on
% % ylim([0 1]);
% % yticks(0:0.2:1);
% xlim([2 8]);
% % xticks(2:0.2:4);



% figure('name','P-Z of S22')
% plot(real(zero{2}),imag(zero{2}),'ro',real(pole{2}),imag(pole{2}),'bx','Linewidth',1.5);
% legend('Zero S22', 'Pole S22', 'NorthWest')
% figure('name','S11&ES11')
% plot(w_low,abs(S11),'y',w_low,abs(ExtrSxx{1}),'b--','Linewidth',1.5);
% figure('name','S22&ES22')
% plot(w_low,abs(S22),'y',w_low,abs(ExtrSxx{2}),'b--','Linewidth',1.5);






% ======================= Unwrap Phase In Matlab Three Ports

New_S11 = S11.*exp(-1i*unwrap_Palfa{1});
New_S22 = S22.*exp(-1i*unwrap_Palfa{2});

New_S12 = S12.*exp(-1i*(unwrap_Palfa{1}*0.5+unwrap_Palfa{2}*0.5));

% figure('name','After')
% plot(w_low,angle(New_S11),'r',w_low,angle(New_S22),'b','LineWidth',1.5);
% title('S11&S22 in Phase(after de embedding)')
% legend('S11', 'S22', 'NorthWest')
% grid on

% ======================= Polt Unwraped Phase

figure('name','Phase loading')
subplot(1,2,1);
yyaxis left
plot(w_low,abs(Extralfa{1}),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa{1},'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S11')
grid on
subplot(1,2,2);
yyaxis left
plot(w_low,abs(Extralfa{2}),'Linewidth',1.5);
ylim([0,1.6]);
yyaxis right
plot(w_low,unwrap_Palfa{2},'Linewidth',1.5);
legend('alfa in mag', 'alfa in phase', 'NorthWest')
title('Alfa in S22')
grid on

end