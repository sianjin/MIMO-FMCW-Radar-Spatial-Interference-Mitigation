clear
%% ROC validation and compare
Nsim = 1e6; % number of simulation
%% Antenna setup
M = 4; % # of Tx antennas
N = 4; % # of Rx antennas
d_r = 0.5; % normalized inter-RX antenna distance
d_t = 0.5*N; % normalized inter-TX antenna distance
%% SNR setup
% Reference: https://www.mathworks.com/help/phased/ug/signal-detection-in-white-gaussian-noise.html
b = 1; % normalized signal power
snrdB = -5;
snr = 10^(snrdB/10);
sigma_square = abs(b)^2/snr;
R_noise = sigma_square * eye(M*N); % convariance of noise
L_noise = chol(R_noise/2);
z = L_noise.'*(randn(M*N,Nsim)+1i*randn(M*N,Nsim));
%% Target steering vectors
phi_tar = 30; % target angle
f_tar_r = d_r*sind(phi_tar); % normalized target spatial frequency at RX
f_tar_t = d_t*sind(phi_tar); % normalized target spatial frequency at TX
a_tar_r = exp(-1j*2*pi*(0:N-1)*f_tar_r).'; % target steering vector at RX
a_tar_t = exp(-1j*2*pi*(0:M-1)*f_tar_t).'; % target steering vector at TX
a_tar_v = kron(a_tar_t,a_tar_r); % virtual target steering vector
%% Interference 1 RX steering vector
phi_int1 = 40; % receive interference 1 angle
f_int_r1 = d_r*sind(phi_int1); % normalized interference 1 spatial frequency at RX
a_int_r1 = exp(-1j*2*pi*(0:N-1)*f_int_r1).'; % interference 1 steering vector at RX
a_mix_v1 = kron(a_tar_t,a_int_r1); % virtual essential interference steering vector 1
%% Interference 2 RX steering vector
phi_int2 = 10; % receive interference 2 angle
f_int_r2 = d_r*sind(phi_int2); % normalized interference 2 spatial frequency at RX
a_int_r2 = exp(-1j*2*pi*(0:N-1)*f_int_r2).'; % interference 2 steering vector at RX
a_mix_v2 = kron(a_tar_t,a_int_r2); % virtual essential interference steering vector 2
%% Interference 1 decoded TX steering vector
inrdB1 = -10;
inr1 = 10^(inrdB1/10);
sigma_int_square1 = sigma_square*inr1;
correlationCoeff1 = 0.6;
R_int_normalized1 = zeros(M,M);
for m_row = 1:M
    for m_col = 1:M
        R_int_normalized1(m_row,m_col) = correlationCoeff1^(abs(m_row-m_col));
    end
end
R_int_t1 = sigma_int_square1*R_int_normalized1; % decoded Tx interference 1 covariance matrix
L_int1 = chol(R_int_t1/2);
a_int_t1 = L_int1'*(randn(M,Nsim)+1i*randn(M,Nsim)); % decoded interference steering vector 1 at TX
a_int_v1 = kron(a_int_t1,a_int_r1); % virtual interference steering vector 1 
b_tuta_int_1 = a_tar_t'* a_int_t1/(a_tar_t'*a_tar_t); % essential interference 1 amplitude
h_square1 = (a_tar_t'*R_int_t1*a_tar_t)/norm(a_tar_t)^4; % essential interference 1 power
%% Interference 2 decoded TX steering vector
inrdB2 = -10;
inr2 = 10^(inrdB2/10);
sigma_int_square2 = sigma_square*inr2;
correlationCoeff2 = 0.5;
R_int_normalized2 = zeros(M,M);
for m_row = 1:M
    for m_col = 1:M
        R_int_normalized2(m_row,m_col) = correlationCoeff2^(abs(m_row-m_col));
    end
end
R_int_t2 = sigma_int_square2*R_int_normalized2; % decoded Tx interference 2 covariance matrix
L_int2 = chol(R_int_t2/2); 
a_int_t2 = L_int2'*(randn(M,Nsim)+1i*randn(M,Nsim)); % decoded interference steering vector 2 at TX
a_int_v2 = kron(a_int_t2,a_int_r2); % virtual interference steering vector 2
b_tuta_int_2 = a_tar_t'* a_int_t2/(a_tar_t'*a_tar_t); % essential interference 2 amplitude
h_square2 = (a_tar_t'*R_int_t2*a_tar_t)/norm(a_tar_t)^4; % essential interference 2 power
%% Sample covariance matrix for adaptive processing 
K = M*N; % number of range-Doppler bins for covariance matrix estimation
Rs = zeros(M*N,M*N,Nsim);
for iter =1:Nsim
    A_int_t1_ref = L_int1'*(randn(M,K)+1i*randn(M,K)); % decoded interference steering vector 1 at TX
    A_int_v1 = kron(A_int_t1_ref,a_int_r1); % virtual interference steering vector 1
    A_int_t2 = L_int2'*(randn(M,K)+1i*randn(M,K)); % decoded interference steering vector 2 at TX
    A_int_v2 = kron(A_int_t2,a_int_r2); % virtual interference steering vector 2
    Yref = A_int_v1 + A_int_v2 + L_noise.'*(randn(M*N,K)+1i*randn(M*N,K));
    Rs(:,:,iter) = Yref*Yref'/K; % sample convariance matrix 
end
%% Reconstructed convariance matrix for AGS
% Capon spatial spectrum
gridSize = 1;
Theta_int = -90:gridSize:90;
rho = 10; % scaling parameter
% Threshold for determine interference region
th = zeros(1,Nsim);
for iter = 1:Nsim
    th(iter) = abs(eigs(Rs(:,:,iter),1,'smallestabs')); 
end
% Convariance matrix reconstruction
Nang = length(Theta_int);
Rr = zeros(M*N,M*N,Nsim);
for angIdx = 1:Nang
    phi_int = Theta_int(angIdx);
    f_int_r = d_r*sind(phi_int); % normalized interference spatial frequency at RX
    a_int_r = exp(-1j*2*pi*(0:N-1)*f_int_r).'; % interference steering vector at RX
    a_mix_v = kron(a_tar_t,a_int_r); % virtual essential interference steering vector
    powAng = abs(1/pagemtimes(pagemrdivide(pagemtimes(a_mix_v',ones(1,1,Nsim)),Rs),a_mix_v));
    for iter = 1:Nsim
        if abs(powAng(1,1,iter)) > th(iter)
            Rr(:,:,iter) = Rr(:,:,iter) + rho*powAng(1,1,iter)*(a_mix_v*a_mix_v');
        end
    end
end
Rr = Rr + pagemtimes(eye(M*N),ones(1,1,Nsim));
%% Monte-Carlo simulation under 2 cases
% Under H_0
y_0 = a_int_v1 + a_int_v2 + z;
% clairvoyant
num_clairvoyant_0 = abs(a_tar_v'*(y_0 - a_mix_v1*b_tuta_int_1 - a_mix_v2*b_tuta_int_2)).^2;
den_clairvoyant = M*N;
T_y_0_clairvoyant = 2/sigma_square*num_clairvoyant_0/den_clairvoyant;
T_y_0_sort_clairvoyant = sort(T_y_0_clairvoyant);
% GS
R = h_square1/sigma_square*(a_mix_v1*a_mix_v1') ...
    + h_square2/sigma_square*(a_mix_v2*a_mix_v2') + eye(M*N); % Covariance matrix of noise + essential interference 
num_GS_0 = abs((R\a_tar_v)'*y_0).^2;
den_GS = abs(a_tar_v'*(R\a_tar_v));
T_y_0_GS = 2/sigma_square*num_GS_0/den_GS;
T_y_0_sort_GS = sort(T_y_0_GS);
% LCMV 
R_LCMV = 1/sigma_square*kron(R_int_t1,a_int_r1*a_int_r1') ...
    + 1/sigma_square*kron(R_int_t2,a_int_r2*a_int_r2') + eye(M*N); % Covariance matrix of noise + interference 
num_LCMV_0 = abs((R_LCMV\a_tar_v)'*y_0).^2;
den_LCMV = abs(a_tar_v'*(R_LCMV\a_tar_v));
T_y_0_LCMV = 2/sigma_square*num_LCMV_0./den_LCMV;
T_y_0_sort_LCMV = sort(T_y_0_LCMV);
% LCMV-SMI
num_LCMV_SMI_0 = zeros(1,Nsim);
den_LCMV_SMI = zeros(1,Nsim);
for iter = 1:Nsim
    num_LCMV_SMI_0(iter) = abs((Rs(:,:,iter)\a_tar_v)'*y_0(:,iter)).^2;
    den_LCMV_SMI(iter) = abs(a_tar_v'*(Rs(:,:,iter)\a_tar_v));
end
T_y_0_LCMV_SMI = 2/sigma_square*num_LCMV_SMI_0./den_LCMV_SMI;
T_y_0_sort_LCMV_SMI = sort(T_y_0_LCMV_SMI);
% IAGS
num_IAGS_0 = zeros(1,Nsim);
den_IAGS = zeros(1,Nsim);
for iter = 1:Nsim
    num_IAGS_0(iter) = abs((Rr(:,:,iter)\a_tar_v)'*y_0(:,iter)).^2;
    den_IAGS(iter) = abs(a_tar_v'*(Rr(:,:,iter)\a_tar_v));
end
T_y_0_IAGS = 2/sigma_square*num_IAGS_0./den_IAGS;
T_y_0_sort_IAGS = sort(T_y_0_IAGS);
% Under H_1
y_1 = kron(b*a_tar_v,ones(1,Nsim)) + a_int_v1 + a_int_v2 + z;
% clairvoyant
num_clairvoyant_1 = abs(a_tar_v'*(y_1 - a_mix_v1*b_tuta_int_1 - a_mix_v2*b_tuta_int_2)).^2;
T_y_1_clairvoyant = 2/sigma_square*num_clairvoyant_1/den_clairvoyant;
% GS
num_GS_1 = abs((R\a_tar_v)'*y_1).^2;
T_y_1_GS = 2/sigma_square*num_GS_1/den_GS;
% LCMV 
num_LCMV_1 = abs((R_LCMV\a_tar_v)'*y_1).^2;
T_y_1_LCMV = 2/sigma_square*num_LCMV_1/den_LCMV;
% LCMV-SMI
num_LCMV_SMI_1 = zeros(1,Nsim);
for iter = 1:Nsim
    num_LCMV_SMI_1(iter) = abs((Rs(:,:,iter)\a_tar_v)'*y_1(:,iter)).^2;
end
T_y_1_LCMV_SMI = 2/sigma_square*num_LCMV_SMI_1./den_LCMV_SMI;
% IAGS
num_IAGS_1 = zeros(1,Nsim);
for iter = 1:Nsim
    num_IAGS_1 (iter) = abs((Rr(:,:,iter)\a_tar_v)'*y_1(:,iter)).^2;
end
T_y_1_IAGS = 2/sigma_square*num_IAGS_1./den_IAGS;
%% Monte-Carlo ROC curve
% false alarm rate
count_fa = [1e2:1e2:1e3,2e3:1e3:1e4,2e4:1e4:1e5,2e5:1e5:1e6];
% count_fa = [1e2:1e2:1e3,2e3:1e3:1e4,2e4:1e4:1e5];
Pfa = count_fa/Nsim;
% clairvoyant
Pd_MC_clairvoyant = zeros(1,length(count_fa));
thresh_clairvoyant = zeros(1,length(count_fa));
% GS
Pd_MC_GS = zeros(1,length(count_fa));
thresh_GS = zeros(1,length(count_fa));
% LCMV
Pd_MC_LCMV = zeros(1,length(count_fa));
thresh_LCMV = zeros(1,length(count_fa));
% LCMV-SMI
Pd_MC_LCMV_SMI = zeros(1,length(count_fa));
thresh_LCMV_SMI = zeros(1,length(count_fa));
% IAGS
Pd_MC_IAGS = zeros(1,length(count_fa));
thresh_IAGS = zeros(1,length(count_fa));
for iter = 1:length(count_fa)
    % clairvoyant
    thresh_iter_clairvoyant = T_y_0_sort_clairvoyant(end-count_fa(iter)+1);
    Pd_MC_clairvoyant(iter) = sum(T_y_1_clairvoyant>thresh_iter_clairvoyant)/Nsim;
    thresh_clairvoyant(iter) = thresh_iter_clairvoyant;
    % GS
    thresh_iter_GS = T_y_0_sort_GS(end-count_fa(iter)+1);
    Pd_MC_GS(iter) = sum(T_y_1_GS>thresh_iter_GS)/Nsim;
    thresh_GS(iter) = thresh_iter_GS;
    % LCMV
    thresh_iter_LCMV = T_y_0_sort_LCMV(end-count_fa(iter)+1);
    Pd_MC_LCMV(iter) = sum(T_y_1_LCMV>thresh_iter_LCMV)/Nsim;
    thresh_LCMV(iter) = thresh_iter_LCMV;
    % LCMV-SMI
    thresh_iter_LCMV_SMI = T_y_0_sort_LCMV_SMI(end-count_fa(iter)+1);
    Pd_MC_LCMV_SMI(iter) = sum(T_y_1_LCMV_SMI>thresh_iter_LCMV_SMI)/Nsim;
    thresh_LCMV_SMI(iter) = thresh_iter_LCMV_SMI;
    % IAGS
    thresh_iter_IAGS = T_y_0_sort_IAGS(end-count_fa(iter)+1);
    Pd_MC_IAGS(iter) = sum(T_y_1_IAGS>thresh_iter_IAGS)/Nsim;
    thresh_IAGS(iter) = thresh_iter_IAGS;
end
figure
plot(Pfa, Pd_MC_clairvoyant,'-.')
hold on
plot(Pfa, Pd_MC_GS,'--*')
hold on
plot(Pfa, Pd_MC_LCMV_SMI,'--v')
hold on
plot(Pfa, Pd_MC_IAGS,':o')
xlabel('P_{FA}')
ylabel('P_{D}')
grid on
set(gca, 'XScale', 'log')
legend('Clairvoyant (Monte-Carlo)','GS (Monte-Carlo)','LCMV-SMI (Monte-Carlo)','AGS (Monte-Carlo)')
%% Theoretical ROC curve
% theoretical detection threshold
gamma = -2*log(Pfa); 
% clairvoyant
lambda_clairvoyant = 2*(abs(b))^2/sigma_square*den_clairvoyant;
Pd_theory_clairvoyant = zeros(1,length(count_fa));
% GS
lambda_GS = 2*(abs(b))^2/sigma_square*den_GS;
Pd_theory_GS = zeros(1,length(count_fa));
% LCMV (ideal convariance matrix)
lambda_LCMV = 2*(abs(b))^2/sigma_square*den_LCMV;
Pd_theory_LCMV = zeros(1,length(count_fa));
for iter = 1:length(count_fa)
    % clairvoyant
    Pd_theory_clairvoyant(iter) = marcumq(sqrt(lambda_clairvoyant),sqrt(gamma(iter))); 
    % GS
    Pd_theory_GS(iter) = marcumq(sqrt(lambda_GS),sqrt(gamma(iter))); 
    % LCMV (ideal convariance matrix)
    Pd_theory_LCMV(iter) = marcumq(sqrt(lambda_LCMV),sqrt(gamma(iter))); 
end
%% Plot ROC curve
figure
plot(Pfa, Pd_theory_clairvoyant,'-.',LineWidth = 2)
hold on
plot(Pfa, Pd_theory_GS,'-*',LineWidth = 1)
hold on
plot(Pfa, Pd_MC_LCMV_SMI,'--v',LineWidth = 1)
hold on
plot(Pfa, Pd_MC_IAGS,':o',LineWidth = 1)
grid on
set(gca, 'XScale', 'log')
xlabel('P_{FA}')
ylabel('P_{D}')
legend('Clairvoyant (theory)','GS (theory)','LCMV-SMI (Monte-Carlo)','AGS (Monte-Carlo)')
