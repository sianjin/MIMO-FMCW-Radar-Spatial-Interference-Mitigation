function Rr = getCovRecon(Rs,M,N,a_tar_t)
% Capon spatial spectrum
gridSize = 1;
Theta_int = -90:gridSize:90;
rho = 10; % scaling parameter
th = abs(eigs(Rs,1,'smallestabs')); % threshold for determine interference region
% Convariance matrix reconstruction
Nang = length(Theta_int);
Rr = eye(M*N);
for angIdx = 1:Nang
    phi_int = Theta_int(angIdx);
    a_int_r = exp(1j*2*pi*(0:N-1)'/2*sind(phi_int));
    a_mix_v = kron(a_tar_t,a_int_r); % virtual essential interference steering vector
    powAng = abs(1/(a_mix_v'/Rs*a_mix_v));
    if powAng > th
        Rr = Rr + rho*powAng*(a_mix_v*a_mix_v');
    end
end
end