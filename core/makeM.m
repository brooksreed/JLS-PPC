function M = makeM(N_CONTROLS,N_p,N_VEH)
% M is (Nv*Np*Nu) square

base_M = zeros(N_p);
base_M(1:N_p-1,2:N_p) = diag(ones(1,N_p-1));
base_M(N_p,N_p) = 0;   %m_f
m = kron(base_M,eye(N_CONTROLS));
M = kron(eye(N_VEH),m);

end