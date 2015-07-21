function M_out = makeM(N_CONTROLS,N_HORIZON,N_VEH)
% M is (Nv*Np*Nu) square

base_M = zeros(N_HORIZON);
base_M(1:N_HORIZON-1,2:N_HORIZON) = diag(ones(1,N_HORIZON-1));
base_M(N_HORIZON,N_HORIZON) = 0;   %m_f
m = kron(sparse(base_M),sparse(eye(N_CONTROLS)));
M = kron(eye(N_VEH),sparse(m));
M_out = sparse(M);

end