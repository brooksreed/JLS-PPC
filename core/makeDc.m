function Dc = makeDc(PI_C,alpha_c_now,N_CONTROLS,N_HORIZON)
% PI_C and alpha_c_now are both (Nv x 1)

dt = diag(PI_C.*alpha_c_now);
Dc = kron(dt,kron(eye(N_HORIZON),eye(N_CONTROLS)));
end
