function Dm = makeDm(PI_M,alpha_m,N_Y_ALL)
% Pi_m and alpha_m are both (Nv x 1)

mt = diag(PI_M.*alpha_m);
Dm = kron(eye(N_Y_ALL),mt);
end