function D_c = makeD_c(Pi_c,alpha_c,Nu,Np)
% Pi_c and alpha_c are both (Nv x 1)

% v1.0 6/13/2015
dt = diag(Pi_c.*alpha_c);
D_c = kron(dt,kron(eye(Np),eye(Nu)));
end
