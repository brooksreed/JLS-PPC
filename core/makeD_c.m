function D_c = makeD(Pi_c,alpha_c,Nu,Np)
% Pi_c and alpha_c are both (Nv x 1)
dt = diag(Pi_c.*alpha_c);
D_c = kron(dt,kron(eye(Np),eye(Nu)));
end
