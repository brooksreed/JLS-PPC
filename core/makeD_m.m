function D_m = makeS(Pi_m,alpha_m,Ny)
% Pi_m and alpha_m are both (Nv x 1)
% more robust -- design S based on mt and C?
mt = diag(Pi_m.*alpha_m);
D_m = kron(eye(Ny),mt);
end