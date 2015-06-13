function D_m = makeD_m(Pi_m,alpha_m,Ny)
% Pi_m and alpha_m are both (Nv x 1)

% v1.0 6/13/2015
% to do: would it be more robust to shape based on C?
mt = diag(Pi_m.*alpha_m);
D_m = kron(eye(Ny),mt);
end