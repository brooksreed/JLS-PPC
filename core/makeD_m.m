function S = makeS(xi,beta,Nz)
% xi and beta are both (Nv x 1)
% NOTE - swapped from original kron(st,eye(Nz)) 
% this matches OP convention with separate Ap, Aq
% more robust -- design S based on st and C
st = diag(xi.*beta);
S = kron(eye(Nz),st);
end