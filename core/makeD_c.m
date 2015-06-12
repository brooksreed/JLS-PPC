function D = makeD(pi,alpha,Nu,Np)
% pi and alpha are both (Nv x 1)
dt = diag(pi.*alpha);
D = kron(dt,kron(eye(Np),eye(Nu)));
end
