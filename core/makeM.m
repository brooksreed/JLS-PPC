function M = makeM(Nu,Np,Nv)
% M is (Nv*Np*Nu) square
baseM = zeros(Np);
baseM(1:Np-1,2:Np) = diag(ones(1,Np-1));
baseM(Np,Np) = 0;   %m_f
m = kron(baseM,eye(Nu));
M = kron(eye(Nv),m);
end