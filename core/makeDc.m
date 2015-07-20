function Dc_out = makeDc(PI_C,alpha_c_now,N_CONTROLS,N_HORIZON)
% PI_C and alpha_c_now are both (Nv x 1)

% (future speedups:)
% faster kron w/ sparse inputs: 
% http://www.mathworks.com/matlabcentral/fileexchange/28889-kronecker-product

dt = diag(PI_C.*alpha_c_now);
Dc = kron(dt,kron(sparse(eye(N_HORIZON)),sparse(eye(N_CONTROLS))));
Dc_out = sparse(Dc);

end
