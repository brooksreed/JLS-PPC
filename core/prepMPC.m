function [X_fwd,p_i] = prepMPC(t_sim,Xh_in,U_in_t,Dc_in,PI_C,...
    A,Bu,E,M,N_STATES,N_VEH,N_CONTROLS,N_HORIZON,TAU_M,TAU_C)
% propagate estimate fwd open-loop TAU_M + TAU_C steps using JLS system
% [X_fwd,p_i] = prepMPC(t_sim,Xh_in,U_in_t,Dc_in,PI_C,...
%    A,Bu,E,M,N_STATES,N_VEH,N_CONTROLS,N_p,TAU_M,TAU_C)
%   
% Xh_in is x_hat{t-TAU_M|t-TAU_M}, b_hat{(t-1)-TAU_M}
% U_in is U{t-TAU_M:t+TAU_C-1}
% Dc_in is Dc_hat{t-TAU_M:t+TAU_C-1}
% 
% X_fwd is the forward-propagated estimate (TAU_M + TAU_C into future)
% p_i is vector of indices to tell schedMPC when to constrain controls to
%   priors (due to scheduling)

N_FWDS = TAU_M+TAU_C;
X_fwd = zeros(N_VEH*N_HORIZON*N_CONTROLS+N_STATES,N_FWDS+1);
X_fwd(:,1) = Xh_in;

% moving fwd
for i = 1:N_FWDS
    
    U_fwd = U_in_t(:,i);
    Dc_fwd = Dc_in(:,:,i);
    
    AA_fwd = [A,Bu*E*M*(eye(N_HORIZON*N_CONTROLS*N_VEH)-Dc_fwd);...
        zeros(N_VEH*N_HORIZON*N_CONTROLS,N_STATES),...
        (eye(N_HORIZON*N_CONTROLS*N_VEH)-Dc_fwd)*M];
    BU_fwd = [Bu*E*Dc_fwd;Dc_fwd];
    X_fwd(:,i+1) = AA_fwd*X_fwd(:,i)+BU_fwd*U_fwd;
    
end

%%%%%%%%%%% compute k_p^i's
p_i = zeros(N_VEH,1);
for i = 1:N_VEH
    i_inds = find(PI_C(i,(t_sim):end)==1);
    p_i(i) = i_inds(1)-1;
end

end
