function [XFwd,p_i] = prepMPC(t_global,Xh_t,U_in_t,D_cIn_t,Pi_c_all,...
    A,Bu,E,M,Nx,Nv,Nu,Np,tm,tc)
% propagate estimate fwd open-loop tm+tc steps using MJLS system
%
% Xh is xHat_{t-tm|t-tm}, bHat_{(t-1)-tm}
% Uin is U_{t-tm} to U_{t+tc-1}
% D_cIn is Dh_{t-tm} to Dh_{t+tc-1}
% 
% XFwd is the forward-propagated estimate (tm+tc into future)
% p_i is vector of indices to tell schedMPC when to constrain controls to
%   priors (due to scheduling)

% v1.0 6/13/2015

nFwds = tm+tc;
XFwd = zeros(Nv*Np*Nu+Nx,nFwds+1);
XFwd(:,1) = Xh_t;

% moving fwd
for i = 1:nFwds
    
    UFwd = U_in_t(:,i);
    D_cFwd = D_cIn_t(:,:,i);
    
    AAFwd = [A,Bu*E*M*(eye(Np*Nu*Nv)-D_cFwd);...
        zeros(Nv*Np*Nu,Nx),(eye(Np*Nu*Nv)-D_cFwd)*M];
    BUFwd = [Bu*E*D_cFwd;D_cFwd];
    XFwd(:,i+1) = AAFwd*XFwd(:,i)+BUFwd*UFwd;
    
end

%%%%%%%%%%% compute k_p^i's
p_i = zeros(Nv,1);
for i = 1:Nv
    iInds = find(Pi_c_all(i,(t_global):end)==1);
    p_i(i) = iInds(1)-1;
end

end
