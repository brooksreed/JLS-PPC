function [XFwd,kpi] = prepMPC(t,Xh,Uin,Din,pi,...
    A,Bu,E1,M,Nx,Nv,Nu,Np,tm,tc)
% propagate estimate fwd open-loop tm+tc steps using MJLS system

% (assume no packet-loss for now...)

% Xh is xHat_{t-tm|t-tm}, bHat_{(t-1)-tm}
% Uin is U_{t-tm} to U_{t+tc-1}
% Din is Dh_{t-tm} to Dh_{t+tc-1}

nFwds = tm+tc;
XFwd = zeros(Nv*Np*Nu+Nx,nFwds+1);
XFwd(:,1) = Xh;

% moving fwd
for i = 1:nFwds
    
    UFwd = Uin(:,i);
    DFwd = Din(:,:,i);
    
    AAFwd = [A,Bu*E1*M*(eye(Np*Nu*Nv)-DFwd);...
        zeros(Nv*Np*Nu,Nx),(eye(Np*Nu*Nv)-DFwd)*M];
    BUFwd = [Bu*E1*DFwd;DFwd];
    XFwd(:,i+1) = AAFwd*XFwd(:,i)+BUFwd*UFwd;
    
end

%%%%%%%%%%% compute k_p^i's
kpi = zeros(Nv,1);
for i = 1:Nv
    iInds = find(pi(i,(t):end)==1);
    kpi(i) = iInds(1)-1;
end

end
