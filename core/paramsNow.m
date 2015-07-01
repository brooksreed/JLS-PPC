function [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t,T_MPC)
% constructs "local window" constraint time series for MPC
% inputs are arrays, 2nd dimension is time series

if((t+T_MPC-1)>size(Umax,2))
    t = t - (size(Umax,2)-(t+T_MPC-1));
end
umax = Umax(:,t:t+T_MPC-1);
umin = Umin(:,t:t+T_MPC-1);
if(isempty(Xmax) || isempty(Xmin))
    xmax = [];
    xmin = [];
else
    xmax = Xmax(:,t+1:t+T_MPC+1);
    xmin = Xmin(:,t+1:t+T_MPC+1);
end
end

