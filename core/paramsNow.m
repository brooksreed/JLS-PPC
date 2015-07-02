function [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t,N_HORIZON)
% constructs "local window" constraint time series for MPC
% inputs are arrays, 2nd dimension is time series

if((t+N_HORIZON-1)>size(Umax,2))
    t = t - (size(Umax,2)-(t+N_HORIZON-1));
end
umax = Umax(:,t:t+N_HORIZON-1);
umin = Umin(:,t:t+N_HORIZON-1);
if(isempty(Xmax) || isempty(Xmin))
    xmax = [];
    xmin = [];
else
    xmax = Xmax(:,t+1:t+N_HORIZON+1);
    xmin = Xmin(:,t+1:t+N_HORIZON+1);
end
end

