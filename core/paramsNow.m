function [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t,Ns)
if((t+Ns-1)>size(Umax,2))
    t = t - (size(Umax,2)-(t+Ns-1));
end
umax = Umax(:,t:t+Ns-1);
umin = Umin(:,t:t+Ns-1);
if(isempty(Xmax) || isempty(Xmin))
    xmax = [];
    xmin = [];
else
    xmax = Xmax(:,t+1:t+Ns+1);
    xmin = Xmin(:,t+1:t+Ns+1);
end
end

