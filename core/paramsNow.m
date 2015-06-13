function [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,k,Ns)
if((k+Ns-1)>size(Umax,2))
    k = k - (size(Umax,2)-(k+Ns-1));
end
umax = Umax(:,k:k+Ns-1);
umin = Umin(:,k:k+Ns-1);
if(isempty(Xmax) || isempty(Xmin))
    xmax = [];
    xmin = [];
else
    xmax = Xmax(:,k+1:k+Ns+1);
    xmin = Xmin(:,k+1:k+Ns+1);
end
end

