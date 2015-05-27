function [U,cost,status,XOut,violateSlack] = schedMPC(xIn,...
            bHatMPC,kpi,Np,A,Bu,M,E1,Q,Qf,R,umax,umin,xmin,xmax,uDB)
% solve deterministic scheduled MPC with cvx 
% [U,cost,status,XOut,violateSlack] = schedMPC(xIn,...
%            bHatMPC,kpi,Np,A,Bu,M,E1,Q,Qf,R,umax,umin,xmin,xmax)
% xIn: (Nx x 1) xHat_{t+tc|t-tm}, bHatMPC: (Nv*Nu*Np x 1) bHat_{t+tc-1}
% kpi: Nv x 1 vector of #steps to constrain control priors 
% Np: horizon; A,Bu: system
% M, E1: buffer shift
% Q, Qf, R: LQR params
% umax, umin: (Nu x Np) control constraints
% xmin and xmax are optional, default is none
% state constraints implemented with slack variable "barrier"
% violations ouput gives i (state), j (time step) and violation size
% XOut output gives predicted state trajectory

% BR, 4/29/2014

%{
- 4/30/2014: added saturation on control priors (rounding error)
- 8/17/2014: added uDB (deadband on control actions) -- symmetric
%}

if(nargin<15)
    uDB = [];
end

if(nargin<14)
    xmin = [];
    xmax = [];
end

if(nargin<12)
    umin = [];
    umax = [];
end

if(isempty(umax) || isempty(umin) )
    controlConstraints=0;
else
    controlConstraints = 1;
    if(size(umin,2)==1)
        umin = repmat(umin,[1,Np]);
    end
    if(size(umax,2)==1)
        umax = repmat(umax,[1,Np]);
    end
end

if(isempty(uDB))
    uDeadband = 0;
else
    uDeadband = 1;
    if(size(uDB,2)==1)
        uDB = repmat(uDB,[1,Np]);
    end
end

if(isempty(xmin) || isempty(xmax))
    stateConstraints = 0;
    violateSlack = [];
else
    stateConstraints = 1;
end

NU = size(Bu,2); 
Nx = size(A,1);
Nv = length(kpi);
Nu = NU/Nv;

blockQ = kron(eye(Np),Q);
blockQ((Np*Nx-Nx+1):Np*Nx,(Np*Nx-Nx+1):Np*Nx) = Qf;
blockR = kron(eye(Np),R);

cvx_clear
cvx_begin 
cvx_solver gurobi
%cvx_solver sedumi
cvx_quiet(true)

if(stateConstraints)
    variable X(Nx,Np+1) 
    variable U(NU,Np) 
    variable violateSlack1(Nx,Np+1) 
    variable violateSlack2(Nx,Np+1) 
    variable d1(NU,Np) binary
    variable d2(NU,Np) binary
else
    variable X(Nx,Np+1) 
    variable U(NU,Np) 
    variable d1(NU,Np) binary
    variable d2(NU,Np) binary
end

% control constraints


if(uDeadband)
    %U >= uDB;
    %U <= -uDB;
    U - (umin - uDB).*(1-d2) >= uDB;
    U - (umax - uDB).*d2 <= uDB;
    U - (umin + uDB).*d1 >= -uDB;
    U - (umax + uDB).*(1-d1) <= -uDB;
    U >= umin.*d1;
    U <= umax.*d2;
    d1+d2 <= 1;
    
elseif(controlConstraints)
    U <= umax; U >= umin;
end
% add extra constraints equal to priors
for i = 1:Nv
    if(kpi(i)>=1)
        E1i = E1((Nu*i-(Nu-1)):(Nu*i),:);
        for k = 1:kpi(i)
            UPrior = E1i*M^(k)*bHatMPC;
            % saturate (in case of rounding error - stay feasible)
            if(UPrior>=umax(i,k))
                UPrior=umax(i,k);
            end
            if(UPrior<=umin(i,k))
                UPrior=umin(i,k);
            end
            U((Nu*i-(Nu-1)):(Nu*i),k) == UPrior;
        end
    end
end

% state constraints
if(stateConstraints)
        
    violateSlack1 >= 0;
    violateSlack2 >= 0;
    X - violateSlack1 <= xmax; X + violateSlack2 >= xmin;
    
end

X(:,2:Np+1) == A*X(:,1:Np)+Bu*U;
X(:,1) == xIn; 

x = reshape(X(:,2:end),Np*Nx,1);
u = reshape(U,Np*NU,1);
if(stateConstraints)        
    minimize( x'*blockQ*x + u'*blockR*u + ...
        1e8*sum(sum(violateSlack1 + violateSlack2)) )
else
    minimize ( x'*blockQ*x + u'*blockR*u )
end

cvx_end
status = cvx_status;

% compute "predicted" cost 
jd = 0;
for kk = 2:Np
    jd = jd + X(:,kk)'*Q*X(:,kk) + U(:,kk-1)'*R*U(:,kk-1);
end
jd = jd + X(:,Np+1)'*Qf*X(:,Np+1);
cost = jd;

% plan doesn't include 1st step (xIn)
XOut = X(:,2:end);

if(stateConstraints)
    % constraint violations don't either
    violateSlack = zeros(Nx,Np,2);
    violateSlack(:,:,1) = violateSlack1(:,2:end);
    violateSlack(:,:,2) = violateSlack2(:,2:end);
end
