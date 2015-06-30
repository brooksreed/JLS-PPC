function [U,cost,status,XOut,violateSlack] = schedMPC(xIn,...
            bHatMPC,p_i,horizon,A,Bu,M,E,Q,Qf,R,...
            umax,umin,xmin,xmax,u_deadband)
% solve deterministic scheduled MPC with cvx 
% [U,cost,status,XOut,violateSlack] = schedMPC(xIn,...
%            bHatMPC,kpi,horizon,A,Bu,M,E,Q,Qf,R,umax,umin,xmin,xmax,uDB)
% xIn: (Nx x 1) xHat_{t+tc|t-tm}, bHatMPC: (Nv*Nu*horizon x 1) bHat_{t+tc-1}
% k_priors: Nv x 1 vector of #steps to constrain control priors 
% horizon: horizon; A,Bu: system
% M, E: buffer shift
% Q, Qf, R: LQR params
% umax, umin: (Nu x horizon) control constraints
% xmin and xmax are optional, default is none
% state constraints implemented with slack variable "barrier"
% violations ouput gives i (state), j (time step) and violation size
% uDB is deadband width (NOTE - this makes MPC very slow, requires MIQP
%   solver such as Gurobi)
%
% U is computed plan
% cost, status from CVX
% XOut output gives predicted state trajectory
% violateSlack is for state constraints (if forced to be violated)

% BR, 4/29/2014
%{
- 4/30/2014: added saturation on control priors (rounding error)
- 8/17/2014: added uDB (deadband on control actions) -- symmetric
%}

% v1.0 6/13/2015

% TO DO: 
% try-catch for cvx solver slow on startup?
% [try a faster solver such as QPOasis?]

if(nargin<15)
    u_deadband = [];
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
        umin = repmat(umin,[1,horizon]);
    end
    if(size(umax,2)==1)
        umax = repmat(umax,[1,horizon]);
    end
end

if(isempty(u_deadband))
    uDeadband = 0;
else
    uDeadband = 1;
    if(size(u_deadband,2)==1)
        u_deadband = repmat(u_deadband,[1,horizon]);
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
Nv = length(p_i);
Nu = NU/Nv;

blockQ = kron(eye(horizon),Q);
blockQ((horizon*Nx-Nx+1):horizon*Nx,(horizon*Nx-Nx+1):horizon*Nx) = Qf;
blockR = kron(eye(horizon),R);

cvx_clear
cvx_begin 
try
    % first choice not installed with CVX
    cvx_solver gurobi
catch 
    cvx_solver sedumi
end
cvx_quiet(true)

if(stateConstraints)
    variable X(Nx,horizon+1) 
    variable U(NU,horizon) 
    variable violateSlack1(Nx,horizon+1) 
    variable violateSlack2(Nx,horizon+1) 
    variable d1(NU,horizon) binary
    variable d2(NU,horizon) binary
else
    variable X(Nx,horizon+1) 
    variable U(NU,horizon) 
    variable d1(NU,horizon) binary
    variable d2(NU,horizon) binary
end

% control constraints


if(uDeadband)
    %U >= uDB;
    %U <= -uDB;
    U - (umin - u_deadband).*(1-d2) >= u_deadband;
    U - (umax - u_deadband).*d2 <= u_deadband;
    U - (umin + u_deadband).*d1 >= -u_deadband;
    U - (umax + u_deadband).*(1-d1) <= -u_deadband;
    U >= umin.*d1;
    U <= umax.*d2;
    d1+d2 <= 1;
    
elseif(controlConstraints)
    U <= umax; U >= umin;
end
% add extra constraints equal to priors
for i = 1:Nv
    if(p_i(i)>=1)
        E1i = E((Nu*i-(Nu-1)):(Nu*i),:);
        for k = 1:p_i(i)
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

X(:,2:horizon+1) == A*X(:,1:horizon)+Bu*U;
X(:,1) == xIn; 

x = reshape(X(:,2:end),horizon*Nx,1);
u = reshape(U,horizon*NU,1);
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
for kk = 2:horizon
    jd = jd + X(:,kk)'*Q*X(:,kk) + U(:,kk-1)'*R*U(:,kk-1);
end
jd = jd + X(:,horizon+1)'*Qf*X(:,horizon+1);
cost = jd;

% plan doesn't include 1st step (xIn)
XOut = X(:,2:end);

if(stateConstraints)
    % constraint violations don't either
    violateSlack = zeros(Nx,horizon,2);
    violateSlack(:,:,1) = violateSlack1(:,2:end);
    violateSlack(:,:,2) = violateSlack2(:,2:end);
end
