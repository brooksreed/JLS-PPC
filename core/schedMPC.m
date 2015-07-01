function [U,cost,status,X_Out,violate_slack] = schedMPC(x_In,...
            b_hat_MPC,p_i,MPC_horizon,A,Bu,M,E,Q,Qf,R,...
            umax,umin,xmin,xmax,u_deadband)
% solve deterministic scheduled MPC with cvx 
% [U,cost,status,X_Out,violate_slack] = schedMPC(x_In,...
%             b_hat_MPC,p_i,MPC_horizon,A,Bu,M,E,Q,Qf,R,...
%             umax,umin,xmin,xmax,u_deadband)
% x_In(N_STATES x 1):  x_hat{t+TAU_C|t-TAU_M},
% b_hat_MPC(N_VEH*N_CONTROLS*horizon x 1): b_hat{t+TAU_C-1}
% p_i (Nv x 1) vector of #steps to constrain control priors 
% horizon: scalar
% A,Bu: system
% M, E: buffer shift
% Q, Qf, R: cost function params (states, terminal states, controls)
% umax, umin: (Nu x horizon) control constraints
% xmin and xmax are optional, default is none
% state constraints implemented with slack variable "barrier"
% uDB is deadband width (NOTE - this makes MPC very slow, requires MIQP
%   solver such as Gurobi)
%
% U is computed plan
% cost, status from CVX
% X_Out output gives predicted state trajectory
% violate_slack is for state constraints (if forced to be violated)
%   gives i (state), j (time step) and violation size

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
        umin = repmat(umin,[1,MPC_horizon]);
    end
    if(size(umax,2)==1)
        umax = repmat(umax,[1,MPC_horizon]);
    end
end

if(isempty(u_deadband))
    use_deadband = 0;
else
    use_deadband = 1;
    if(size(u_deadband,2)==1)
        u_deadband = repmat(u_deadband,[1,MPC_horizon]);
    end
end

if(isempty(xmin) || isempty(xmax))
    state_constraints = 0;
    violate_slack = [];
else
    state_constraints = 1;
end

N_U = size(Bu,2); 
N_x = size(A,1);
N_v = length(p_i);
N_u = N_U/N_v;

blockQ = kron(eye(MPC_horizon),Q);
blockQ((MPC_horizon*N_x-N_x+1):MPC_horizon*N_x,...
    (MPC_horizon*N_x-N_x+1):MPC_horizon*N_x) = Qf;
blockR = kron(eye(MPC_horizon),R);

cvx_clear
cvx_begin 
try
    % first choice not installed with CVX
    cvx_solver gurobi
catch 
    cvx_solver sedumi
end
cvx_quiet(true)

if(state_constraints)
    variable X(Nx,horizon+1) 
    variable U(NU,horizon) 
    variable violate_slack1(Nx,horizon+1) 
    variable violate_slack2(Nx,horizon+1) 
    variable d1(NU,horizon) binary
    variable d2(NU,horizon) binary
else
    variable X(Nx,horizon+1) 
    variable U(NU,horizon) 
    variable d1(NU,horizon) binary
    variable d2(NU,horizon) binary
end

% control constraints


if(use_deadband)
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
for i = 1:N_v
    if(p_i(i)>=1)
        Ei = E((N_u*i-(N_u-1)):(N_u*i),:);
        for k = 1:p_i(i)
            U_prior = Ei*M^(k)*b_hat_MPC;
            % saturate (in case of rounding error - stay feasible)
            if(U_prior>=umax(i,k))
                U_prior=umax(i,k);
            end
            if(U_prior<=umin(i,k))
                U_prior=umin(i,k);
            end
            U((N_u*i-(N_u-1)):(N_u*i),k) == U_prior;
        end
    end
end

% state constraints
if(state_constraints)
        
    violate_slack1 >= 0;
    violate_slack2 >= 0;
    X - violate_slack1 <= xmax; X + violate_slack2 >= xmin;
    
end

X(:,2:MPC_horizon+1) == A*X(:,1:MPC_horizon)+Bu*U;
X(:,1) == x_In; 

x = reshape(X(:,2:end),MPC_horizon*N_x,1);
u = reshape(U,MPC_horizon*N_U,1);
if(state_constraints)        
    minimize( x'*blockQ*x + u'*blockR*u + ...
        1e8*sum(sum(violateSlack1 + violateSlack2)) )
else
    minimize ( x'*blockQ*x + u'*blockR*u )
end

cvx_end
status = cvx_status;

% compute "predicted" cost 
jd = 0;
for kk = 2:MPC_horizon
    jd = jd + X(:,kk)'*Q*X(:,kk) + U(:,kk-1)'*R*U(:,kk-1);
end
jd = jd + X(:,MPC_horizon+1)'*Qf*X(:,MPC_horizon+1);
cost = jd;

% plan doesn't include 1st step (xIn)
X_Out = X(:,2:end);

if(state_constraints)
    % constraint violations don't either
    violate_slack = zeros(N_x,MPC_horizon,2);
    violate_slack(:,:,1) = violate_slack1(:,2:end);
    violate_slack(:,:,2) = violate_slack2(:,2:end);
end
