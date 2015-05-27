function [U,cost,status,violateSlack,XOut] = detMPC(xIn,...
    Np,A,B,Q,Qf,R,umax,umin,xmin,xmax,Xref)        
% solve deterministic MPC with cvx 
% [U,cost,status,violateSlack,XOut] = detMPC(xIn,...
%     T,A,B,Q,Qf,R,umax,umin,xmin,xmax,Xref)  
% T is horizon
% umax, umin: (Nu x T) control constraints
%       if (Nu x 1), converted via repmat
%       Or, set to [] 
%  xmin and xmax are optional, default is none
%  Xref is optional, default is zeros (regulator)
%  state constraints implemented with slack variable "barrier"
%  violations ouput gives i (state), j (time step) and violation size
%  X output gives predicted state trajectory

% BR, 2/17/2014

% 3/25/2014: vectorized U and X constraints
% added checker and specific display of violated constraints
% 3/26/2014: changed to reshaped vectors of x, u and quadratic cost 
%       (seemed like 'fro' norm of Qhalf, Rhalf caused inaccuracy)
% 8/12/2014: added Xref (optional) input

if(nargin<11)
    Xref = [];
end

if(nargin<10)
    xmin = [];
    xmax = [];
end
if(nargin<8)
    umin = [];
    umax = [];
end


if(isempty(xmin) || isempty(xmax))
    stateConstraints = 0;
    violateSlack = [];
else
    stateConstraints = 1;
end

Nu = size(B,2); Nx = size(A,1);

if(isempty(Xref))
   Xref = zeros(Nx,Np+1);
   
end

if(isempty(umax) || isempty(umin) )
    controlConstraints=0;
else
    controlConstraints = 1;
    if(size(umin,1)==1)
        umin = umin*ones(Nu,1);
    end
    if(size(umax,1)==1)
        umax = umax*ones(Nu,1);
    end
    
    if(size(umin,2)==1)
        umin = repmat(umin,[1,Np]);
    end
    if(size(umax,2)==1)
        umax = repmat(umax,[1,Np]);
    end
end

blockQ = kron(eye(Np),Q);
blockQ((Np*Nx-Nx+1):Np*Nx,(Np*Nx-Nx+1):Np*Nx) = Qf;
blockR = kron(eye(Np),R);

cvx_begin 
cvx_solver gurobi
%cvx_solver sedumi
cvx_quiet(true)

if(stateConstraints)
    variables X(Nx,Np+1) U(Nu,Np) violateSlack1(Nx,Np+1) violateSlack2(Nx,Np+1)
else
    variables X(Nx,Np+1) U(Nu,Np) 
end

% control constraints
if(controlConstraints)
    U <= umax; U >= umin;
end

% state constraints
if(stateConstraints)
    violateSlack1 >= 0;
    violateSlack2 >= 0;
    X - violateSlack1 <= xmax; X + violateSlack2 >= xmin;
end

X(:,2:Np+1) == A*X(:,1:Np)+B*U;
X(:,1) == xIn; 

x = reshape(X(:,2:end),Np*Nx,1);
Xr = reshape(Xref(:,2:end),Np*Nx,1);
u = reshape(U,Np*Nu,1);
if(stateConstraints)        
    minimize( (x-Xr)'*blockQ*(x-Xr) + u'*blockR*u + ...
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

