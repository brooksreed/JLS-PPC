function [results] = simNoComms(Ns,Np,Ap,ApModel,Aq,Bu,Bw,C,CModel,...
    Q,Qf,R,W,V,umax,umin,Xmax,Xmin,xi,xIC,P1,xHat1,w,v)
% runs simulation of system with no comms
% possibly using a (deterministic) measurement schedule
% toggle (INSIDE THIS FCN) whether vehicle meas. follow schedule
% intended for Loners+MPC (and maybe coupled LowerBound) sims


% currently restricts to time-invariant constraint input
% (code in the sim converts to fcn of time)

% BR, 8/19/2014


vehicleNavSched = 1; % if vehicle nav follows meas schedule
fprintf('Vehicle nav sched: %d\n',vehicleNavSched)

% initialization stuff

% construct A for propogating true system
Nxq = size(Aq,1);
Nxp = size(Ap,1);
Nxpm = size(ApModel,1);
Nv = size(Aq,1);

A = [Ap,zeros(Nxp,Nxq);zeros(Nxq,Nxp),Aq];
Acontrol = [ApModel,zeros(Nxpm,Nxq);zeros(Nxq,Nxpm),Aq];

Nx = size(A,1);
Nxm = size(Acontrol,1);
NU = size(Bu,2);    % all controls
NZ = size(C,1);
Nz = NZ/Nv;

% model B
Bum = [zeros(Nxpm,NU);eye(NU)];

Umax = repmat(umax,1,Ns+Np);
Umin = repmat(umin,1,Ns+Np);

% construct A for estimation and control

X = zeros(Nx,Ns);
y = zeros(NZ,Ns);           % true z (could be reconstructed Cx+v)
Xh = NaN*zeros(Nxm,Ns);      % includes xHat and bHat
P = zeros(Nxm,Nxm,Ns);
u = zeros(NU,Ns);
S = zeros(NZ,NZ,Ns);

Jcomp = NaN*zeros(1,Ns);
XPlan = NaN*zeros(Nxm,Np,Ns);
MPCtime = NaN*zeros(Ns,1);
looptime = zeros(Ns,1);
MPCFail = zeros(Ns,1);

% initial state x_1 and z_1
X(:,1) = xIC;
y(:,1) = C*xIC+v(:,1);
P(:,:,1) = P1;
Xh(:,1) = xHat1;

% first step propagation - gives x_2
X(:,2) =  A*xIC + w(:,1);
y(:,2) = C*X(:,1) + v(:,2);

% Loop starts at "estimator" at time t=2, when first meas (t=1) is RX'd
% (future: make init such that loop starts at tm+1)
for t = 2:(Ns-1)
    looptic = tic;
      
    %%%%%%%
    % run estimator using MODEL
    %%%%%%%
    
    % meas schedule
    if(vehicleNavSched)
        S(:,:,t) = kron(eye(Nz),diag(xi(:,t)));
    else
        S(:,:,t) = blkdiag(diag(xi(:,t)),eye(Nv));
    end

    % covariance prior (standard)
    Ppre = Acontrol*P(:,:,t-1)*Acontrol' + W ; %P_{t+1|t}

    % Kalman gain = fcn of Ppre
    L = ((Ppre*CModel')/(CModel*Ppre*CModel'+V))*S(:,:,t);   % L_{t}
    P(:,:,t) = Ppre - L*(CModel*Ppre);

    % propagate system
    I = eye(size(Acontrol));
    Xh(:,t) = (I-L*CModel)*(Acontrol*Xh(:,t-1) + Bum*u(:,t-1)) + ...
        L*S(:,:,t)*y(:,t);
        
    %%%%%%%

    MPCtic = tic;

    % grab constraints (to accommodate time-varying)
    [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t+1,Np);

    solveStatus = 0;
    counter = 1;
    TMPC = Np;
    while(solveStatus==0)

        [Umpc,Jcomp(t),status,~,XP] = detMPC(Xh(:,t),TMPC,...
            Acontrol,Bum,Q,Qf,R,umax,umin,xmin,xmax,[]);
        
        if(strfind(status,'Solved'))
            solveStatus=1;
            fprintf('\nStep %d: %s\n',t,status)

        elseif( strcmp(status,'Failed') )
            disp('FAILED')

        elseif( strcmp(status,'Infeasible') )
            disp('INFEASIBLE')
            disp(counter)
            TMPC = TMPC+4;
            [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,...
                t+1,TMPC);
        end

        counter = counter+1;
        if(counter>2)
            disp('MAXCOUNTER')
            Umpc = zeros(NU,Np);
            MPCFail(t) = 1;
            break
        end

    end

    XPlan(:,:,t) = XP(:,1:Np);
    MPCtime(t) = toc(MPCtic);
    u(:,t) = Umpc(:,1);

    % Propogate system using true A
    X(:,t+1) = A*X(:,t) + Bu*u(:,t) + Bw*w(:,t);
    y(:,t+1) = C*X(:,t+1) + v(:,t+1);

    looptime(t) = toc(looptic);
    
    %%%%%%%
    %if( (t>3) && printouts)
    %    printoutsMJLSPPC
    %end
    
end

% xF = X(:,end);
% % compute "actual" cost
% jj = u(:,1)'*R*u(:,1);  % control at step 1
% for t = 2:Ns
%     % states at step 2 (affected by u(1)), through Ns
%     jj = jj + X(:,t)'*Q*X(:,t) + u(:,t)'*R*u(:,t);
% end
% % final state (affected by final u)
% Jsim = jj + xF'*Q*xF;

results.X = X;
results.u = u;
results.Xh = Xh;
%results.Jsim = Jsim/Ns;
results.P = P;
results.Jcomp = Jcomp;
results.XPlan = XPlan;
results.MPCtime = MPCtime;
results.looptime = looptime;
results.MPCFail = MPCFail;

end

