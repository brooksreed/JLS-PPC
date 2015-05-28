function [results] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tap,...
    alphaBar,Pi,Xi,Lambda,umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,...
    w,v,alpha,beta,gamma,nACKHistory)
% runs simulation of MJLS/scheduled PPC

% currently restricts to time-invariant constraint input
% (code in the sim converts to fcn of time)

% BR, 4/23/2014
% modifying for delayed ACKs, 6/13/2014

printouts = 0;

% currently always uses alphaBar state prior adjustment
covPriorAdj = 0;
if(covPriorAdj)
    disp('COV PRIOR ADJUST ON')
end

% initialization stuff
Nx = size(A,1);
Nv = length(alphaBar);
NU = size(Bu,2);    % all controls
Nu = NU/Nv;         % per-vehicle
Nw = size(Bw,2);
NY = size(C,1);
Ny = NY/Nv;

% initially no steps since dropped ACK
tNoACK = zeros(Nv,1);

Umax = repmat(umax,1,Ns+Np);
Umin = repmat(umin,1,Ns+Np);

if(~ischar(codebook))
    % Nearest-Neighbor quantizer: make partition
    partition = ( codebook(2:end)+codebook(1:(end-1)) )/2;
end

M = makeM(Nu,Np,Nv);
e1t = [1,zeros(1,Np-1)];
E1 = kron(eye(Nv),kron(e1t,eye(Nu)));
BW = [Bw;zeros(Nu*Np*Nv,Nw)];

X = zeros(Nx+Nu*Np*Nv,Ns);  % includes x and b
y = zeros(NY,Ns);           % true y (could be reconstructed Cx+v)
yh = zeros(NY,Ns);          % y into estimator (could be reconstructed)
Xh = NaN*zeros(Nx+Nu*Np*Nv,Ns); % includes xHat and bHat
P = zeros(Nx,Nx,Ns);
u = zeros(NU,Ns);
U = zeros(Nu*Np*Nv,Ns);

% initialize Dhat and Dyes
Dh = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns+tc);   % D hat for estimator
Ds = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns+tc);   % used to track utilde for cov. prior
alphat = repmat(alphaBar,[1 Ns]);
for t = (tc+1):Ns
    Dh(:,:,t) = makeD(Pi(:,t-tc),alphat(:,t-tc),Nu,Np);
    Ds(:,:,t) = makeD(Pi(:,t-tc),Pi(:,t-tc),Nu,Np);
end

D = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns);    % true D
S = zeros(Nv*Ny,Nv*Ny,Ns);
a = zeros(Nv,Nv,Ns);

utilde = zeros(NU,Ns);
bYes = zeros(Nu*Np*Nv,1);

XhMPC = NaN*zeros(Nv*Np*Nu+Nx,tm+tc+1,Ns);	% includes xHatMPC and bHatMPC
Jcomp = NaN*zeros(1,Ns);
XPlan = NaN*zeros(Nx,Np,Ns);
MPCtime = NaN*zeros(Ns,1);
looptime = zeros(Ns,1);
MPCFail = zeros(Ns,1);

% initial state x_1 and z_1
X(1:Nx,1) = xIC;
y(:,1) = C*xIC+v(:,1);

% initial buffer - all zeros
X(Nx+1:end,1) = zeros(Nu*Np*Nv,1);

% first step propagation - gives x_2
X(1:Nx,2) =  A*xIC + w(:,1);
y(:,2) = C*X(1:Nx,1) + v(:,2);
u(:,1) = E1*X(Nx+1:end,1);

% Loop starts at "estimator" at time t=tm+1, when first meas (t=1) is RX'd
for t = (tm+1):(Ns-1)
    looptic = tic;
    
    % determine which measurements are available at estimator at this step
    % at step t, meas. are sent at t-tm
    S(:,:,t-tm) = makeS(Xi(:,t-tm),beta(:,t-tm),Ny);
    yh(:,t-tm) = S(:,:,t-tm)*y(:,t-tm);     % available tm steps after sent

    % determine ACKs available at this step
    % update Dh (and alphat), and KFstart
    [Dh,alphat,a,KFstart,tNoACK] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
        Lambda,gamma,t,tm,tc,ta,tap,Nu,Np,tNoACK,nACKHistory);    
    
    %%%%%%%
    % run estimator
    %%%%%%%
    
    % first compute utilde (control if packets all successful)
    if( (t-tm-1) >= 1)
        bYes = M*(eye(Np*Nu*Nv)-Ds(:,:,t-tm-1))*bYes + ...
            Ds(:,:,t-tm-1)*U(:,t-tm-1);
        utilde(:,t-tm-1) = E1*bYes;
    end
    
    % run KF from KFstart up until time of recent measurement
    for td = KFstart:(t-tm)
        if(td<=1)
            AKF = eye(size(A));
            DKFh = makeD(zeros(Nv,1),zeros(Nv,1),Nu,Np);
            XhIn = [xHat1;zeros(Nu*Np*Nv,1)];
            Pin = P1;
            Uin = zeros(Nu*Np*Nv,1);
            ut = zeros(NU,1);
        else
            AKF = A;
            DKFh = Dh(:,:,td-1);
            XhIn = Xh(:,td-1);
            Pin = P(:,:,td-1);
            Uin = U(:,td-1);
            ut = utilde(:,td-1);
        end
        
        % Xh(:,t-tm): xHat_{t-tm|t-tm},bHat_{t-tm-1}
        [Xh(:,td),P(:,:,td)] = JLSKF(XhIn,Pin,yh(:,td),Uin,DKFh,...
            Nx,Nv,Nu,Np,S(:,:,td),AKF,Bu,E1,M,C,W,V,...
            ut,alphaBar,covPriorAdj);
        
    end
    
    %%%%%%%
    % if control is to be computed/sent this step
    if(max(Pi(:,t)))
        
        %%%%%%%
        
        MPCtic = tic;
        
        % grab constraints (to accommodate time-varying)
        [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t+1,Np);
        
        % Forward propagation: XhMPC and k_p^i's
        % starts with Xh: xHat_{t-tm|t-tm}, bHat_{t-tm-1}
        % compute xHatMPC_{t-tm+1:t+tc|t-tm}, bHatMPC_{t-tm:t+tc-1}
        % first step: uses uHat_{t-tm} <-- Dh(:,:,t-tm)
        % *NOTE* Goal of XhMPC(:,end,:) is to match true X(:,:) 
        
        Ufwd = U(:,(t-tm):(t+tc-1));
        Dfwd = Dh(:,:,(t-tm):(t+tc-1));
        [XhMPC(:,:,t+tc),kpi] = prepMPC(t,Xh(:,t-tm),Ufwd,Dfwd,Pi,...
            A,Bu,E1,M,Nx,Nv,Nu,Np,tm,tc);
        
        solveStatus = 0;
        counter = 1;
        TMPC = Np;
        while(solveStatus==0)
            
            % compute U_{t+tc}^i, forall i s.t. {Pi(i,t) = 1}
            [Umpc,Jcomp(t),status,XP,~] = schedMPC(XhMPC(1:Nx,end,t+tc),...
                XhMPC((Nx+1):end,end,t+tc),kpi,TMPC,A,Bu,M,E1,Q,Qf,R,...
                umax,umin,xmin,xmax,[]);
            
            if(strfind(status,'Solved'))
                solveStatus=1;
                fprintf('\nStep %d: %s\n',t,status)
                
            elseif( strcmp(status,'Failed') )
                disp('FAILED')
                
            elseif( strcmp(status,'Infeasible') )
                disp('INFEASIBLE')
                disp(counter)
                disp(XhMPC(:,end,t+tc))
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
        Umpc = Umpc(:,1:Np);    % truncate if TMPC>Np
        XPlan(:,:,t) = XP(:,1:Np);
        
        MPCtime(t) = toc(MPCtic);
        
    else % (for saving - make clear not set)
        
        Umpc = zeros(NU,Np);
        
    end
    
    % translate (NU x Np Umpc) into buffer shape Nu x Np x Nv
    Uvec = reshape(Umpc',[Nu*Np*Nv,1]);
    if(isempty(strfind(codebook,'none')))
        [~,U(:,t+tc)] = quantiz(Uvec,partition,codebook);
    elseif(strfind(codebook,'none'))
        U(:,t+tc) = Uvec;
    end
    
    % true system propagation
    % reshape into MJLS form (for step t), uses D_t(pi(t-tc),alpha(t-tc))
    if(t<=tc)
        D(:,:,t) = makeD(zeros(1,Nv),zeros(1,Nv),Nu,Np);
    else
        D(:,:,t) = makeD(Pi(:,t-tc),alpha(:,t-tc),Nu,Np);
    end
    I = eye(size(D(:,:,1)));
    AA = [A,Bu*E1*M*(I-D(:,:,t));zeros(Nv*Np*Nu,Nx),M*(I-D(:,:,t))];
    BU = [Bu*E1*D(:,:,t);D(:,:,t)];
    
    % propagate system x_{t+1} = f(D_t,x_t,U_t,w_t)
    X(:,t+1) = AA*X(:,t) + BU*U(:,t) + BW*w(:,t);
    
    % actual applied u_t for saving
    u(:,t) = E1*X(Nx+1:end,t+1);
    
    % measurement z_{t+1}:
    y(:,t+1) = C*X(1:Nx,t+1) + v(:,t+1);
    
    looptime(t) = toc(looptic);
    
    %%%%%%%
    if( (t>3) && printouts)
        try
            printoutsJLSPPC
        catch
            disp('printout error')
        end
    end
    
end

xF = X(1:Nx,end);
% compute "actual" cost
jj = u(:,1)'*R*u(:,1);  % control at step 1
for t = 2:Ns
    % states at step 2 (affected by u(1)), through Ns
    jj = jj + X(1:Nx,t)'*Q*X(1:Nx,t) + u(:,t)'*R*u(:,t);
end
% final state (affected by final u)
Jsim = jj + xF'*Q*xF;

results.X = X;
results.u = u;
results.utilde = utilde;
results.Xh = Xh;
results.Jsim = Jsim/Ns;
results.U = U;
results.P = P;
results.XhMPC = XhMPC;
results.Jcomp = Jcomp;
results.XPlan = XPlan;
results.MPCtime = MPCtime;
results.looptime = looptime;
results.MPCFail = MPCFail;
results.rmsEstError = nanrms(X(1:Nv,:) - Xh(1:Nv,:),2);
results.rmsPosError = nanrms(X(1,:),2);

end





