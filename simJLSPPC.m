function [results] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tac,...
    alpha_cBar,Pi_c,Pi_m,Pi_a,Ts,umax,umin,codebook,Xmax,Xmin,xIC,P1,...
    xHat1,w,v,alpha_c,alpha_m,alpha_a,covPriorAdj,nACKHistory)
% runs simulation of MJLS/scheduled PPC

% currently restricts to time-invariant constraint input
% (code in the sim converts to fcn of time)

% BR, 4/23/2014
% modifying for delayed ACKs, 6/13/2014
% v1.0 6/13/2015
% v1.1 6/16/2015

% TO DO:
% Clean up printouts?
% Make sure prep for KF is ready for MIMO and multiple P*s
% Add logging of P*, P**, etc.?
% Add more automated tests/checks?

printDebug = 1;

% currently always uses alphaBar state prior adjustment
if(covPriorAdj)
    disp('COV PRIOR ADJUST ON')
end

% initialization stuff
Nx = size(A,1);
Nv = length(alpha_cBar);
NU = size(Bu,2);    % all controls
Nu = NU/Nv;         % per-vehicle
Nw = size(Bw,2);
NY = size(C,1);
Ny = NY/Nv;

% initialize tNoACK to zeros (no ACK for 'step 0')
tNoACK = zeros(Nv,Ns);

% If initial controls unknown
tNoACK(:,1:Ns) = 1:Ns;

% after first planned control RX, unknown
% initialize based on 2nd pi_c + tau_c ?? (need tau_ac too?)
% tmp = find(pi_c);
% start = tmp(2)+tc
% tNoACK(:,start:Ns) = 1:(Ns-start);



% tNoACK *AS KNOWN EACH STEP*
tNoACKSave = cell(Ns);

Umax = repmat(umax,1,Ns+Np);
Umin = repmat(umin,1,Ns+Np);

if(~ischar(codebook))
    % Nearest-Neighbor quantizer: make partition
    partition = ( codebook(2:end)+codebook(1:(end-1)) )/2;
end

M = makeM(Nu,Np,Nv);
etmp = [1,zeros(1,Np-1)];
E = kron(eye(Nv),kron(etmp,eye(Nu)));
BW = [Bw;zeros(Nu*Np*Nv,Nw)];

X = zeros(Nx+Nu*Np*Nv,Ns);  % includes x and b
y = zeros(NY,Ns);           % true y (could be reconstructed Cx+v)
yh = zeros(NY,Ns);          % y into estimator (could be reconstructed)
Xh = NaN*zeros(Nx+Nu*Np*Nv,Ns); % includes xHat and bHat
P = zeros(Nx,Nx,Ns);
u = zeros(NU,Ns);
U = zeros(Nu*Np*Nv,Ns);

% initialize Dhat and Dyes
D_cHat = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns+tc);   % D hat for estimator
D_cNoLoss = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns+tc);   % used to track utilde for cov. prior
alphaHat = repmat(alpha_cBar,[1 Ns]);
for t = (tc+1):Ns
    D_cHat(:,:,t) = makeD_c(Pi_c(:,t-tc),alphaHat(:,t-tc),Nu,Np);
    D_cNoLoss(:,:,t) = makeD_c(Pi_c(:,t-tc),Pi_c(:,t-tc),Nu,Np);
end

D_c = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns);    % true D
D_m = zeros(Nv*Ny,Nv*Ny,Ns);
D_a = zeros(Nv,Nv,Ns);

uNoLoss = zeros(NU,Ns);
bNoLoss = zeros(Nu*Np*Nv,Ns);

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
u(:,1) = E*X(Nx+1:end,1);

% Loop starts at "estimator" at time t=tm+1, when first meas (t=1) is RX'd
for t = (tm+1):(Ns-1)
    looptic = tic;
    
    % determine which measurements are available at estimator at this step
    % at step t, meas. are sent at t-tm
    D_m(:,:,t-tm) = makeD_m(Pi_m(:,t-tm),alpha_m(:,t-tm),Ny);
    yh(:,t-tm) = D_m(:,:,t-tm)*y(:,t-tm);     % available tm steps after sent
    
    if(printDebug)
        fprintf('\n~~~STEP t=%d AT ESTIMATOR~~~\n',t)
        for i = 1:Nv
            if(Pi_m(i,t-tm)*alpha_m(i,t-tm)==1)
                fprintf('\nt=%d, t-%d Meas %d RX success\n',t,i,tm)
            end
        end
    end
    
    % determine ACKs available at this step
    % update Dh (and alphaHat), KFstart
    % uses tNoACK(t-ta-1) to determine how far back to update
    % updates tNoACK(t-ta) and history based on ACKs RX'd now
    % also increments a lookahead of tNoACK(t-ta+1 --> future)
    [D_cHat,alphaHat,D_a,KFstart,tNoACK,~] = JLSJumpEstimator(D_cHat,...
        Pi_c,D_a,alpha_c,alphaHat,Pi_a,Ts,alpha_a,t,tm,tc,ta,tac,...
        Nu,Np,tNoACK,nACKHistory,printDebug);
    tNoACKSave{t} = tNoACK;
    
    if(printDebug)
        disp(tNoACK)
    end
    
    %%%%%%%%%%%%%%%
    % run estimator
    %%%%%%%%%%%%%%%
    
    % compute buffer and control action if all control packets through
    if( (t-tm-1) >= 1)
        if(t-tm-1>1)
            bPrev = bNoLoss(:,t-tm-2);
        else
            bPrev = zeros(Nu*Np*Nv,1);
        end
        bNoLoss(:,t-tm-1) = M*(eye(Np*Nu*Nv)-D_cNoLoss(:,:,t-tm-1))*bPrev + ...
            D_cNoLoss(:,:,t-tm-1)*U(:,t-tm-1);
        uNoLoss(:,t-tm-1) = E*bNoLoss(:,t-tm-1);
    end
    
    % run KF from KFstart up until time of recent measurement
    for tKF = KFstart:(t-tm)
        
        if(covPriorAdj)
                        
            % determine tNoACK vector for specific filter step
            if(Nv==1)
                if(tKF-1-tac>0)
                    tNoACK_KF = tNoACK(i,tKF+tac-1);
                else
                    tNoACK_KF = 1;
                end
            else
                
%                 tNoACK_KF = zeros(1,Nv);
%                 for i = 1:Nv
%                     if(td-1-tac(i)>0)
%                         tNoACK_KF(i) = tNoACK(i,td+tac(i)-1);
%                     end
%                 end
            
                % ambiguity: ACKs don't ACK all channels...
                % which P*, P**, etc. to use?  
                
                % for now, use Pstar each step? 
                % (only 1-step formulated anyway)
                % (artificially constrain dU to zero for ACK'd channels?)
                tNoACK_KF = 1;
                
            end
            
            % prepare control options for use in cov. prior adj.
            UOptions = zeros(NU,tNoACK_KF);
            for k = 1:tNoACK_KF
                if(tKF-k<1)
                    bTMP = zeros(Nu*Np*Nv,1);
                else
                    bTMP = bNoLoss(:,tKF-k);
                end
                UOptions(:,k) = E*M^k*bTMP;
            end
            
        else
            
            UOptions = [];
            tNoACK_KF = [];
            
        end
        
        if(tKF<=1)
            AKF = eye(size(A));
            DKFh = makeD_c(zeros(Nv,1),zeros(Nv,1),Nu,Np);
            XhIn = [xHat1;zeros(Nu*Np*Nv,1)];
            Pin = P1;
            Uin = zeros(Nu*Np*Nv,1);
            yIn = yh(:,1);
            SIn = D_m(:,:,1);
        else
            AKF = A;
            DKFh = D_cHat(:,:,tKF-1);
            XhIn = Xh(:,tKF-1);
            Pin = P(:,:,tKF-1);
            Uin = U(:,tKF-1);
            yIn = yh(:,tKF);
            SIn = D_m(:,:,tKF);
        end
        
        % Xh(:,t-tm): xHat_{t-tm|t-tm},bHat_{t-tm-1}
        [Xh(:,tKF),P(:,:,tKF)] = JLSKF(XhIn,Pin,yIn,Uin,DKFh,...
            Nx,Nv,Nu,Np,SIn,AKF,Bu,E,M,C,W,V,...
            t,tKF,tNoACK_KF,covPriorAdj,UOptions,alpha_cBar);
        
        if(printDebug)
            if(covPriorAdj)
                if(Nv==1)
                    fprintf('\nt=%d, KF tKF=%d, tNoACK_KF(%d)=%d\n',...
                        t,tKF,tKF+tac-1,tNoACK_KF)
                else
                    %   fprintf('\nt=%d, KF td=%d, tNoACK_KF(%d)=%d\n',...
                    %      t,td,td+tac(i)-1,tNoACK_KF(i))
                    % (add in once resolve ambiguity)
                end
            end
            if(size(A,1)==1)
                % only print estimate for scalar sys
                % (ADD A WINDOW?)
                dispStart = tKF-5;
                if(dispStart<1);dispStart=1;end
                fprintf('Xh(:,%d:%d)=\n',dispStart,tKF);disp(Xh(:,dispStart:tKF)')
            end
        end
        
    end
    
    %%%%%%%
    % if control is to be computed/sent this step
    if(max(Pi_c(:,t)))
        
        %%%%%%%
        
        MPCtic = tic;
        
        % grab constraints (to accommodate time-varying)
        [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t+1,Np);
        
        % Forward propagation: XhMPC and k_p^i's
        % starts with Xh: xHat_{t-tm|t-tm}, bHat_{t-tm-1}
        % compute xHatMPC_{t-tm+1:t+tc|t-tm}, bHatMPC_{t-tm:t+tc-1}
        % first step: uses uHat_{t-tm} <-- Dh(:,:,t-tm)
        % *DEBUG NOTE* Goal of XhMPC(:,end,:) is to match true X(:,:)
        
        Ufwd = U(:,(t-tm):(t+tc-1));
        Dfwd = D_cHat(:,:,(t-tm):(t+tc-1));
        [XhMPC(:,:,t+tc),kpi] = prepMPC(t,Xh(:,t-tm),Ufwd,Dfwd,Pi_c,...
            A,Bu,E,M,Nx,Nv,Nu,Np,tm,tc);
        
        solveStatus = 0;
        counter = 1;
        TMPC = Np;
        while(solveStatus==0)
            
            % compute U_{t+tc}^i, forall i s.t. {Pi_c(i,t) = 1}
            [Umpc,Jcomp(t),status,XP,~] = schedMPC(XhMPC(1:Nx,end,t+tc),...
                XhMPC((Nx+1):end,end,t+tc),kpi,TMPC,A,Bu,M,E,Q,Qf,R,...
                umax,umin,xmin,xmax,[]);
            
            if(strfind(status,'Solved'))
                solveStatus=1;
                fprintf('\nt=%d, MPC: %s\n',t,status)
                
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
    % reshape into MJLS form (for step t), uses D_t(pi(t-tc),alpha_c(t-tc))
    if(t<=tc)
        D_c(:,:,t) = makeD_c(zeros(1,Nv),zeros(1,Nv),Nu,Np);
    else
        D_c(:,:,t) = makeD_c(Pi_c(:,t-tc),alpha_c(:,t-tc),Nu,Np);
    end
    I = eye(size(D_c(:,:,1)));
    AA = [A,Bu*E*M*(I-D_c(:,:,t));zeros(Nv*Np*Nu,Nx),M*(I-D_c(:,:,t))];
    BU = [Bu*E*D_c(:,:,t);D_c(:,:,t)];
    
    % propagate system x_{t+1} = f(D_t,x_t,U_t,w_t)
    X(:,t+1) = AA*X(:,t) + BU*U(:,t) + BW*w(:,t);
    
    % actual applied u_t for saving
    u(:,t) = E*X(Nx+1:end,t+1);
    
    % measurement z_{t+1}:
    y(:,t+1) = C*X(1:Nx,t+1) + v(:,t+1);
    
    looptime(t) = toc(looptic);
    
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


% add to results struct for saving:

results.X = X;
results.u = u;
results.U = U;

results.uNoLoss = uNoLoss;
results.bNoLoss = bNoLoss;
results.tNoACK = tNoACK;
results.tNoACKSave = tNoACKSave;

results.Xh = Xh;
results.P = P;

results.XhMPC = XhMPC;
results.Jcomp = Jcomp;
results.XPlan = XPlan;
results.MPCtime = MPCtime;
results.MPCFail = MPCFail;

results.looptime = looptime;
results.Jsim = Jsim/Ns;
results.rmsEstError = nanrms(X(1:Nv,:) - Xh(1:Nv,:),2);
results.rmsPosError = nanrms(X(1,:),2);

% run-specific parameters
results.P1 = P1;
results.v = v;
results.w = w;
results.xIC = xIC;
results.xHat1 = xHat1;
results.alpha_c = alpha_c;
results.alpha_m = alpha_m;
results.alpha_a = alpha_a;
results.Pi_c = Pi_c;
results.Pi_m = Pi_m;
results.Pi_a = Pi_a;

% system
sys.A = A;
sys.Bu = Bu;
sys.C = C;
sys.umax = umax;
sys.umin = umin;
sys.Q = Q;
sys.Qf = Qf;
sys.R = R;
sys.Np = Np;
sys.W = W;
sys.V = V;
sys.alpha_cBar = alpha_cBar;
sys.ta = ta;
sys.tc = tc;
sys.tm = tm;
sys.nACKHistory = nACKHistory;

results.sys = sys;

end





