function [results] = simJLSPPC(Ns,Np,Asys,Busys,Bwsys,Csys,...
    QMPC,QfMPC,RMPC,WKF,VKF,tm,tc,ta,tac,acBar,Pi_c,Pi_m,Pi_a,Ts,...
    umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,w,v,a_c,a_m,a_a,...
    covPriorAdj,nACKHistory,printDebug)
% runs simulation of MJLS/scheduled PPC



% BR, 4/23/2014
% modifying for delayed ACKs, 6/13/2014
% v1.0 6/13/2015
% v1.1 6/16/2015

% TO DO:
% Clean up printouts?
% Make sure prep for KF is ready for MIMO and multiple P*s
% Add logging of P*, P**, etc.?
% Add more automated tests/checks?


% initialization stuff
Nsysx = size(Asys,1);
Nv = length(acBar);
NU = size(Busys,2);    % all controls
Nu = NU/Nv;         % per-vehicle
Nw = size(Bwsys,2);
NY = size(Csys,1);
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

if(covPriorAdj)
    disp('COV PRIOR ADJUST ON')
    % compute Pstar terms for specific alpha_cBar
    if(length(acBar)==1)
        
        % load in the symbolic lookup table - faster
        [PstarCoefficients,PstarFinalCoefficients] = ...
            evaluatePstars(acBar);
        cv.PstarCoefficients = PstarCoefficients;
        cv.PstarFinalCoefficients = PstarFinalCoefficients;
        
        % optionally, if nStars>11 desired, can compute numerically for a
        % specific alpha_cBar (may take a while)
        %{
        nStars = 12;
        PstarCoefficients = computePstars(nStars,alpha_cBar,'NUMERIC');
        %}
        
    else
        %
        disp('USING MULTIVEHICLE ONE-STEP P* APPROX')
    end
    
end


Umax = repmat(umax,1,Ns+Np);
Umin = repmat(umin,1,Ns+Np);

if(~ischar(codebook))
    % Nearest-Neighbor quantizer: make partition
    partition = ( codebook(2:end)+codebook(1:(end-1)) )/2;
end

M = makeM(Nu,Np,Nv);
etmp = [1,zeros(1,Np-1)];
E = kron(eye(Nv),kron(etmp,eye(Nu)));
BW = [Bwsys;zeros(Nu*Np*Nv,Nw)];

X = zeros(Nsysx+Nu*Np*Nv,Ns);  % includes x and b
y = zeros(NY,Ns);           % true y (could be reconstructed Cx+v)
yh = zeros(NY,Ns);          % y into estimator (could be reconstructed)
Xh = NaN*zeros(Nsysx+Nu*Np*Nv,Ns); % includes xHat and bHat
P = zeros(Nsysx,Nsysx,Ns);
u = zeros(NU,Ns);
U = zeros(Nu*Np*Nv,Ns);

% initialize Dhat and Dyes
D_cHat = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns+tc);   % D hat for estimator
D_cNoLoss = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns+tc);   % used to track utilde for cov. prior
alphaHat = repmat(acBar,[1 Ns]);
for t = (tc+1):Ns
    D_cHat(:,:,t) = makeD_c(Pi_c(:,t-tc),alphaHat(:,t-tc),Nu,Np);
    D_cNoLoss(:,:,t) = makeD_c(Pi_c(:,t-tc),Pi_c(:,t-tc),Nu,Np);
end

D_c = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns);    % true D
D_m = zeros(Nv*Ny,Nv*Ny,Ns);
D_a = zeros(Nv,Nv,Ns);

uNoLoss = zeros(NU,Ns);
bNoLoss = zeros(Nu*Np*Nv,Ns);

XhMPC = NaN*zeros(Nv*Np*Nu+Nsysx,tm+tc+1,Ns);	% includes xHatMPC and bHatMPC
Jcomp = NaN*zeros(1,Ns);
XPlan = NaN*zeros(Nsysx,Np,Ns);
MPCtime = NaN*zeros(Ns,1);
looptime = zeros(Ns,1);
MPCFail = zeros(Ns,1);

% initial state x_1 and z_1
X(1:Nsysx,1) = xIC;
y(:,1) = Csys*xIC+v(:,1);

% initial buffer - all zeros
X(Nsysx+1:end,1) = zeros(Nu*Np*Nv,1);

% first step propagation - gives x_2
X(1:Nsysx,2) =  Asys*xIC + w(:,1);
y(:,2) = Csys*X(1:Nsysx,1) + v(:,2);
u(:,1) = E*X(Nsysx+1:end,1);

% Loop starts at "estimator" at time t=tm+1, when first meas (t=1) is RX'd
for t = (tm+1):(Ns-1)
    looptic = tic;
    
    % determine which measurements are available at estimator at this step
    % at step t, meas. are sent at t-tm
    D_m(:,:,t-tm) = makeD_m(Pi_m(:,t-tm),a_m(:,t-tm),Ny);
    yh(:,t-tm) = D_m(:,:,t-tm)*y(:,t-tm);     % available tm steps after sent
    
    if(printDebug)
        fprintf('\n~~~STEP t=%d AT ESTIMATOR~~~\n',t)
        for i = 1:Nv
            if(Pi_m(i,t-tm)*a_m(i,t-tm)==1)
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
        Pi_c,D_a,a_c,alphaHat,Pi_a,Ts,a_a,t,tm,tc,ta,tac,...
        Nu,Np,tNoACK,nACKHistory,printDebug);
    tNoACKSave{t} = tNoACK;
    
    if(printDebug)
        dispStart = t-8;
        dispEnd = t+Ts+2;
        if(dispStart<1);dispStart=1;end
        if(dispEnd>Ns);dispEnd=Ns;end
        fprintf('\n tNoACK(%d:%d):',dispStart,dispEnd)
        disp(tNoACK(:,dispStart:dispEnd))
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
                    tNoACK_KF = tNoACK(1,tKF+tac-1);
                else
                    tNoACK_KF = 1;
                end
            else
                
                % KEEP THIS?  OR MODIFY FOR JUST >0 OR NOT
                %{
                tNoACK_KF = zeros(1,Nv);
                for i = 1:Nv
                    if(td-1-tac(i)>0)
                        tNoACK_KF(i) = tNoACK(i,td+tac(i)-1);
                    end
                end
                %} 

                % (artificially constrain dU to zero for ACK'd channels?)
                % (DO THIS HERE? OR INSIDE KF?)
                
                tNoACK_KF = 1;
                
            end
            
            % prepare control options for use in cov. prior adj.
            % Uoptions is indexed backwards
            % Uoptions(:,1) is most recent, (:,tNoACK_KF) is furthest back
            uOptions = zeros(NU,tNoACK_KF);
            for k = 1:tNoACK_KF
                if(tKF-k<1)
                    bTMP = zeros(Nu*Np*Nv,1);
                else
                    bTMP = bNoLoss(:,tKF-k);
                end
                uOptions(:,k) = E*M^k*bTMP;
            end
            
        else
            
            uOptions = [];
            tNoACK_KF = [];
            
        end
        
        if(tKF<=1)
            AKF = eye(size(Asys));
            DKFh = makeD_c(zeros(Nv,1),zeros(Nv,1),Nu,Np);
            XhIn = [xHat1;zeros(Nu*Np*Nv,1)];
            Pin = P1;
            Uin = zeros(Nu*Np*Nv,1);
            yIn = yh(:,1);
            SIn = D_m(:,:,1);
        else
            AKF = Asys;
            DKFh = D_cHat(:,:,tKF-1);
            XhIn = Xh(:,tKF-1);
            Pin = P(:,:,tKF-1);
            Uin = U(:,tKF-1);
            yIn = yh(:,tKF);
            SIn = D_m(:,:,tKF);
        end
        
        if(covPriorAdj)
            cv.uOptions = uOptions;
            cv.tNoACK = tNoACK_KF;
        else
            cv=0;
        end
        
        printDebugKF=0;
        if(printDebugKF~=0)
            printDebugKF.t = t;printDebugKF.tKF = tKF;
        end
            
        % Xh(:,t-tm): xHat_{t-tm|t-tm},bHat_{t-tm-1}
        [Xh(:,tKF),P(:,:,tKF)] = JLSKF(XhIn,Pin,yIn,Uin,DKFh,...
            Nsysx,Nv,Nu,Np,SIn,AKF,Busys,E,M,Csys,WKF,VKF,acBar,...
            cv,printDebugKF);
        
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
            if(size(Asys,1)==1)
                % only print estimate for scalar sys
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
        [XhMPC(:,:,t+tc),p_i] = prepMPC(t,Xh(:,t-tm),Ufwd,Dfwd,Pi_c,...
            Asys,Busys,E,M,Nsysx,Nv,Nu,Np,tm,tc);
              
        solveStatus = 0;
        counter = 1;
        TMPC = Np;
        while(solveStatus==0)
            
            % compute U_{t+tc}^i, forall i s.t. {Pi_c(i,t) = 1}
            [Umpc,Jcomp(t),status,XP,~] = schedMPC(XhMPC(1:Nsysx,end,t+tc),...
                XhMPC((Nsysx+1):end,end,t+tc),p_i,TMPC,Asys,Busys,M,E,QMPC,QfMPC,RMPC,...
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
        D_c(:,:,t) = makeD_c(Pi_c(:,t-tc),a_c(:,t-tc),Nu,Np);
    end
    I = eye(size(D_c(:,:,1)));
    AA = [Asys,Busys*E*M*(I-D_c(:,:,t));zeros(Nv*Np*Nu,Nsysx),M*(I-D_c(:,:,t))];
    BU = [Busys*E*D_c(:,:,t);D_c(:,:,t)];
    
    % propagate system x_{t+1} = f(D_t,x_t,U_t,w_t)
    X(:,t+1) = AA*X(:,t) + BU*U(:,t) + BW*w(:,t);
    
    % actual applied u_t for saving
    u(:,t) = E*X(Nsysx+1:end,t+1);
    
    % measurement z_{t+1}:
    y(:,t+1) = Csys*X(1:Nsysx,t+1) + v(:,t+1);
    
    looptime(t) = toc(looptic);
    
end

xF = X(1:Nsysx,end);
% compute "actual" cost
jj = u(:,1)'*RMPC*u(:,1);  % control at step 1
for t = 2:Ns
    % states at step 2 (affected by u(1)), through Ns
    jj = jj + X(1:Nsysx,t)'*QMPC*X(1:Nsysx,t) + u(:,t)'*RMPC*u(:,t);
end
% final state (affected by final u)
Jsim = jj + xF'*QMPC*xF;


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
results.rmsEstError = nanrms(X(1:Nsysx,:) - Xh(1:Nsysx,:),2);
results.rmsStateError = nanrms(X(1:Nsysx,:),2);

% run-specific parameters
results.P1 = P1;
results.v = v;
results.w = w;
results.xIC = xIC;
results.xHat1 = xHat1;
results.alpha_c = a_c;
results.alpha_m = a_m;
results.alpha_a = a_a;
results.Pi_c = Pi_c;
results.Pi_m = Pi_m;
results.Pi_a = Pi_a;

% system
sys.A = Asys;
sys.Bu = Busys;
sys.C = Csys;
sys.umax = umax;
sys.umin = umin;
sys.Q = QMPC;
sys.Qf = QfMPC;
sys.R = RMPC;
sys.Np = Np;
sys.W = WKF;
sys.V = VKF;
sys.alpha_cBar = acBar;
sys.ta = ta;
sys.tc = tc;
sys.tm = tm;
sys.nACKHistory = nACKHistory;

results.sys = sys;

end





