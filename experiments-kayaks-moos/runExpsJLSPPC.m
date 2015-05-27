% Run JLS-PPC pursuit experiments
% Standalone script, adapted from simJLSPPC

% 3 vehicles, fake field measurements
% perturbs generated before, and loaded here

% MOOS integration through Josh's timing code
% (preload origin and direction of C2 lines)

% Grab measurements (vehicle nav through acomms)
% (make field meas from nav)
% run estimator, controller
% Post quantized velocity command trajectories

% BR, 9/15/2014

clear variables global
clear mexmoos
%close all
clc
format compact

printouts = 1;  % some debugging/monitoring to terminal
ifDiary = 1;    % save a diary of terminal
ifSave = 1;     % auto save workspace at end of run
    
%-----MEXMOOS-----%
serverhost = '192.168.1.100';
%serverhost = 'localhost';
mexmoos('init','SERVERHOST',serverhost,'SERVERPORT','9000',...
    'MOOSNAME','mexmoos');
pause(1);
mexmoos('REGISTER','PURSUIT_VEHICLE_REPORT',0);

mexmoos('REGISTER','PURSUIT_X_NOSTROMO',0)
mexmoos('REGISTER','PURSUIT_Y_NOSTROMO',0)
mexmoos('REGISTER','PURSUIT_X_SILVANA',0)
mexmoos('REGISTER','PURSUIT_Y_SILVANA',0)
mexmoos('REGISTER','PURSUIT_X_KESTREL',0)
mexmoos('REGISTER','PURSUIT_Y_KESTREL',0)

mexmoos('REGISTER','PURSUIT_COMMAND_RECEIVED_NOSTROMO',0)
mexmoos('REGISTER','PURSUIT_COMMAND_RECEIVED_SILVANA',0)
mexmoos('REGISTER','PURSUIT_COMMAND_RECEIVED_KESTREL',0)

% OP system model and specific perturbations
% LOAD: 
%   Ap             6x6        
%   ApInt          6x6      
%   ApLoners      12x12           
%   Cp             3x6     
%   CpInt          3x6            
%   CpLoners       3x12             
%   Ns             1x1           
%   Nv             1x1              
%   Wn             6x1                 
%   Wp             6x6              
%   WpInt          6x6              
%   WpLoners      12x12              
%   Z              6x1                
%   gb             2x3x220               
%   kiCross        1x1               
%   pb             2x3x220         
%   wp             1x1             
%   wpt            6x220           
%   xp             6x220       

cd('/home/josh/Dropbox (MIT)/JLS MPC code/MOOS exps')

%load('expRefs_test3slow_Ns220_2014-Sep-19-14-39'); % slow perturbs

%load('expRefs_test1_Ns220_kc405e-4,wp35e-2_2014-Sep-23-11-51')
load('expRefs_test1_parallel_Ns220_kc405e-4,wp35e-2_2014-Sep-23-11-51')

figure
stairs((Cp*xp)')
xlabel('time step')
drawnow

%% SET UP SYSTEM PARAMS

%{
3 vehicle MX;  108 = 3 codebook entries
--> use multiples of this
108,216,324, etc.


MX 324 codebook = linspace(-8.4,8.4,9);
IL 324 codebook = linspace(-8.4,8.4,17);

%}

%sched = 'MXpiggyback_PS108';
%sched = 'MXpiggyback_PS324';
sched = 'ILpiggyback_PS108';
%sched = 'ILpiggyback_PS324';

addedLossProb = .2;

% currently always uses alphaBar state prior adjustment
covPriorAdj = 0;
if(covPriorAdj)
    disp('COV PRIOR ADJUST ON')
end

% vehicle disturbance variance
wq = .1*eye(Nv);

% sensor noise
% (NOTE: this sys uses gradients of 1)
Vphi = .1*eye(Nv);  % artificially added, AND used in KF
Vq = .1*eye(Nv);    % Vq for KF (vehicle nav error)
VqAdd = 0*eye(Nv);  % Vq to add onto qt projection
% Vstr = 'Vphi1_q25e-2';

% MPC params
Qr = blkdiag(10*eye(Nv),zeros(Nv));
Rr = 1*eye(Nv);
Qf_factor = 10;
NpMult = 3;

% constraints (speed in m/s * slot length in s)
umaxq = 1.2*7;
uminq = -1.2*7;

% packet success ESTIMATE (for jump estimator/KF)
alphaBar = [0.95,0.95,0.95]';

% delays
tm = 1;tc = 1;ta = 1;

ds = dateString('DHMS');
if(ifDiary)
    dfname = sprintf('JLSExpLog_%s_%s.txt',sched,ds);
    diary(dfname);
end

% initial states for estimator
P0s = 10;
xHat0p = 0*ones(Nv,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vehicle system model
Aq = eye(Nv);
Nxq = size(Aq,1);

% make block system
A = [Ap,zeros(Nv*2,Nv);zeros(Nv,Nv*2),Aq];
Nx = size(A,1);
Bu = [zeros(2*Nv,Nv);eye(Nv)];NU = size(Bu,2);
C = [eye(Nv),zeros(Nv),-eye(Nv);zeros(Nv,2*Nv),eye(Nv)];NZ = size(C,1);
V = blkdiag(Vphi,Vq);
W = blkdiag(Wp,wq);
NU = size(Bu,2);    % all controls
Nu = NU/Nv;         % per-vehicle
NZ = size(C,1);
Nz = NZ/Nv;

% set up for LQRY (output weighting)
tmp = [C', zeros(Nx,Nv);zeros(Nv,2*Nv) eye(Nv)]...
    *[Qr,zeros(2*Nv,Nv);zeros(Nv,2*Nv),Rr]*...
    [C,zeros(2*Nv,Nv);zeros(Nv,Nx),eye(Nv)];
Q = tmp(1:Nx,1:Nx);
R = tmp(Nx+1:end,Nx+1:end);
Qf = Q*Qf_factor;

umax = umaxq*ones(Nv,1);
umin = uminq*ones(Nv,1);
Xmax = [];Xmin = [];

P0 = P0s*eye(size(C,2));
xHat0 = [Cp'*xHat0p;zeros(Nv,1)]; % hat uses CpLoners

% read in comms params from sched
tmp = textscan(sched,'%*s PS%s','Delimiter','_');
if(strfind(cell2mat(tmp{1}),'Inf'))
    packetSize = Inf;
else
    packetSize = str2num(cell2mat(tmp{1}));
end

% RIGHT NOW ONLY WORKS WITH PIGGYBACK SCHEDS
if(strfind(sched,'IL'))
    
    [Pi,Xi,Lambda,tap,TsIL] = createSchedule(sched,Nv,Ns+2*Nv,tc);
    Np = NpMult*TsIL;  % MPC horizon
    
    if(isinf(packetSize))
        codebook = 'none';
    else
        nCommandsIL = Np*(NU/Nv);
        nLevels = packetSize/nCommandsIL;   
        % make sure nLevels is odd (so zero in codebook)
        if(mod(nLevels,2)==0);nLevels = nLevels -1;end
        fprintf('\n %s nLevels = %d \n',sched,nLevels)
        codebook = linspace(uminq,umaxq,nLevels);
    end
    
elseif(strfind(sched,'MX'))
    
    [Pi,Xi,Lambda,tap,TsMX] = createSchedule(sched,Nv,Ns+2*Nv,tc);
    Np = NpMult*TsMX;  % MPC horizon
    if(isinf(packetSize))
        codebook = 'none';
    else
        nCommandsMX = Np*NU;
        nLevels = packetSize/nCommandsMX;   
        % make sure nLevels is odd (so zero in codebook)
        if(mod(nLevels,2)==0);nLevels = nLevels -1;end
        fprintf('\n %s nLevels = %d \n',sched,nLevels)
        codebook = linspace(uminq,umaxq,nLevels);
    end
    
end

Umax = repmat(umax,1,Ns+Np);
Umin = repmat(umin,1,Ns+Np);

if(~ischar(codebook))
    % Nearest-Neighbor quantizer: make partition
    partition = ( codebook(2:end)+codebook(1:(end-1)) )/2;
end

M = makeM(Nu,Np,Nv);
e1t = [1,zeros(1,Np-1)];
E1 = kron(eye(Nv),kron(e1t,eye(Nu)));

% preallocate
vehRepSave = cell(1,Ns);
XhMPC = NaN*zeros(Nv*Np*Nu+Nx,tm+tc+1,Ns);	% includes xHatMPC and bHatMPC
Jcomp = NaN*zeros(1,Ns);
XPlan = NaN*zeros(Nx,Np,Ns);
MPCtime = NaN*zeros(Ns,1);
looptime = zeros(Ns,1);
MPCFail = zeros(Ns,1);
trueNavSave = zeros(2,3,Ns);
cmdAckSave = zeros(3,Ns);

% packet loss sequences (will be filled in real-time)
beta = zeros(Nv,Ns);alpha = zeros(Nv,Ns);gamma = zeros(Nv,Ns);
addedMeasLossesSave = zeros(Nv,Ns);
addedControlLossesSave = zeros(Nv,Ns);

% vehicle positions (as known to estimator)
qxy = zeros(2,Nv,Ns);
% vehicle perturbs
qt = zeros(Nv,Ns);
% true field meas (no noise added)
zPhi = zeros(Nv,Ns);
% true full meas
y = zeros(NZ,Ns);           
% y into estimator 
yh = zeros(NZ,Ns);          
% estimate
Xh = NaN*zeros(Nx+Nu*Np*Nv,Ns); % includes xHat and bHat
P = zeros(Nx,Nx,Ns);
% control action, plan, and encoded plan
u = zeros(NU,Ns);
U = zeros(Nu*Np*Nv,Ns);
UInd = zeros(Nu*Np*Nv,Ns);

% initialize Dhat and Dyes
Dh = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns);   % D hat for estimator
Ds = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns);   % used to track utilde for cov. prior
alphat = repmat(alphaBar,[1 Ns]);   % estimate of alpha, used in Dh
for t = (tc+1):Ns
    Dh(:,:,t) = makeD(Pi(:,t-tc),alphat(:,t-tc),Nu,Np);
    Ds(:,:,t) = makeD(Pi(:,t-tc),Pi(:,t-tc),Nu,Np);
end
D = zeros(Nu*Np*Nv,Nu*Np*Nv,Ns);    % true D
S = zeros(Nv*Nz,Nv*Nz,Ns);          % true S (which is known to est.)
% ack indicator a
a = zeros(Nv,Nv,Ns);
% utilde and bYes for updating copy of planned buffers
utilde = zeros(NU,Ns);
bYes = zeros(Nu*Np*Nv,1);

% white sensor noise
Vphys = blkdiag(Vphi,VqAdd);
v = real(sqrtm(Vphys))*randn(NZ,Ns);

% vehicle disturbances -- NOT SIMULATED
% wqt = real(sqrtm(wq))*randn(Nv,Ns);
% w = [wpt;wqt];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

simstart = tic;

P(:,:,1) = P0;
Xh(:,1) = [xHat0;zeros(Nu*Np*Nv,1)];

% post reference points to pMarineViewer
% (change if reference is moving)
pbLabel = {'pb1','pb2','pb3'};
for i = 1:Nv
    postStrRef = sprintf('type=circle,color=blue,x=%f,y=%f,label=%s',...
        pb(1,i,1),pb(2,i,1),pbLabel{i});
    mexmoos('NOTIFY','VIEW_MARKER',postStrRef)
end

% loop starts at "estimator" at time t=2 when first meas (t=1) is RX'd
for t=2:(Ns-1)
    looptic = tic;
    
    % read in measurements from vehicles
        
        % sync start of mission -- wait until vID=1
        % t=2 --> t-tm=1, and Xi(1,t-tm)==1
        if(t==2)
            gotStart = 0;
            while(~gotStart)
                [vehicleReport,vehRepOut,trueNav,cmdAck] = grabVehicleReport;
                vID = vehRepOut{1};
                fprintf('\n Waiting to start: vID=%d\n',vID)
                if(vID==1); gotStart=1; simstartTime = datestr(now); end
            end
        else
            % always syncs on vehicle report
            [vehicleReport,vehRepOut,trueNav,cmdAck] = grabVehicleReport;
            vID = vehRepOut{1};
            if(vID==-1)
                vID = [];   % command slot vID -- set to blank
                disp('veh report: command slot')
            end
            vehRepSave{t} = vehicleReport;
            trueNavSave(:,:,t) = trueNav;
            cmdAckSave(:,t) = cmdAck';
            fprintf('\ncmdAck:')
            disp(cmdAck)
        end
        
    
    %%% GRAB beta(:,t-tm) and vehicle nav from MOOS
    % construct y(:,t-tm) from vehicle nav
    % set gamma(:,t-ta) = beta(:,t-tm)
    
    if(vID)
        beta(vID,t-tm) = vehRepOut{2};
        
        % ADD SOME RANDOM PACKET LOSS
        addedLosses = rand(3,1)<addedLossProb;
        disp('meas added losses')
        disp(addedLosses)
        addedMeasLossesSave(:,t-tm) = addedLosses;
        
        beta(:,t-tm) = beta(:,t-tm).*(1-addedLosses);
        
        gamma(vID,t-tm) = beta(vID,t-tm);
        ACK = vehRepOut{3};
        
        if( (t-ta-tap(vID)-tc) > 0 )
            alpha(vID,t-ta-tap(vID)-tc) = ACK;
        end
        
        cycleCount = vehRepOut{4};
        
        % cycle should equal t (cycle @ estimator) - 1
        if(printouts)
            fprintf('\nCYCLE COUNT NOW (AT ESTIMATOR): %d \n',t)
            fprintf('PURSUIT_VEHICLE_REPORT:\n%s\n',vehicleReport)
            fprintf('VEHICLE %d BETA(t-tm): %d\n',vID,beta(vID,t-tm))
            fprintf('VEHICLE %d ACK(t-tm): %d\n',vID,ACK)
        end
        
        if(beta(vID,t-tm))
            tmp = textscan(char(vehRepOut{5}),'%d,%f,%f');
            qxy(1,vID,t-tm) = tmp{2};
            qxy(2,vID,t-tm) = tmp{3};
        end
        
    end
    
    % convert to scalar qtilde
    % use pbar, gbar, ptilde, and qtilde to make (fake) field measurement
    for i = 1:Nv
        [qt(i,t-tm),zPhi(i,t-tm),errq] = pursuitExpMeas(...
            qxy(1,i,t-tm),qxy(2,i,t-tm),pb(:,i,t-tm),gb(:,i,t-tm),...
            Cp(i,:)*xp(:,t-tm));
    end
    y(:,t-tm) = [zPhi(:,t-tm);qt(:,t-tm)]+v(:,t-tm);

    % at step t, meas. are sent at t-tm
    S(:,:,t-tm) = makeS(Xi(:,t-tm),beta(:,t-tm),Nz);
    yh(:,t-tm) = S(:,:,t-tm)*y(:,t-tm);     % available tm steps after sent
    
    % determine ACKs available at this step
    % update Dh (and alphat), and KFstart
    [Dh,alphat,a,KFstart] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
        Lambda,gamma,t,tm,tc,ta,tap,Nu,Np);    
    
    %%%%%%%
    % run estimator
    %%%%%%%
    
    % first compute utilde (control if packets all successful)
    % used if no ACK available
    if( (t-tm-1) >= 1)
        bYes = M*(eye(Np*Nu*Nv)-Ds(:,:,t-tm-1))*bYes + ...
            Ds(:,:,t-tm-1)*U(:,t-tm-1);
        utilde(:,t-tm-1) = E1*bYes;
    end
    
    % run KF from KFstart up until time of recent measurement
    for td = KFstart:(t-tm)
        if(td<=1)
            % first step - no dynamics
            AKF = eye(size(A));
            DKFh = makeD(zeros(Nv,1),zeros(Nv,1),Nu,Np);
            XhIn = [xHat0;zeros(Nu*Np*Nv,1)];
            Pin = P0;
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
    
    if(printouts)
        
        disp('True  qtilde for t-tm')
        disp(beta(:,t-tm).*qt(:,t-tm))
        disp('Est. qtildeHat for t-tm')
        disp(Xh((2*Nv+1):(3*Nv),t-tm))
        
        disp('True ptilde for t-tm')
        disp(Cp*xp(:,t-tm))
        disp('Est. ptildeHat for t-tm')
        disp(Cp*Xh(1:(2*Nv),t-tm))
        
    end

    % POST TO PMARINEVIEWER
    ptLabel = {'pt1','pt2','pt3'};
    qtLabel = {'qt1','qt2','qt3'};
    ptHatLabel = {'pHat1','pHat2','pHat3'};
    for i = 1:Nv
        pt = Cp(i,:)*xp(:,t-tm);
        ptHat = Cp(i,:)*Xh(1:(Nx-Nxq),t-tm);
        [xi,xhi,qi] = pursuitExpPost(pb(:,i,t-tm),gb(:,i,t-tm),...
            pt,ptHat,qt(i,t-tm));
        postStrQ = sprintf('type=circle,color=magenta,x=%f,y=%f,label=%s',...
            qi(1),qi(2),qtLabel{i});
        postStrX = sprintf('type=circle,color=red,x=%f,y=%f,label=%s',...
            xi(1),xi(2),ptLabel{i});
        postStrXh = sprintf('type=diamond,color=silver,x=%f,y=%f,label=%s',...
            xhi(1),xhi(2),ptHatLabel{i});

        % only post the actually updated qtilde
        if(beta(i,t-tm))
            mexmoos('NOTIFY','VIEW_MARKER',postStrQ)
        end
        mexmoos('NOTIFY','VIEW_MARKER',postStrX)
        mexmoos('NOTIFY','VIEW_MARKER',postStrXh)

    end
       
    %%%%%%%
    % if control is to be computed/sent this step
    if(max(Pi(:,t)))
        MPCtic = tic;
        
        % grab constraints (to accommodate time-varying)
        [umax,umin,xmax,xmin] = paramsNow(Umax,Umin,Xmax,Xmin,t+1,Np);
        
        % Forward propagation: XhMPC and k_p^i's
        % starts with Xh: xHat_{t-tm|t-tm}, bHat_{t-tm-1}
        % compute xHatMPC_{t-tm+1:t+tc|t-tm}, bHatMPC_{t-tm:t+tc-1}
        % first step: uses uHat_{t-tm} <-- Dh(:,:,t-tm)
        
        Ufwd = U(:,(t-tm):(t+tc-1));
        Dfwd = Dh(:,:,(t-tm):(t+tc-1));
        [XhMPC(:,:,t+tc),kpi] = prepMPC(t,Xh(:,t-tm),Ufwd,Dfwd,Pi,...
            A,Bu,E1,M,Nx,Nv,Nu,Np,tm,tc);
        
        if(printouts)
            
            disp('Est buffers b_(t-tm-1)')
            buffDisp = Xh( (3*Nv+1):end,t-tm);
            buffDisp = reshape(buffDisp,[Nu*Np,3]);
            disp(buffDisp)
            
            %disp('Fwd propagation XhMPC for t+tc')
            %disp(XhMPC(:,:,t+tc))
            
        end
        
        solveStatus = 0;
        counter = 1;
        TMPC = Np;
        while(solveStatus==0)
            
            % compute U_{t+tc}^i, forall i s.t. {pi(i,t) = 1}
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
        
        % translate (NU x Np Umpc) into buffer shape Nu x Np x Nv
        % quantize UMat --> UIndMat
        Uvec = reshape(Umpc',[Nu*Np*Nv,1]);
        if(isempty(strfind(codebook,'none')))
            [UInd(:,t+tc),U(:,t+tc)] = quantiz(Uvec,partition,codebook);
            UIndMat = reshape(UInd(:,t+tc),[Nu*Np,Nv])-floor(length(codebook)/2);
            UMat = reshape(U(:,t+tc),[Nu*Np,Nv]);
            if(printouts)
                disp('Control plans')
                disp(UMat)
            end
        elseif(strfind(codebook,'none'))
            U(:,t+tc) = Uvec;
        end
        
        % POST TO MOOS FOR SENDING OVER ACOMMS
        % your command:id,comma delimited integers
        
        % find which vehicles to send commands to this step
        sendInds = find(Pi(:,t));
        
        addedLosses = rand(length(sendInds),1)<addedLossProb;
        disp('control added losses')
        disp(addedLosses)
        addedControlLossesSave(:,t) = addedLosses;
        
        sendInds = sendInds(logical(1-addedLosses))
        vID
        
        % send commands
        for i = 1:length(sendInds)
            % post command to MOOS (grabbed by timing script)
            cmdStr = sprintf(',%d',UIndMat(:,sendInds(i)));
            postStr = sprintf('%d%s',sendInds(i),cmdStr);
            if(printouts)
                disp('POSTING')
                disp(postStr)
            end
            mexmoos('NOTIFY','PURSUIT_VEHICLE_COMMAND',postStr)
        end

        % (post to MOOS log)
        % Xh(:,t-tm), XhMPC(1:Nx,end,t+tc)
        XhStr = sprintf(',%0.3f',Xh(:,t-tm));
        XhMPCStr = sprintf(',%0.3f',XhMPC(:,end,t+tc));
        estSaveString = sprintf('t,%d,Xh%s,XhMPC%s',t,XhStr,XhMPCStr);
        mexmoos('NOTIFY','PURSUIT_ESTIMATES',estSaveString);

    else % (for saving - make clear not set)
        Umpc = zeros(NU,Np);
    end
    
    looptime(t) = toc(looptic);
    % no timing here... synced to MOOS
    
end


%% cleanup

% save separate diary for each trial
if(ifDiary)
    diary off
end

if(ifSave)
    resfname = sprintf('JLSpursuitResults_%s_%s.mat',sched,ds);
    save(resfname)
end

