% Script to run LONERS/LOWERBOUND pursuit experiments
% 3 vehicles, fake field measurements

% BR, 9/12/2014

% perturbs generated before, and loaded here

% grab a basic vehicleReport
% run estimator and controller
% post commands in same way


clear variables global
clear mexmoos
%close all
clc
format compact

cd('/home/josh/Dropbox (MIT)/JLS MPC code/MOOS exps')

printouts = 1;
ifDiary = 1;
ifSave = 1;

%-----MEXMOOS-----%
serverhost = '192.168.1.100';
%serverhost = 'localhost';
mexmoos('init','SERVERHOST',serverhost,'SERVERPORT','9000',...
    'MOOSNAME','mexmoos');
pause(1);
mexmoos('REGISTER','PURSUIT_VEHICLE_REPORT',0);


slotLength = 7;

% OP system model and specific perturbations
% load Ap,Cp,Wp,ApLoners,CpLoners,WpLoners,Wn,Z,ApInt,CpInt,WpInt
% and xp, wpt
% set up origin (reference) and direction (unit vector) of C2 lines
%load('/home/josh/Dropbox (MIT)/JLS MPC code/MOOS exps/expRefs_test1_Ns220_2014-Sep-17-14.mat')
%load('expRefs_test1_Ns220_kc405e-4,wp35e-2_2014-Sep-23-11-51')
load('expRefs_test1_parallel_Ns220_kc405e-4,wp35e-2_2014-Sep-23-11-51')



figure
stairs((Cp*xp)')
xlabel('time step')
drawnow

%% SET UP SYSTEM PARAMS

%sched = 'Loners';
%sched = 'LonersInt';
%sched = 'LowerBound';

%sched = 'Loners_sched';
sched = 'LonersInt_sched';
%sched = 'LowerBound_sched';

vehicleNavSched = 1; % if vehicle nav follows meas schedule
fprintf('\n Vehicle nav sched: %d\n',vehicleNavSched)

% vehicle disturbance variance
wq = .1*eye(Nv);

% sensor noise
% (NOTE: this sys uses gradients of 1)
Vphi = .1*eye(Nv);
Vq = .1*eye(Nv);
% Vstr = 'Vphi1_q25e-2';

% MPC params
Qr = blkdiag(10*eye(Nv),zeros(Nv));
Rr = 1*eye(Nv);
Qf_factor = 10;
Np = 10; % make shorter eventually?

% constraints
%umax = 12.8*ones(Nv,1);
%umin = -12.8*ones(Nv,1);
umax = 1.2*7*ones(Nv,1);
umin = -1.2*7*ones(Nv,1);

Xmax = [];Xmin = [];

ds = dateString('DHMS');
if(ifDiary)
    dfname = sprintf('JLSExpLog_NoComms_%s_%s.txt',sched,ds);
    diary(dfname);
end

% initial states for estimator
P0s = 10;
xHat0p = 0*ones(Nv,1);

% vehicle system model
Aq = eye(Nv);
Nxq = size(Aq,1);

% make block system/LQRY etc
if(strfind(sched,'Loners'))
    
    % LONERS system
    CLoners = [CpLoners,-eye(Nv);zeros(Nv,size(CpLoners,2)),eye(Nv)];
    NxpL = size(CpLoners,2);
    QrL = blkdiag(10*eye(Nv),zeros(Nv));
    tmp = [CLoners', zeros(NxpL+Nv,Nv);zeros(Nv,2*Nv) eye(Nv)]...
        *[QrL,zeros(2*Nv,Nv);zeros(Nv,2*Nv),Rr]*...
        [CLoners,zeros(2*Nv,Nv);zeros(Nv,NxpL+Nv),eye(Nv)];
    QLoners = tmp(1:(NxpL+Nv),1:(NxpL+Nv));
    RLoners = tmp((NxpL+Nv+1):end,(NxpL+Nv+1):end);
    QfLoners = QLoners*Qf_factor;
    WLoners = blkdiag(WpLoners,wq);
    
    CInt = [CpInt,-eye(Nv);zeros(Nv,size(CpInt,2)),eye(Nv)];
    NxpI = size(CpInt,2);
    QrI = blkdiag(10*eye(Nv),zeros(Nv));
    tmp = [CInt', zeros(NxpI+Nv,Nv);zeros(Nv,2*Nv) eye(Nv)]...
        *[QrI,zeros(2*Nv,Nv);zeros(Nv,2*Nv),Rr]*...
        [CInt,zeros(2*Nv,Nv);zeros(Nv,NxpI+Nv),eye(Nv)];
    QInt = tmp(1:(NxpI+Nv),1:(NxpI+Nv));
    RInt = tmp((NxpI+Nv+1):end,(NxpI+Nv+1):end);
    QfInt = QInt*Qf_factor;
    WInt = blkdiag(WpInt,wq);
    
    if(strfind(sched,'sched'))
        [~,XiNC,~,~,~] = createSchedule('RRmeas',Nv,Ns+2*Nv,0);
    else
        XiNC = ones(Nv,Ns+2*Nv);
    end
    
    % LONERS SIM
    if(strfind(sched,'Int'))
        %ApL = ApInt;
        ApModel = ApInt;
        CModel = CInt;
        CpModel = CpInt;
        W = WInt;
        Q = QInt;
        Qf = QfInt;
        R = RInt;
    else
        %ApL = ApLoners;
        ApModel = ApLoners;
        CModel = CLoners;
        CpModel = CpLoners;
        W = WLoners;
        Q = QLoners;
        Qf = QfLoners;
        R = RLoners;
    end
    
    P0 = P0s*eye(size(CModel,2));  % est. uses CLoners size
    
elseif(strfind(sched,'LowerBound'))
    
    % make block system
    ApModel = Ap;
    Nx = size(ApModel,1)+Nv;
    %A = [Ap,zeros(Nv*2,Nv);zeros(Nv,Nv*2),Aq];
    %Bu = [zeros(2*Nv,Nv);eye(Nv)];
    %NU = size(Bu,2);
    CModel = [eye(Nv),zeros(Nv),-eye(Nv);zeros(Nv,2*Nv),eye(Nv)];
    %NZ = size(C,1);
    Bw = eye(Nx);nw = size(Bw,2);
    
    % set up for LQRY (output weighting)
    tmp = [CModel', zeros(Nx,Nv);zeros(Nv,2*Nv) eye(Nv)]...
        *[Qr,zeros(2*Nv,Nv);zeros(Nv,2*Nv),Rr]*...
        [CModel,zeros(2*Nv,Nv);zeros(Nv,Nx),eye(Nv)];
    Q = tmp(1:Nx,1:Nx);
    R = tmp(Nx+1:end,Nx+1:end);
    Qf = Q*Qf_factor;
    W = blkdiag(Wp,wq);
    
    % MPC horizon: multiple of schedule
    if(strfind(sched,'sched'))
        [~,XiNC,~,~,~] = createSchedule('RRmeas',Nv,Ns+2*Nv,0);
    else
        XiNC = ones(Nv,Ns+2*Nv);
    end
    
    % Sim coupled model with perfect comms
    P0 = P0s*eye(size(CModel,2));
    
end

V = blkdiag(Vphi,Vq);

Nxq = size(Aq,1);
Nxp = size(Ap,1);
Nxpm = size(ApModel,1);
Nv = size(Aq,1);

Acontrol = [ApModel,zeros(Nxpm,Nxq);zeros(Nxq,Nxpm),Aq];

Nxm = size(Acontrol,1);
NU = Nv;    % all controls
Nz = 2;
NZ = 2*Nv;

% model B
Bum = [zeros(Nxpm,NU);eye(NU)];

Umax = repmat(umax,1,Ns+Np);
Umin = repmat(umin,1,Ns+Np);


% vehicle true positions
qxy = zeros(2,Nv,Ns);
% vehicle perturbs
qt = zeros(Nv,Ns);
% true field meas (no noise added)
zPhi = zeros(Nv,Ns);
y = zeros(NZ,Ns);           % true (could be reconstructed Cx+v)

Xh = NaN*zeros(Nxm,Ns);      % includes xHat and bHat
P = zeros(Nxm,Nxm,Ns);
u = zeros(NU,Ns);
S = zeros(NZ,NZ,Ns);

Jcomp = NaN*zeros(1,Ns);
XPlan = NaN*zeros(Nxm,Np,Ns);
MPCtime = NaN*zeros(Ns,1);
looptime = zeros(Ns,1);
MPCFail = zeros(Ns,1);

% white sensor noise
v = real(sqrtm(V))*randn(NZ,Ns);

vehicleNames = {'NOSTROMO','SILVANA','KESTREL'};
simstart = tic;
simstartTime = datestr(now)

% Loop starts at "estimator" at time t=1, when first meas (t=1) is RX'd

% post reference points to pMarineViewer
% (change if reference is moving)
pbLabel = {'pb1','pb2','pb3'};
for i = 1:Nv
    postStrRef = sprintf('type=circle,color=blue,x=%f,y=%f,label=%s',...
        pb(1,i,1),pb(2,i,1),pbLabel{i});
    mexmoos('NOTIFY','VIEW_MARKER',postStrRef)
end

for t = 1:(Ns-1)
    fprintf('\n\nStep %d\n',t)
    looptic = tic;
    
    % READ IN VEHICLE REPORT - VEHICLE NAV
    % id,x,y:id,x,y:id,x,y

    % always waiting on report
    gotRep = 0;
    while(~gotRep)
        % grab mail from mexmoos
        % (most recent messages are lower indices in the buffer)
        msgs=mexmoos('FETCH');
        for k = 1:length(msgs)
            msgStr = msgs(k).KEY;
            if(strcmp(msgStr,'PURSUIT_VEHICLE_REPORT'))
                vehicleReport = msgs(k).STR;
                gotRep = 1;
                fprintf('\nGot Vehicle Report')
                vehRepOut = textscan(vehicleReport,...
                    '%s %s %s','Delimiter',':');
            end
        end
    end
    
    % pull out qxy and save
    for i = 1:length(vehRepOut)
        tmp = textscan(char(vehRepOut{i}),'%d %f %f','Delimiter',',');
        vID = tmp{1};
        qxy(1,vID,t) = tmp{2};
        qxy(2,vID,t) = tmp{3};   
    end
        
    % convert to scalar q
    % use G to make field measurement
    for i = 1:Nv
        [qt(i,t),zPhi(i,t),errq] = pursuitExpMeas(...
            qxy(1,i,t),qxy(2,i,t),pb(:,i,t),gb(:,i,t),Cp(i,:)*xp(:,t));
    end
    y(:,t) = [zPhi(:,t);qt(:,t)]+v(:,1);
    
    if(printouts)
        disp('Nav (x,y)')
        disp(qxy(:,:,t))
        disp('qt')
        disp(qt(:,t))
        disp('zPhi')
        disp(zPhi(:,t))
    end
    
    %%%%%%%
    % run estimator 
    %%%%%%%
    
    % meas schedule
    if(vehicleNavSched)
        S(:,:,t) = kron(eye(Nz),diag(XiNC(:,t)));
    else
        S(:,:,t) = blkdiag(diag(XiNC(:,t)),eye(Nv));
    end
    
    % covariance prior (standard)
    if(t==1)
        PP = P0;
        XXh = zeros(Nxm,1);
        uu = zeros(NU,1);
    else
        PP = P(:,:,t-1);
        XXh = Xh(:,t-1);
        uu = u(:,t-1);
    end
    Ppre = Acontrol*PP*Acontrol' + W ; %P_{t+1|t}
    
    % Kalman gain = fcn of Ppre
    L = ((Ppre*CModel')/(CModel*Ppre*CModel'+V))*S(:,:,t);   % L_{t}
    P(:,:,t) = Ppre - L*(CModel*Ppre);
    
    % propagate system
    I = eye(size(Acontrol));
    Xh(:,t) = (I-L*CModel)*(Acontrol*XXh + Bum*uu) + ...
        L*S(:,:,t)*y(:,t);
    
    %%%%%%%
    
    if(printouts)
        disp('True p for t')
        disp(Cp*xp(:,t))
        disp('Est. pHat for t')
        if(strfind(sched,'Int'))
            disp(CpInt*Xh(1:Nxpm,t))
        elseif(strfind(sched,'Loners'))
            disp(CpLoners*Xh(1:Nxpm,t))
        else
            disp(Cp*Xh(1:size(Ap,1),t))
        end
    end
    
    % POST TO PMARINEVIEWER
    %   VIEW_MARKER = "type=efield,x=100,y=20,scale=4.3,label=alpha,color=red,width=4.5"

    ptLabel = {'pt1','pt2','pt3'};
    qtLabel = {'qt1','qt2','qt3'};
    ptHatLabel = {'pHat1','pHat2','pHat3'};
    for i = 1:Nv
        pt = Cp(i,:)*xp(:,t);
        if(strfind(sched,'Int'))
            ptHat = CpInt*Xh(1:Nxpm,t);
        elseif(strfind(sched,'Loners'))
            ptHat = CpLoners*Xh(1:Nxpm,t);
        else
            ptHat = Cp*Xh(1:size(Ap,1),t);
        end
        [xi,xhi,qi] = pursuitExpPost(pb(:,i,t),gb(:,i,t),pt,...
            ptHat(i),qt(i,t));
        postStrQ = sprintf('type=circle,color=magenta,x=%f,y=%f,label=%s',...
            qi(1),qi(2),qtLabel{i});
        postStrX = sprintf('type=circle,color=red,x=%f,y=%f,label=%s',...
            xi(1),xi(2),ptLabel{i});
        postStrXh = sprintf('type=diamond,color=silver,x=%f,y=%f,label=%s',...
            xhi(1),xhi(2),ptHatLabel{i});

        % only post the actually updated qtilde
        mexmoos('NOTIFY','VIEW_MARKER',postStrQ)
        mexmoos('NOTIFY','VIEW_MARKER',postStrX)
        mexmoos('NOTIFY','VIEW_MARKER',postStrXh)


    end
    
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
    
%   For your command, just post as one string I guess:
%   id,speed:id,speed:id,speed
    cmdString = sprintf('1,%f:2,%f:3,%f',u(1,t),u(2,t),u(3,t));
    mexmoos('NOTIFY','PURSUIT_VEHICLE_COMMAND',cmdString);
    
    % log Xh to MOOS
    XhStr = sprintf(',%0.3f',Xh(:,t));
    estSaveString = sprintf('t,%d,Xh%s',t,XhStr);
    mexmoos('NOTIFY','PURSUIT_ESTIMATES',estSaveString);

    if(printouts)
        disp('COMMANDED q(x,y)')
        u(:,t)
    end
    
    % TIMING... WAIT ONE SLOT LENGTH
    looptime(t) = toc(looptic);
     
    % (Now, timing done by MOOS)
    %{
    waitTime = slotLength - looptime(t)
    while(toc(looptic)<slotLength)
        msgs=mexmoos('FETCH');
    end
    %}
end



%%

% save separate diary for each trial
if(ifDiary)
    diary off
end

if(ifSave)
    resfname = sprintf('pursuitResults_%s_%s.mat',sched,ds);
    save(resfname)
end


%{
    for i = 1:Nv
        
        % COORD XFORM... VELOCITY TO WAYPOINT USING C2 LINE AND DT
        %qcmd(:,i,t) = qxy(:,i,t) + u(i,t)*slotLength*gb(:,i,t);
        %qcmdStr = sprintf('points=%f,%f',qcmd(1,i,t),qcmd(2,i,t));
%         if(realtime)
%             mexmoos('NOTIFY',sprintf('PURSUIT_WAYPOINT_UPDATES_%s',vehicleNames{i}),qcmdStr)
%             mexmoos('NOTIFY',sprintf('PURSUIT_WAYPOINT_UPDATES_%s',vehicleNames{i}),sprintf('speed=%f',abs(u(i,t))))
%             mexmoos('NOTIFY',sprintf('PURSUIT_ACTION_%s',vehicleNames{i}),'WAYPOINT')
%         end
    end

%}