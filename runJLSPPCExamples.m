% Script to run a simple simulation of JLS-PPC
% Calls simJLSPPC (which calls functions in core)

% Basic setup:
% scalar or double integrator (SISO or MIMO) systems
% a few options for schedules
% handles varying lengths of ACK histories

% some simple plots for SISO systems at end

% BR, 5/27/2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all
clc

print_debug = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIM LENGTH
SIM_LENGTH = 25; % sim length

% MPC HORIZON:
N_HORIZON_MULT = 4; % the MPC horizon Np = Ts*NpMult

%%%%%%%%%%%
% SYSTEM (set up in setupSystemJPLPPC)

system = 'SISO_DOUBLE_INTEGRATOR';
%system = 'MIMO_DOUBLE_INTEGRATOR';
%system = 'SCALAR';

%%%%%%%%%%%
% SCHEDULE

% right now, "*_piggyback or *_noACK" are only options

if( ~isempty(strfind(system,'SISO')) || ~isempty(strfind(system,'SCALAR')))
    
    % 'SISOALL' - [1],[1] - std. discrete time
    % 'SISO2' - [1 0], [0 1] for pi, xi
    % 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi
    
    %sched = 'SISO4_piggyback';
    %sched = 'SISO4_noACK';
    
    %sched = 'SISO2_noACK';
    sched = 'SISO2_piggyback';
    
    %sched = 'SISOALL_piggyback';
    %sched = 'SISOALL_noACK';
    
elseif( ~isempty(strfind(system,'MIMO')))
    
    % 'MIMO' options: MX, IL
    
    sched = 'MX_piggyback';
    %sched = 'MX_noACK';
    
    %sched = 'IL_piggyback';
    %sched = 'IL_noACK';
    
end

% DELAYS
TAU_C = 1; % control delay
TAU_M = 1; % meas delay (Also ACK delay with piggyback unless overwritten)
TAU_A = 1; % ACK delay

% ACK SETTINGS
% Length of ACK history sent
% Should be 1 + a multiple of schedule length
% eg...  to ACK the newest control command, make N_ACKHISTORY = 1
%        to ACK newest and previous commands, make N_ACKHISTORY = 1 + T_S
%        to ACK the newest and 2 prev. cmds, make N_ACKHISTORY = 1 + 2*T_S
N_ACKHISTORY = 3;

% adjustment to covariance priors due to no ACKs/control losses:
cov_prior_adj = 1;

%%%%%%%%%%%%%
% PACKET SUCCESS PROBABILITIES

% Nv is number of control AND meas channels (asymmetric # not implemented)
if( ~isempty(strfind(system,'SISO')) || ~isempty(strfind(system,'SCALAR')))
    
    N_VEH = 1;
    ALPHAC_BAR = .5; % controls
    ALPHAM_BAR = .5;  % measurements
    ALPHAA_BAR = .5; % ACKs (if piggyback used, betaBar overrides gammaBar)
    
else
    
    N_VEH = 2;
    
    %ALPHAC_BAR = [.75;.75];
    %ALPHAM_BAR = [.75;.75];
    %ALPHAA_BAR = [.75;.75];
    
    ALPHAC_BAR = .75*ones(N_VEH,1);
    ALPHAM_BAR = .75*ones(N_VEH,1);
    ALPHAA_BAR = .75*ones(N_VEH,1);
    
end

%% system setup

% (system params are in setupSystemJLSPPC script)
setupSystemJLSPPC

% initial conditions
x_IC = 5*randn(size(A,1),1);
if(size(A,1)==2)
    % position only, no initial velocity (so like step resp from rest)
    x_IC(2)=0;x_IC(1)=5;
end

% (IF WANT TO DEBUG CONTROLLER - INIT ESTIMATOR PERFECTLY)
% xHat1 = xIC;P1 = 1*eye(2);

% (DEBUG SEQUENCES FOR ACKS)
%{
% tests w/ 1 then 2 then 3 missed ACKs
if(T_S==1)
    alpha_m(1:12) = [1 1 1 0 1 0 0 1 0 0 0 1];
elseif(T_S==2)
    alpha_m(1:22) = [0 1 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 1];
elseif(T_S==4)
    alpha_m(1:23) = [0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1];
end
%}

if(strfind(sched,'piggyback'))
    % ACK piggybacked to measurement
    alpha_a = alpha_m;   % overwrite
end

%% call sim fcn

[r] = simJLSPPC(SIM_LENGTH,N_HORIZON,A,Bu,Bw,C,Q,Qf,R,W,V,TAU_M,TAU_C,...
    TAU_A,TAU_AC,ALPHAC_BAR,PI_C,PI_M,PI_A,T_S,U_MAX,U_MIN,CODEBOOK,...
    X_MAX,X_MIN,x_IC,P_1,x_hat_1,w,v,alpha_c,alpha_m,alpha_a,...
    cov_prior_adj,N_ACKHISTORY,print_debug);
r.sys.sched = sched;
r.sys.system =system;
r.sys.ALPHAM_BAR = ALPHAM_BAR;
r.sys.ALPHAA_BAR = ALPHAA_BAR;

% (save r struct here if want)
%{
fname = sprintf('results_%s',dateString('DHM'));
save(fname,'r')

%or rename r
uniquer = r;
save(fname,'uniquer')

%}


%% plots

% load results struct named "r"

if(size(r.sys.C,1)==1 && size(r.sys.Bu,2)==1)
    
    plotJLSPPC_SISO(r)
    
else
    
    % very simple MIMO plot:
    
    NX_SYS = size(r.P,1);    % underlying system states (no buffer)
    SIM_LENGTH = size(r.X,2);
    
    CPlot = r.sys.C;
    figure
    subplot(3,1,[1 2])
    hx = plot(0:SIM_LENGTH-1,CPlot*r.X(1:NX_SYS,:));
    hold on
    hxh = plot(0:SIM_LENGTH-1,CPlot*r.Xh(1:NX_SYS,:),':');
    legend([hx(1) hxh(1)],'X','XHat')
    title('MIMO System (colors are i/o channels)')
    
    subplot(3,1,3)
    hu = stairs(repmat(0:SIM_LENGTH-1,[2,1])',r.u');
    xlabel('time step')
    ylabel('u')
    
end





