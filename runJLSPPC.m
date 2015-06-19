% Script to run a simple simulation of JLS-PPC
% double integrator or scalar systems
% a few options for schedules
% handles varying lengths of ACK histories

% calls simJLSPPC (which calls functions in core)
% some simple plots for SISO systems at end

% BR, 5/27/2015

% v1.0 6/13/2015
% v1.1 6/16/2015

% TO DO: 

% clean up/separate inputs vs. other? 
% one toggle for 'debugging mode'?
% more detailed plots for SISO
% more Pstars

% MIMO cleaner system setup, make Nc and Nm vs. Nv? 
% MIMO basic plots
% MIMO Pstar - ambiguity with partial ACKs

% automated saving/logging?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIM LENGTH
Ns = 30; % sim length

% SYSTEM (set up in setupSystemJPLPPC)
%system = 'SISO_DOUBLE_INTEGRATOR';
%system = 'MIMO_DOUBLE_INTEGRATOR';
system = 'SCALAR';

% SCHEDULE
% 'SISO_' options for _piggyback or _noACK:
% 'SISOALL' - [1],[1] - std. discrete time
% 'SISO2' - [1 0], [0 1] for pi, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi 

sched = 'SISO4_piggyback';
%sched = 'SISO4_noACK';

%sched = 'SISO2_noACK';
%sched = 'SISO2_piggyback';

%sched = 'SISO2ALLCONTROL_noACK';
%sched = 'SISO2ALLCONTROL_piggyback';
%sched = 'SISOALL_piggyback';
%sched = 'SISOALL_noACK';

% 'MIMO' options: MX, IL
%sched = 'MX_piggyback';
%sched = 'IL_piggyback';
%sched = 'MX_noACK';
%sched = 'IL_noACK';

% DELAYS
tc = 1; % control delay
tm = 1; % meas delay (Also ACK delay with piggyback unless overwritten)
ta = 1; % ACK delay

% ACK SETTINGS
% # ACK Histories sent (makes most sense to be multiple of schedule length)
% For 'SINGLE ACK': sys nACKHistory = Ts (schedule length)
nACKHistory = 5;
% adjustment to covariance priors due to no ACKs/control losses:
covPriorAdj = 1;

% MPC HORIZON:
NpMult = 4; % the MPC horizon Np = Ts*NpMult 

% packet success probabilities
if( ~isempty(strfind(system,'SISO')) || ~isempty(strfind(system,'SCALAR')))
    
    Nv = 1;
    alpha_cBar = .75; % controls
    alpha_mBar = .7;  % measurements
    alpha_aBar = .7; % ACKs (if piggyback used, betaBar overrides gammaBar)
    
else
   
    Nv = 2;

%     alpha_cBar = [.75;.75];
%     alpha_mBar = [.75;.75];
%     alpha_aBar = [.75;.75];

    alpha_cBar = 1*ones(Nv,1);
    alpha_mBar = 1*ones(Nv,1);
    alpha_aBar = 1*ones(Nv,1);

end

%% system setup
setupSystemJLSPPC

% initial conditions
xIC = 5*randn(size(A,1),1);
if(size(A,1)==2)
    % position only, no initial velocity (so like step resp from rest)
    xIC(2)=0;xIC(1)=5;
end

% (IF WANT TO DEBUG CONTROLLER - INIT ESTIMATOR PERFECTLY)
% xHat1 = xIC;P1 = 1*eye(2);

% hardcode ta
%ta = 2;disp('OVERWRITING ta')

% hardcoded sequences for consistent debugging with SISOALL
%alpha_c(:,1:11) = [1 1 0 0 1 1 0 1 0 0 1];
%alpha_m(:,1:15) =  [1 1 1 0 1 1 0 0 1 1 0 0 0 1 1];

% with SISO2
% Pi_c          =  [1 0 1 0 1 0 1 0 1 0 1 0];
% Pi_m          =  [0 1 0 1 0 1 0 1 0 1 0 1];
%alpha_m(:,1:15) =  [0 1 0 0 0 1 0 0 0 0 0 1 0 1 0];
                        % one missed, 2 missed

% with SISO4
% Pi_c          =  [1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0];
% Pi_m          =  [0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1];
alpha_m(:,1:18) =  [0 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 1];
%alpha_m(:,1:18) =  [0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 1];

                                    % one missed                     
if(strfind(sched,'piggyback'))
    % ACK piggybacked to measurement
    alpha_a = alpha_m;   % overwrite
end

%% call sim fcn

% hack to overwrite covPriorAdj and rerun with same everything else
%covPriorAdj = 0;

% (pull-out run-specific parameters for re-running?)

[r] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tac,...
    alpha_cBar,Pi_c,Pi_m,Pi_a,Ts,umax,umin,codebook,Xmax,Xmin,xIC,P1,...
    xHat1,w,v,alpha_c,alpha_m,alpha_a,covPriorAdj,nACKHistory);

r.sys.sched = sched;
r.sys.system =system;
r.sys.alpha_mBar = alpha_mBar;
r.sys.alpha_aBar = alpha_aBar;


% (save r struct here if want)
%{

old = cd('C:\Brooks\Dropbox\Research Dropbox\MATLAB Code\JLS-PPC local');
fname = sprintf('results_%s',dateString('DHM'));
save(fname,'r')

%or rename r 
uniquer = r;
save(fname,'uniquer')
cd(old)

%}


%% plots
if(size(C,1)==1 && size(Bu,2)==1)
    plotJLSPPC_SISO(r)
end




