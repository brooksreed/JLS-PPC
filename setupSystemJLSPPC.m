% sets up JLSPPC system for simulations
% 'inputs' are system, sched
% system params are set inside here

% pulls these settings outside of run script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM INIT/SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch system
    case 'SISO_DOUBLE_INTEGRATOR'
        
        % controller settings
        U_MAX = 10;      
        U_MIN = -10;
        N_QUANT_LEVELS = 15;   % quantization levels
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = 1;
        
        % process/measurement noise
        W = 1;
        V = 1;
                
        % cov... uncertain position but better-known velocity (closer to zero)
        P_1 = [25,0;0,9];     
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = [0.5;1];
        Bw = Bu;
        C = [1 0];      % position output
        Xmax = [];Xmin = [];
        CODEBOOK = linspace(U_MIN,U_MAX,N_QUANT_LEVELS);
        
 case 'MIMO_DOUBLE_INTEGRATOR'
     
        % controller settings
        U_MAX = [10;10];
        U_MIN = [-10;-10];
        N_QUANT_LEVELS = 15;   % quantization levels
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = eye(2);
               
        % process/measurement noise
        W = eye(2);
        V = eye(2);
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P_1 = [25,0;0,9];  
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = eye(2);
        Bw = Bu;
        C = eye(2);    % full state output
        Xmax = [];Xmin = [];
        CODEBOOK = linspace(min(U_MIN),max(U_MAX),N_QUANT_LEVELS);
        
    case 'SCALAR'
        
        % controller settings        
        U_MAX = 1;
        U_MIN = -1;
        N_QUANT_LEVELS = 33;   % quantization levels
        Q = 10;
        Qf = 10*Q;
        R = 1;
        
        % process/measurement noise
        W = 1;
        V = 1;
        
        % initial covariance
        P_1 = 9;
        
        % Scalar integrator
        A =1;
        Bu = 1;
        Bw = Bu;
        C = 1;      % position output
        Xmax = [];Xmin = [];
        CODEBOOK = linspace(U_MIN,U_MAX,N_QUANT_LEVELS);
        
end

% check
if( length(ALPHAC_BAR)~=N_VEH )
    disp('WARNING: alpha_c, Nv mismatch')
end
if( length(ALPHAM_BAR)~=N_VEH )
    disp('WARNING: alpha_m, Nv mismatch')
end
if( rem(size(C,1),N_VEH) )
    disp('WARNING: Nv not into C')
end
if( rem(size(Bu,2),N_VEH) )
    disp('WARNING: Nv not into Bu')
end    
    
% estimation init
x_hat_1 = zeros(size(A,1),1);

% random time-series realizations
w = sqrt(W)*randn(size(Bw,2),SIM_LENGTH);
v = sqrt(V)*randn(size(C,1),SIM_LENGTH);

% packet loss sequences
alpha_m = zeros(N_VEH,SIM_LENGTH);
alpha_c = zeros(N_VEH,SIM_LENGTH);
alpha_a = zeros(N_VEH,SIM_LENGTH);
for k = 1:SIM_LENGTH
    alpha_m(:,k) = (sign(rand(N_VEH,1) - (1-(ALPHAM_BAR)))*0.5 + 0.5);
    alpha_c(:,k) = (sign(rand(N_VEH,1) - (1-(ALPHAC_BAR)))*0.5 + 0.5);
    alpha_a(:,k) = (sign(rand(N_VEH,1) - (1-(ALPHAA_BAR)))*0.5 + 0.5);
end

%%%%%%%%%%%%%%%%%%

if(strfind(sched,'piggyback'))
    if(TAU_A~=TAU_M)
        disp('warning - resetting ta to tm')
        TAU_A = TAU_M;
    end
    ALPHAA_BAR = ALPHAM_BAR;
elseif(strfind(sched,'NoACK'))
    TAU_A = 0;
end

% schedule time series
[PI_C,PI_M,PI_A,TAU_AC,T_S] = createSchedule(sched,N_VEH,SIM_LENGTH,TAU_C);
N_HORIZON = N_HORIZON_MULT*T_S;

