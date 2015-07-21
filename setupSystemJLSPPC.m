% sets up JLSPPC system for simulations
% 'inputs' are system, sched
% system params are set inside here

% pulls these settings outside of run script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM INIT/SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE -- following notation and descriptions in paper, 
% N_VEH sets the number of comms channels 
% (even though systems may not represent vehicles)
% must have equal number comms channels for controls and measurements
% (each "agent" has an input channel and an output channel)

switch system
    
    case 'SCALAR'
                
        % controller settings        
        U_MAX = 1;
        U_MIN = -1;
        N_QUANT_LEVELS = 33;   % quantization levels
        Q = 10;
        Qf = 10*Q;
        R = 1;
        
        % process/measurement noise
        W_gen = 1;  % cov for process noise input with Bw
        V = 1;
        
        % initial covariance
        P_1 = 9;
        
        % Scalar integrator
        A = 1;
        Bu = 1;
        Bw = Bu;
        C = 1;      % position output
        X_MAX = [];X_MIN = [];
        CODEBOOK = linspace(U_MIN,U_MAX,N_QUANT_LEVELS);
        
 case 'SISO_DOUBLE_INTEGRATOR'
        
     % 1 input, 1 output

        % controller settings
        U_MAX = 10;      
        U_MIN = -10;
        N_QUANT_LEVELS = 15;   % quantization levels
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = 1;
        
        % process/measurement noise
        W_gen = 1;  % cov for process noise input with Bw
        V = 1;
                
        % cov... uncertain position but better-known velocity (closer to zero)
        P_1 = [25,0;0,9];     
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = [0.5;1];
        Bw = Bu;
        C = [1 0];      % position output
        X_MAX = [];X_MIN = [];
        CODEBOOK = linspace(U_MIN,U_MAX,N_QUANT_LEVELS);
        
        
 case 'MIMO_DOUBLE_INTEGRATOR'
     
        % 2 inputs, 2 outputs
        
        % controller settings
        U_MAX = [1;10];
        U_MIN = [-1;-10];
        N_QUANT_LEVELS = 15;   % quantization levels
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = eye(2);
               
        % process/measurement noise
        W_gen = eye(2); % cov for process noise input with Bw
        V = eye(2);
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P_1 = [25,0;0,9];  
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = eye(2);
        Bw = Bu;
        C = eye(2);    % full state output
        X_MAX = [];X_MIN = [];
        CODEBOOK = linspace(min(U_MIN),max(U_MAX),N_QUANT_LEVELS);
        
    case 'MIMO_RELATIVE_MEASUREMENTS'

        % controller settings
        U_MAX = 10*ones(N_VEH,1);
        U_MIN = -10*ones(N_VEH,1);
        N_QUANT_LEVELS = 9;   % quantization levels
        Q1 = [10,0;0,1];
        Q = kron(eye(N_VEH),Q1);
        Qf = 10*Q;
        R = eye(N_VEH);
               
        % process/measurement noise
        W_gen = .1*eye(N_VEH); % cov for process noise input with Bw
        V = 1*eye(N_VEH);
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P_11 = [25,0;0,9];  
        P_1 = kron(eye(N_VEH),P_11);
        
        % A, B: set of double integrators
        % Double Integrator
        A1 =[1,1;0,1];
        Bu1 = [0.5;1];
        
        A = kron(eye(N_VEH),A1);
        Bu = kron(eye(N_VEH),Bu1);
        Bw = Bu;
        
        % C: relative position measurements
        C = zeros(N_VEH,size(A,1));
        C(1,1) = 1;
        for i = 2:N_VEH
            C(i,(2*i-3)) = -1;
            C(i,(2*i-1)) = 1;
        end
        
        X_MAX = [];X_MIN = [];
        CODEBOOK = linspace(min(U_MIN),max(U_MAX),N_QUANT_LEVELS);
          
end

% make W for KF based on W_gen (KF uses APA' + W form)
W = Bw*W_gen*Bw';

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
% note -- use W_gen here (Bw*w is the process noise input)
w = sqrt(W_gen)*randn(size(Bw,2),SIM_LENGTH);
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


