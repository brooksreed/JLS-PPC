% Script to directly test the MPC solver
% (useful for debugging slow or infeasible solver calls)

% (Load a set of inputs for schedMPC)
% (can use the snippet below while debugging to save inputs)

% Can compare performance of different CVX solvers on specific instances
solver = 'Mosek';
%solver = 'Gurobi';

print_debug_cvx = 1;

T_MPC = 40;

%% call the solver

tic
[U,cost,status,X_Out,~] = schedMPC(x_In,b_hat_MPC,p_i,T_MPC,A,Bu,M,E,...
    Q,Qf,R,umax,umin,xmin,xmax,u_deadband,solver,print_debug_cvx);
toc

%% plot the predicted plan and computed control trajectory 

N_VEH = size(Bu,2);

% printouts:
%full(U)
%X_Out

% (Use this with MIMO_RELATIVE_MEASUREMENTS)
%CPlot = kron(eye(N_VEH),[1 0]);

figure
subplot(3,1,[1 2])
plot((CPlot*X_Out)')
title('X Predicted')
grid on

subplot(3,1,3)
plot(full(U)')
ylabel('u')
grid on

xlabel('horizon (steps)')


break

%% snippet for saving debug inputs

% eg set breakpoint at simJLSPPC Line 431: 
% elseif( strcmp(status,'Failed') )

% and run this to save inputs to schedMPC for analysis with this script

x_In = XhMPC(1:NX_SYS,end,t+TAU_C)
b_hat_MPC = XhMPC((NX_SYS+1):end,end,t+TAU_C)
A=A_SYS
Bu = Bu_SYS
Q = QMPC
Qf = QfMPC
R = RMPC
umax = U_MAX
umin = U_MIN
xmax = [];
xmin = [];
u_deadband = [];

savename = 'test_inputs_MPC.mat';
save(savename,'x_In','b_hat_MPC','p_i','T_MPC','A','Bu','M','E',...
    'Q','Qf','R','umax','umin','xmax','xmin','u_deadband')