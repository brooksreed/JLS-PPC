# JLS-PPC
Jump Linear System framework for networked control with schedules, delays, packet loss

Author: Brooks Reed, brooksr8@gmail.com
Release 1.0 7/21/2015
Release 1.1 7/28/2015 - minor updates to readme, added testMPC script

Link to paper (describes notation, JLS system, algorithms):(coming soon)

runJLSPPCExamples is the main file for user input -- runs a simulation

User input: 
- underlying system
    - (some options for these are provided in setupSystemJLSPPC)
        - Scalar integrator system
        - SISO double integrator (eg mass with position measurements and force inputs)
        - MIMO double integrator (eg mass with position and velocity measurements, and force and velocity inputs)
        - Chain of double integrators with relative position measurements
    - system parameters (set manually or using setupSystemJLSPPC)
        - A, Bu, Bw, C
        - process noise W, measurement noise V
        - x_IC initial condition
- MPC parameters and settings
    - control constraints U_MAX, U_MIN
    - control quantization (linear, with N_QUANT_LEVELS)
    - (optional) state constraints X_MAX, X_MIN
    - MPC horizon N_HORIZON
    - MPC weights on states Q, control R, and terminal state Qf
- schedule for control, measurement, and control ACK packets
    - (some options constructed in createSchedule)
    - PI_C for controls, PI_M for measurements, PI_A for ACKs
- packet success probabilities (bernoulli) 
    - ALPHAC_BAR, ALPHAM_BAR, ALPHAA_BAR 
- communication delays (integer steps)
    - TAU_C, TAU_M, TAU_A 
- estimation initialization 
    - P_1 initial covariance
    - x_hat_1 initial estimate
- (experimental): option to turn on/off covariance prior adjustments in KF    

simJLSPPC runs a simulation using all elements of the framework
- outputs a struct of results

plotJLSPPC_SISO makes some simple plots for SISO systems
- runs after a sim in the example script, can also be run after loading a saved results struct

/core directory contains the core JLS functions 
- computePstars: utility function for generating covariance prior adjustment (P*) coefficients
- createSchedule: constructs time series scheduling variables based on a few prototype schedule options
- evaluatePstars: fast lookup table of analytic P* coefficient formulas
- JLSJumpEstimator: handles delayed/lossy ACKs and updates jump variables accordingly
- JLSKF: Kalman Filter that handles missed measurements
    - state priors for control use the best information for the jump variable (depends on ACKs received)
    - option for adjusting the covariance prior when control packets are uncertain
- makeDc, makeDm: utilities for generating jump matrices
- makeM: utility for generating buffer shift matrix
- paramsNow: grabs a receding-horizon window of constraints for use by the MPC
- prepMPC: computes state estimate for the time control is scheduled to be applied via a forward propagation to bridge the round-trip communication delay
- schedMPC: solves the MPC optimization with awareness of control packet scheduling 
    - The MPC solver requires CVX (http://cvxr.com/cvx/) 
    - Currently Mosek (available with CVX academic license) is recommended, but any QP solver hooked into CVX should work

testMPC is a debugging tool for running a single solver instance of MPC

Notes: 
- More detailed comments in code
- "help <functionname>" should help with i/o
- Variable names in code closely match the paper in most cases.  
- Constants are in CAPS. 



