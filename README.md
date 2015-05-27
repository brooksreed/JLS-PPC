# JLS-MPC
Jump Linear System framework for networked control with schedules, delays, packet loss

simJLSPPC runs a simulation using all elements of the framework
simNoComms runs a simulation using deterministic scheduling only (no packet loss or delays)

core contains the core JLS functions (modified KF, jump estimator, schedule setup, MPC solver, etc.)
The MPC solver requires CVX and is currently set up to use gurobi (but any QP solver hooked into CVX will do)

experiments-kayaks-moos contains scripts to use the JLS estimation and control in real-time within the MOOS-IvP framework as implemented on the MIT Hovergroup autonomous kayaks.  (note - this code may need updating as the JLS core changes to add new features)



