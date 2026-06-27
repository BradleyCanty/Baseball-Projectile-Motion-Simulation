# Baseball-Projectile-Motion-Simulation
<img width="3600" height="2400" alt="1_Baseball_Traj" src="https://github.com/user-attachments/assets/6eaa0f51-27d5-4c6a-9cab-ffa74e05d33b" />


This simulates the motion of a baseball when...
- subjected to vacuum
- subjected to Newton drag

Drag is simulated using a piecewise continous curve fit to baseball wind tunnel data, derived using linear-least squares regression (see 'cd_vs_re_curve_fitting_notebook.ipynb').


Various numerical integration methods are used for computing the Newton drag trajectory:
- Euler method
- Tangential method
- 4th order Runge-Kutta method


Additionally, the secant method has been implemented for root finding, used in finding the peak and impact times.
The main program is 'projectile_motion_newton_drag.py'.

**An overview of this program, its outputs, and an in-depth account of the theory of projectile motion with and without air drag is found in the attached powerpoint (see 'BaseballTrajectorySimulation.pptx')**
