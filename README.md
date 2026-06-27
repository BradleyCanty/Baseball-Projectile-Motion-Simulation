# Baseball-Projectile-Motion-Simulation
<img width="3600" height="1957" alt="1_Baseball_Traj" src="https://github.com/user-attachments/assets/4b96070a-0802-4f71-aba7-be6dd3bb0534" />



This simulates the motion of a baseball when...
- subjected to vacuum
- subjected to Newton drag

The user inputs...
- the Batted Ball Speed
- the launch angle
- the time step used in computing the numerical solution


The variation of drag coefficient with Reynolds number is found using a piecewise continous curve fit to baseball wind tunnel data, derived using linear-least squares regression (see 'cd_vs_re_curve_fitting_notebook.ipynb').
<img width="3600" height="2400" alt="5_Baseball_Cd_vs_Re_Data" src="https://github.com/user-attachments/assets/e552edcf-be74-4ff2-80ac-b977e3853af8" />

Various numerical integration methods are used for computing the Newton drag trajectory:
- Euler method
- Tangential method
- 4th order Runge-Kutta method


Additionally, the secant method has been implemented for root finding, used in finding the peak and impact times.


The main program is 'projectile_motion_newton_drag.py'.

**An overview of this program, its outputs, and an in-depth account of the theory of projectile motion with and without air drag is found in the attached powerpoint (see 'BaseballTrajectorySimulation.pptx')**
