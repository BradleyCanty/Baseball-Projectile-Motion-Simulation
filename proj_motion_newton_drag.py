# -*- coding: utf-8 -*-
""" 
proj_motion_newton_drag.py

Description:
    Simulation of projectile motion with air drag (i.e. "Newton drag"), where 
    modelled projectile is a baseball. Plots are output into an 'Images' folder
    in the current working directory.
    This script uses the following libraries:
        - 'math' 
        - 'numpy' for arrays 
        - 'os' for finding the local directory and forming path strings
        - 'matplotlib' for plotting

Created by:
    Bradley Canty
    2026/06/25

References:
    [1] http://baseball.physics.illinois.edu/PointC.html
    [2] http://baseball.physics.illinois.edu/DragTPTMay2014.pdf
    [3] https://en.wikipedia.org/wiki/Projectile_motion
    [4] https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html#scipy.integrate.solve_ivp
    [5] https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html
    [6] https://baseballaero.com/2019/05/03/baseball-drag-crisis/
    [7] http://baseball.physics.illinois.edu/KensrudThesis.pdf
    

"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
from constants import g
from numerical_tools import get_drag_coeff, get_reynolds_number, propagate_odes, secant_method_for_sequence
from cd_vs_re_curve_fitting import Re_data, Cd_data,Re1,Re2,ReEqual,x1,x2
from console_utilities import get_user_input,deg_symbol


#Inital conditions
x0 = 0.0
y0 = 1.0
v0_imp,v0,alpha0_deg,alpha0,time_step_size = get_user_input(debug=False)

#Initial velocity components
vx0 = v0 * math.cos(alpha0)
vy0 = v0 * math.sin(alpha0)

#Solve for position and velocity of projectile in a vacuum
t_vac_impact = (vy0 + math.sqrt(vy0**2 + 2*g*y0))/g
t_vac = np.linspace(0,t_vac_impact,101)
x_vac = vx0*t_vac + x0
y_vac = -0.5*g*t_vac**2 + vy0*t_vac + y0
vx_vac = vx0 * np.ones(len(t_vac))
vy_vac = -g*t_vac + vy0
y_vac_peak = vy0**2/(2*g) + y0

#Solve for peak time (i.e. find time when v_y = 0) when baseball subjected to air drag
t_peak_guess1 = 0;
t_peak_guess2 = 1;
t_peak = secant_method_for_sequence(lambda t,dt: 
                                    propagate_odes((x0,y0,vx0,vy0), t, dt, method='rk4')[4][-1],
                                    t_peak_guess1,
                                    t_peak_guess2,
                                    tol=0.005,
                                    step_size=0.0005)

#Solve for peak height (i.e., find max projectile altitude) when baseball subjected to air drag
y_peak = propagate_odes((x0,y0,vx0,vy0), tspan=t_peak, dt=0.001, method='rk4')[2][-1]

#Solve for impact time (i.e., find time when y = 0) using Newton-Raphson root solving method
t_impact_guess1 = 2 * t_peak
t_impact_guess2 = 2 * t_peak + 1
t_impact = secant_method_for_sequence(lambda t,dt: 
                                      propagate_odes((x0,y0,vx0,vy0), t, dt, method='rk4')[2][-1],
                                      t_impact_guess1,
                                      t_impact_guess2,
                                      tol=0.005,
                                      step_size=0.0005)

tspan = t_impact
dt = time_step_size
t_euler, rx_euler, ry_euler, vx_euler, vy_euler = propagate_odes((x0,y0,vx0,vy0), tspan, dt, 'euler')
t_mid, rx_mid,ry_mid,vx_mid,vy_mid = propagate_odes((x0,y0,vx0,vy0), tspan, dt, 'midpoint')
t_rk4,rx_rk4,ry_rk4,vx_rk4,vy_rk4 = propagate_odes((x0,y0,vx0,vy0), tspan, dt, 'rk4')

#Solve for Reynolds number of baseball throughout trajectory
Re = np.zeros(len(t_rk4))
for i in range(0,len(t_rk4)):
    v = math.sqrt(vx_rk4[i]**2 + vy_rk4[i]**2)
    Re[i] = get_reynolds_number(v)
Re_max = max(Re)
Re_min = min(Re)

#Display inputs and outputs in the console
print('\nINPUTS:')
print(f'Batted ball speed: {v0_imp:.1f} [mi/hr] = {v0:.1f} [m/s]')
print(f'Launch angle: {alpha0_deg:.1f} [{deg_symbol}]')
print(f'Time step size [seconds]: {time_step_size}')
print('\nOUTPUTS:') #where BBS = 100 [mi/hr] and launch angle = 38 [deg]
print('Vacuum trajectory:')
print(f'\tTime of flight: {t_vac_impact:.2f} [s]') # returns 5.65 [s]
print(f'\tHorizontal range: {x_vac[-1]:.1f} [m]')  # returns 198.9 [m]
print(f'\tMaximum height: {y_vac_peak:.1f} [m]')   # returns 38.6 [m]
print('Newton drag trajectory:')
print(f'\tTime of flight (RK4): {t_impact:.2f} [s]')  # returns 4.83 [s] 
print(f'\tHorizontal range (RK4): {rx_rk4[-1]:.1f} [m]')   # returns 111.3 [m]
print(f'\tMaximum height (RK4): {y_peak:.1f} [m]')    # returns 29.1 [m]
print(f'\tReynolds number range (RK4): [{Re_min:.1f},{Re_max:.1f}]')
img_dir = os.getcwd() + '\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)
print('\nImage output location: ', img_dir)

# Plot trajectory
fig, ax = plt.subplots()
fig.canvas.draw()

plt.plot(rx_euler, ry_euler, 'r.', label = 'Newton drag, Euler method')
plt.plot(rx_mid, ry_mid, 'g^', label = 'Newton drag, midpoint method')
plt.plot(rx_rk4, ry_rk4, 'bv', label = 'Newton drag, 4th order Runge-Kutta')
plt.plot(x_vac,y_vac,'m-', label = 'vacuum')

plt.title(f'BBS = {v0_imp:.1f} [mi/hr], launch angle = {alpha0_deg:.1f}{deg_symbol}',color='gray')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

plt.grid(True)
plt.gca().set_aspect('equal')

top_ax_y = ax.get_position().y1
bottom_ax_y = ax.get_position().y0
ax_y_delta = top_ax_y - bottom_ax_y
plt.suptitle('Baseball Trajectory',fontsize=16, y = top_ax_y + .1 + .05 * ax_y_delta, verticalalignment = 'bottom')
plt.legend(loc="upper center", bbox_to_anchor=(0.5, bottom_ax_y - 1.2 + ax_y_delta), ncol=1)

plt.savefig(img_dir + '/1_Baseball_Traj.png', dpi = 600)
plt.tight_layout()
plt.show()

#Plot horizontal velocity vs time
plt.figure()
plt.plot(t_euler, vx_euler, 'r.', label = 'Newton drag, Euler method')
plt.plot(t_mid, vx_mid, 'g^', label = 'Newton drag, midpoint method')
plt.plot(t_rk4, vx_rk4, 'bv', label = 'Newton drag, 4th order Runge-Kutta')
plt.plot(t_vac, vx_vac, 'm-', label = 'vacuum')
plt.suptitle('Baseball\'s Horizontal Velocity Component',fontsize=16,y=1)
plt.title(f'BBS = {v0_imp:.1f} [mi/hr], launch angle = {alpha0_deg:.1f}{deg_symbol}',color='gray')
plt.xlabel('t [s]')
plt.ylabel('$v_x$ [m/s]')
plt.legend(loc="upper right") 
plt.grid(True)
plt.savefig(img_dir + '/2_Baseball_Horiz_Vel.png', dpi = 600)

#Plot vertical velocity vs time
plt.figure()
plt.plot(t_euler, vy_euler, 'r.', label = 'Newton drag, Euler method')
plt.plot(t_mid, vy_mid, 'g^', label = 'Newton drag, midpoint method')
plt.plot(t_rk4, vy_rk4, 'bv', label = 'Newton drag, 4th order Runge-Kutta')
plt.plot(t_vac, vy_vac, 'm-', label = 'vacuum')
plt.suptitle('Baseball\'s Vertical Velocity Component',fontsize=16,y=1)
plt.title(f'BBS = {v0_imp:.1f} [mi/hr], launch angle = {alpha0_deg:.1f}{deg_symbol}',color='gray')
plt.xlabel('t [s]')
plt.ylabel('$v_y$ [m/s]')
plt.legend(loc="upper right") 
plt.grid(True)
plt.savefig(img_dir + '/3_Baseball_Vert_Vel.png', dpi = 600)

#Plot Reynolds number of baseball vs time
plt.figure()
plt.plot(t_rk4,Re,c='g')
plt.suptitle('Baseball\'s Reynolds Number vs Time',fontsize=16,y=1)
plt.title(f'BBS = {v0_imp:.1f} [mi/hr], launch angle = {alpha0_deg:.1f}{deg_symbol}',color='gray')
plt.xlabel('time [s]')
plt.ylabel('Reynolds number')
plt.grid(True)

#Place text box in upper right corner
textstr = r"$Re =  \frac{\rho \, v \, D}{\mu}$" "\n" "where \n" r"$\rho = air\ density$" "\n" "$v=velocity$" "\n" r"$D=baseball\ diameter$" "\n" "$\mu=dynamic\ viscosity$"
props = dict(boxstyle='round',facecolor = (.96, .87, .7))
plt.text(0.4, 0.9, textstr, transform= plt.gca().transAxes, fontsize = 12,
        verticalalignment = 'top', bbox = props)

plt.savefig(img_dir + '/4_Baseball_Reynolds_Num.png', dpi = 600)

#Plot the Cd vs Re data and the piecewise linear curve fits
plt.figure()
plt.plot(Re_data,Cd_data,'b.')
plt.plot([Re1[0],ReEqual],[float(Re1[0]*x1[0] + x1[1]),float(ReEqual*x1[0] + x1[1])],'c')   #1st curve fit
plt.plot([ReEqual,Re2[-1]],[float(ReEqual*x2[0] + x2[1]),float(Re2[-1]*x2[0] + x2[1])],'m') #2nd curve fit
plt.title('Baseball\'s Empirically-derived $C_D$ Curve *')
plt.subplots_adjust(bottom=0.15)
plt.figtext(0.5,0,"* 'Determining Aerodynamic Properties of Sports Balls In-Situ' by Jeff Kensrud", ha="center",fontsize=9)
plt.xlabel('Re')
plt.ylabel('Cd')
plt.grid(True)

#Place text box in upper right corner
textstr = r"$Re =  \frac{\rho \, v \, D}{\mu}$" "\n" "where \n" r"$\rho = air\ density$" "\n" "$v=velocity$" "\n" r"$D=baseball\ diameter$" "\n" "$\mu=dynamic\ viscosity$"
props = dict(boxstyle='round',facecolor = (.96, .87, .7))
plt.text(0.5, 0.9, textstr, transform= plt.gca().transAxes, fontsize = 12,
        verticalalignment = 'top', bbox = props)

plt.savefig(img_dir + '/5_Baseball_Cd_vs_Re_Data.png', dpi = 600)

#Plot drag coefficient range used in the simulation
Re_range = np.linspace(Re_data[0],Re_data[-1],201)
Cd_range = np.zeros(len(Re_range))
for i in range(0,len(Re_range)):
    Cd_range[i] = get_drag_coeff(Re_range[i])
        
Re_observed = np.linspace(Re_min,Re_max,201)
Cd_observed = np.zeros(len(Re_observed))
for i in range(0,len(Re_observed)):
    Cd_observed[i] = get_drag_coeff(Re_observed[i])
        
plt.figure()
plt.plot(Re_range,Cd_range,'b-')
plt.plot(Re_observed,Cd_observed,'r-', label = 'used in simulation')
plt.suptitle('Drag Coefficient vs Reynolds Number*',fontsize=16,y=1.005)
plt.title(f'BBS = {v0_imp:.1f} [mi/hr], launch angle = {alpha0_deg:.1f}{deg_symbol}',color='gray')
plt.xlabel('Re')
plt.ylabel('$C_D$')
plt.subplots_adjust(bottom=0.15)
plt.figtext(0.5,0,"* 'Determining Aerodynamic Properties of Sports Balls In-Situ' by Jeff Kensrud", ha="center",fontsize=9)
plt.legend()
plt.grid(True)
plt.xlim([min(min(Re_range),min(Re_observed)),max(max(Re_range),max(Re_observed))])
plt.ylim([min(min(Cd_data),min(Cd_observed)),max(max(Cd_data),max(Cd_observed))])

#Place text box in upper right corner
textstr = r"$Re =  \frac{\rho \, v \, D}{\mu}$" "\n" "where \n" r"$\rho = air\ density$" "\n" "$v=velocity$" "\n" r"$D=baseball\ diameter$" "\n" "$\mu=dynamic\ viscosity$"
props = dict(boxstyle='round',facecolor = (.96, .87, .7))
plt.text(0.5, 0.9, textstr, transform= plt.gca().transAxes, fontsize = 12,
        verticalalignment = 'top', bbox = props)

plt.savefig(img_dir + '/6_Baseball_Cd_vs_Re_in_Sim.png', dpi = 600)


