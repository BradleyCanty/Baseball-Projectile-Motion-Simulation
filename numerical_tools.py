# -*- coding: utf-8 -*-
"""
numerical_tools.py
"""

import math
import numpy as np
import ode_tools as ot
from cd_vs_re_curve_fitting import ReEqual, x1, x2
from constants import g, rho_air, mu_air, r, D, m

def get_reynolds_number(v):
    return rho_air * v * D / mu_air #Reynolds number of projectile

def get_drag_coeff(Re):
    """
    Form the drag coefficient functions which correspond to the baseball's Reynolds Number (Re)
    Reference: 'Determining Aerodynamic Properties of Sports Balls In-Situ' by Jeff Kensrud
    The mapping of Re to Cd is modeled using a piecewise continuous fit consisting of two lines.
    See 'cd_vs_re_curve_fitting_notebook' Jupyter notebook for details.
    """
    #Line 1
    Cd_Re_slope_1 = x1[0]
    Cd0_1 = x1[1]
    Cd1 = lambda Re: Cd_Re_slope_1 * Re + Cd0_1 #Drag coefficient curve when Re < 200000

    #Line2
    Cd_Re_slope_2 = x2[0]
    Cd0_2 = x2[1]
    Cd2 = lambda Re: Cd_Re_slope_2 * Re + Cd0_2 #Drag coefficient curve when Re >= 200000
    
    if Re < ReEqual:
        Cd = Cd1(Re)
    else:
        Cd = Cd2(Re)
    return Cd

def sys_of_odes(t,state):
    '''
    Establish system of ODEs to solve, i.e., the time derivative of the
    state vector
    
    Parameters:
    -----------
    t : float array
        array of times at which to evaluate the function
    state:
        state vector given by (x, y, vx, vy)
    
    Returns:
    -------
    float array
        state vector derivative given by (vx,vy,ax,ay)  
    '''
    x,y,vx,vy = state
    
    v = math.sqrt(vx ** 2 + vy ** 2) #velocity magnitude
    Re = get_reynolds_number(v)
    Cd = get_drag_coeff(Re)
    A = math.pi * r ** 2 #projected frontal area of projectile
    k_Newton = .5 * rho_air * Cd * A #Newton drag parameter
    zeta = k_Newton / m #specific Newton drag parameter
    
    ax = -zeta * v * vx
    ay = -g - zeta * v * vy
    
    state_dot = np.array([vx, vy, ax, ay])
    
    return state_dot

def propagate_odes(state0,tspan,dt,method='rk4'):
    ode         = sys_of_odes
    func        = ot.methods[method]
    times       = np.arange(0,tspan+dt/2,dt) #dt/2 enables including the end point in the array
    steps       = len(times)
    states      = np.zeros((steps, len(state0)))
    states[0]   = state0
    
    for step in range(steps - 1):
        states[step + 1] = func(ode,times[step], states[step], dt)
        rx = states[:,0]
        ry = states[:,1]
        vx = states[:,2]
        vy = states[:,3]
        
    return times, rx, ry, vx, vy

def secant_method(f,x0,x1,tol,max_iter = 10):
    '''
    Used when function f is a continous function
    
    Inputs:
        f: function which takes x as input
        x0: first starting guess
        x1: second starting guess
        tol: error tolerance, such that f(x_soln) < tol
    Outputs:
        the root of function f
    '''
    i = 0
    x_prev = x0
    x_curr = x1
    x_next = 0

    while i < max_iter:

        x_next = x_curr - f(x_curr) * (x_curr - x_prev)/(f(x_curr)-f(x_prev))
        
        if abs(f(x_next)) < tol:
            #print(f'returned root in {i} iterations')
            return x_next
        
        x_prev = x_curr
        x_curr = x_next
        i=i+1
        
        if i > max_iter:
            print('WARNING: number of iterations of secant method exceeded max_iter')
            
def secant_method_for_sequence(f,x0,x1,tol,step_size,max_iter = 10):
    '''
    This is used when function f produces a sequence of values, i.e., f is not a continous function
    
    Inputs:
        f: function which takes x as input
        x0: first starting guess
        x1: second starting guess
        tol: error tolerance, such that f(x_soln) < tol
    Returns:
        the root of function f
    '''
    if x0 == 0:
        x0 = step_size
    if x1 == 0:
        x1 = step_size
    
    i = 0
    x_prev = x0
    x_curr = x1
    x_next = 0

    while i < max_iter:
        if f(x_curr,step_size)-f(x_prev,step_size) == 0:
            print('WARNING: division by zero.')
            print(f'Consider decreasing step size to less than |x_curr-x_prev| = {abs(x_curr-x_prev)} or decreasing tolerance')
            print(f'At i = {i}, f(x_next,step_size) = {f(x_next,step_size)} while tol = {tol}')
            #print(f'Number of points in time sequence is {x_curr/step_size+1}')
            break
        x_next = x_curr - f(x_curr,step_size) * (x_curr - x_prev)/(f(x_curr,step_size)-f(x_prev,step_size))
        
        if abs(f(x_next,step_size)) < tol:
            #print(f'returned root in {i} iterations')
            return x_next
        
        x_prev = x_curr
        x_curr = x_next
        i=i+1
        
        if i > max_iter:
            print('WARNING: number of iterations of secant method exceeded max_iter')
    
if __name__ == '__main__':
    
    f=lambda x: 9 - x**2 
    x_root = secant_method(f,1,2,0.001)
    
    
    #Inital conditions
    v0 = 44.7
    alpha0 = 0.6283
    x0 = 0.0
    y0 = 1.0

    #Initial velocity components
    vx0 = v0 * math.cos(alpha0)
    vy0 = v0 * math.sin(alpha0)
    
    v = lambda t,dt: propagate_odes((x0,y0,vx0,vy0), t, dt, 'rk4')[4][-1]
    t_peak = secant_method_for_sequence(v,1,2,tol=0.005,step_size=0.001)
    









