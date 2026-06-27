# -*- coding: utf-8 -*-
"""
console_utilities.py
"""
import math
from cd_vs_re_curve_fitting import vmin_data,vmax_data

deg_symbol = u'\N{DEGREE SIGN}'

def get_user_input(debug=False):

    #Prompt the user for launch metrics
    print('Simulate the projectile motion of a baseball...')

    if debug:
        v0_imp = 100 #MLB homerun mean batted ball speed [1]
        v0 = v0_imp * 5280 * 12 * 2.54 / 100 / 3600 #Initial velocity [m/s]
        print(f'Input the batted ball speed [mi/hr]: {v0_imp}')
        alpha0_deg = 36
        print(f'Input the launch angle [{deg_symbol}]: {alpha0_deg}')
        time_step_size = .5
        print(f'Input the time step size [seconds] (0.1 to 1 suggested): {time_step_size}')
        
    else:
        while(True):
            v0_imp_str = input('Input the batted ball speed [mi/hr]: ') #Initial velocity [mi/hr]
            v0_imp = float(v0_imp_str)
            v0 = v0_imp * 5280 * 12 * 2.54 / 100 / 3600 #Initial velocity [m/s]
            if (v0 < vmin_data or v0 > vmax_data):
                vmin_data_imp = vmin_data/(5280 * 12 * 2.54 / 100 / 3600)
                vmax_data_imp = vmax_data/(5280 * 12 * 2.54 / 100 / 3600)
                print(f'ERROR: please enter a BBS between {vmin_data_imp:.1f} and {vmax_data_imp:.1f} [mi/hr]')
                continue
            else:
                break
        
        while (True):
            alpha0_deg_str = input(f'Input the launch angle [{deg_symbol}]: ') #Launch angle [deg]
            alpha0_deg = float(alpha0_deg_str)
            
            if (alpha0_deg < 0 or alpha0_deg > 90):
                print('ERROR: please input a launch angle between 0 and 90 degrees')
                continue
            else:
                break
        
        time_step_size = float(input('Input the time step size [seconds]: '))
            
    alpha0 = alpha0_deg*math.pi/180 #Launch angle [rad]
    
    return v0_imp,v0,alpha0_deg,alpha0,time_step_size

def output_data_to_console():
    return
