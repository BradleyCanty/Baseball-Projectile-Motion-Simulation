# -*- coding: utf-8 -*-
"""
cd_vs_re_curve_fitting.py

Description:
    This script applies a piecewise linear curve fit to the Cd vs Re data of 
    a MLB baseball using linear least squares regression
    See 'cd_vs_re_curve_fitting_notebook' for details
    
Created by:
    Bradley Canty
    2024/01/27
    
"""

from matplotlib import pyplot as plt
import numpy as np
from constants import rho_air, mu_air, D

Re_data = np.array([8.367321716000000015e+04,
    8.587514393000000564e+04,
    9.027899745999999868e+04,
    9.468285099999999511e+04,
    9.688477775999999722e+04,
    9.798574115000000165e+04,
    1.034905580999999947e+05,
    1.045915215000000026e+05,
    1.056924848000000056e+05,
    1.078944116000000067e+05,
    1.100963383999999933e+05,
    1.122982651000000042e+05,
    1.167021187000000064e+05,
    1.211059722000000038e+05,
    1.244088623999999982e+05,
    1.266107890999999945e+05,
    1.343175327999999863e+05,
    1.365194596000000020e+05,
    1.376204230000000098e+05,
    1.387213863000000129e+05,
    1.409233130999999994e+05,
    1.431252398999999859e+05,
    1.431252398999999859e+05,
    1.475290934000000125e+05,
    1.497310201999999990e+05,
    1.519329469000000099e+05,
    1.530339102999999886e+05,
    1.563368005000000121e+05,
    1.574377638999999908e+05,
    1.618416173999999883e+05,
    1.651445076000000117e+05,
    1.706493244999999879e+05,
    1.739522146000000066e+05,
    1.805579948999999906e+05,
    1.849618484999999928e+05,
    1.882647386000000115e+05,
    1.915676288000000059e+05,
    1.937695554999999877e+05,
    2.003753358000000007e+05,
    2.025772625999999873e+05,
    2.113849696999999869e+05,
    2.135868963999999978e+05,
    2.223946034999999974e+05,
    2.312023105999999971e+05,
    2.323032740000000049e+05,
    2.389090542999999889e+05,
    2.422119444000000076e+05,
    2.444138711999999941e+05,
    2.488177246999999916e+05,
    2.532215781999999890e+05,
    2.576254317999999912e+05,
    2.642312121000000043e+05,
    2.708369923999999883e+05,
    2.730389192000000039e+05,
    2.774427726999999722e+05,
    2.818466261999999988e+05,
    2.928562601000000141e+05])

Cd_data = np.array([4.899999999999999911e-01,
    4.500000000000000111e-01,
    4.299999999999999933e-01,
    4.699999999999999734e-01,
    4.650000000000000244e-01,
    4.450000000000000067e-01,
    4.249999999999999889e-01,
    4.349999999999999978e-01,
    4.199999999999999845e-01,
    4.000000000000000222e-01,
    4.400000000000000022e-01,
    4.199999999999999845e-01,
    3.900000000000000133e-01,
    4.149999999999999800e-01,
    4.050000000000000266e-01,
    3.950000000000000178e-01,
    4.000000000000000222e-01,
    3.950000000000000178e-01,
    3.900000000000000133e-01,
    3.649999999999999911e-01,
    3.499999999999999778e-01,
    3.750000000000000000e-01,
    3.900000000000000133e-01,
    3.900000000000000133e-01,
    3.750000000000000000e-01,
    3.800000000000000044e-01,
    3.599999999999999867e-01,
    3.649999999999999911e-01,
    3.800000000000000044e-01,
    3.599999999999999867e-01,
    3.350000000000000200e-01,
    3.449999999999999734e-01,
    3.150000000000000022e-01,
    3.300000000000000155e-01,
    2.800000000000000266e-01,
    3.300000000000000155e-01,
    2.700000000000000178e-01,
    2.899999999999999800e-01,
    2.750000000000000222e-01,
    2.899999999999999800e-01,
    2.750000000000000222e-01,
    2.899999999999999800e-01,
    2.800000000000000266e-01,
    2.949999999999999845e-01,
    2.800000000000000266e-01,
    2.949999999999999845e-01,
    2.800000000000000266e-01,
    2.979999999999999871e-01,
    2.800000000000000266e-01,
    2.750000000000000222e-01,
    2.650000000000000133e-01,
    2.750000000000000222e-01,
    3.099999999999999978e-01,
    3.049999999999999933e-01,
    2.999999999999999889e-01,
    2.600000000000000089e-01,
    2.800000000000000266e-01])

'''
Apply linear least-squares fit to data. Create two lines: one corresponding to
data having Re < 202195, and another having Re > 202195. This corresponds to
Kensrud's approach. Stage all entries in the dataframe having Re < 202195 in a 
list, and stage all other entries in another list. Each list consists of two
lists: one corresponding to the Re column, and one corresponding to the Cd column
'''
Re1 = []
Cd1 = []

Re2 = []
Cd2 = []

for i in range(0,Re_data.shape[0]):
    if Re_data[i] < 200000:
        Re1.append(Re_data[i])
        Cd1.append(Cd_data[i])
    else:
        Re2.append(Re_data[i])
        Cd2.append(Cd_data[i])

'''
Apply least-squares linear regression on each list pair to find the slope and 
displacement of the line
x = (A^T*A)^(-1)*A^T*b = [a;Cd0]
where 
A = matrix of known values = [Re1,1;Re2,1;...;ReN,1]
b = solution vector = [Cd1;Cd2;...;CdN]
a = slope
Cd0 = initial displacement
'''

#First curve fit
A1 = np.matrix([Re1,np.ones(len(Re1))]).T
b1 = np.matrix(Cd1).T
x1 = np.dot(np.dot(np.linalg.inv(np.dot(A1.T,A1)),A1.T),b1) #1st entry is slope, 2nd entry is displacement
x1 = x1.A1 #A1 is shorthand attribute to convert matrix to flattened 1D ndarray

#Second curve fit
A2 = np.matrix([Re2,np.ones(len(Re2))]).T
b2 = np.matrix(Cd2).T
x2 = np.dot(np.dot(np.linalg.inv(np.dot(A2.T,A2)),A2.T),b2) #1st entry is slope, 2nd entry is displacement
x2 = x2.A1 #A1 is shorthand attribute to convert matrix to flattened 1D ndarray

#Solve for point where lines intersect
#i.e. find Re where x1[0]*Re + x1[1] == x2[0]*Re + x2[1]
ReEqual = float((x2[1]-x1[1]) / (x1[0] - x2[0]))

vmax_data = max(Re_data)*mu_air/(rho_air*D) #m/s, max velocity represented in data
vmin_data = min(Re_data)*mu_air/(rho_air*D) #m/s, min velocity represented in data

if __name__ == '__main__':
    #Create scatter plot of Cd vs Re
    plt.figure()
    plt.plot(Re_data,Cd_data,'b.')
    plt.title('Drag Coefficient vs Reynolds Number*')
    plt.figtext(0.5,0,"* 'Determining Aerodynamic Properties of Sports Balls In-Situ' by Jeff Kensrud", ha="center",fontsize=9)
    plt.subplots_adjust(bottom=0.15)
    plt.xlabel('Re')
    plt.ylabel('Cd')
    plt.grid(True)
    
    
    #Plot the 1st curve fit
    plt.plot([Re1[0],ReEqual],
             [float(Re1[0]*x1[0] + x1[1]),float(ReEqual*x1[0] + x1[1])],'c')
    
    #Plot the 2nd curve fit
    plt.plot([ReEqual,Re2[-1]],
             [float(ReEqual*x2[0] + x2[1]),float(Re2[-1]*x2[0] + x2[1])],'m')

















