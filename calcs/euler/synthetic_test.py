"""
Synthetic test 1

A Python program to compute the Synthetic test 1 
Distinct SIs and strong nonlinear magnetic base level 

This code is released from the paper: 
Reliable Euler deconvolution estimates throughout the
vertical derivatives of the total-field anomaly

The program is under the conditions terms in the file README.txt

authors: Felipe F. Melo and Valeria C.F. Barbosa, 2019
email: felipe146@hotmail.com, valcris@on.br
"""

"""
Input:
    
input/synthetic_data.dat - 2d-array with "n" rows by 4 columns: 
    x-coordinate, y-coordinate, z-coordinate, anomaly. Where "n" rows
    correspond to the size of the data.

Parameters:

Size of the moving data window:    
    winsize - an odd integer number. 
              Ex.: for a moving data window of 5 x 5 grid points -> winsize = 5
                                  
Percentage of the solutions that will be keep:
    filt - a float number ranging from 0.0 to 1.0. 
           Ex.: to keep 10% of the solutions -> filt = 0.1

Structural indices used:
    SI_vet - an array that can store any of the four SIs.
             Ex.: to test only the SI = 1 -> SI_vet = [1]
                  to test the four SIs -> SI_vet = [0.01,1,2,3]

The areas to compute the statistics about the mean of the northing,
easting and depth estimates:
    area_cla  - array defining the four vertices of a polygon 
                [south,north,west,east]
"""

import numpy as np
import plot_functions as plt_fc
import euler_python as euler
import estimates_statistics as est_stats

# Input data
data_input=np.loadtxt('input/synthetic_data.dat')
shape = (120, 140)
area = [0, 24000, 0, 28000]
xi=data_input[:,0]
yi=data_input[:,1]
zi=data_input[:,2]
data=data_input[:,3]

'''
Plot input data - Figure 2d
'''
plt_fc.plot_input_data(data,xi,yi,zi,shape)
'''
These are the two parameters of our methodology for Euler deconvolution:
window size and the percentage of solutions to keep
'''
#moving data window size
winsize=7
#percentage of the solutions that will be keep
filt=0.1
#empty array for multiple SIs
est_classic=[]
#Define below the SIs to be tested
SI_vet=[0.001,1,2,3]
'''
Euler deconvolution for multiple SIs
'''
for SI in (SI_vet):
    classic = euler.euler_deconv(data,xi,yi,zi,shape,area,SI,winsize,
                                       filt)
    est_classic.append(classic)     
#Here finishes Euler deconvolution 
'''
Plot Figures 4 and 7 - Selected depth and base level estimates for all SIs
'''
plt_fc.plot_classic(data,est_classic,xi,yi,zi,shape)
''' 
Areas used to get the statistics - Defined after the classic plot
south,north,west,east
''' 
area_cla0=[0.,25000,24000,28000]
area_cla1=[9200,25000,15000,20000]
area_cla2=[14000,18000,5000,10000]
area_cla3=[5000,8000,5000,8000]

est_stats.classic(est_classic,area_cla0,SI_vet,'classic_plt0')
est_stats.classic(est_classic,area_cla1,SI_vet,'classic_plt1')
est_stats.classic(est_classic,area_cla2,SI_vet,'classic_plt2')
est_stats.classic(est_classic,area_cla3,SI_vet,'classic_plt3')


"""
To save the estimates in a txt file:

#convert the list to an array  
output=np.asarray(est_classic)  
#save the estimates in distinct files according to the SI  
for i in range(len(SI_vet)):  
	np.savetxt("estimates_SI_" + str(i) + ".txt",output[i],delimiter=",")
"""
