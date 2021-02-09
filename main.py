"""
Created on Sun Feb  7 23:10:41 2021
@author: necbal

Generates a dynamical system in 2d of a  given by diffusion (LaplaceR)
plus reaction terms, that when start outside equilibrium
start generating non-trivial patterns

Follows https://www.youtube.com/watch?v=BV9ny785UNc
"""
import matplotlib.pyplot as plt
from aux_functions import LaplaceR, cmap, norm, subtract

# Parameters for diffusion
d_B = .5
d_A = 1
# Parameters for reaction
f = .055
k = .062


steps = 5000 #time steps in the dynamics
steps_plot = 1 #steps between plots shown
gs = 100 #gridsize

#initial condition of chemicals A and B
Grid_B = [[ 0 for i in range(0,gs)]
         for i in range(0,int(gs))]
Grid_A = [[0 for i in range(0,gs)]
        for i in range(0,int(gs))]

#Add two squares of chemical B as a perturbation of the initial conditon
for i in range(10,40):
    for j in range(10,40):
        Grid_B[i][j] = 1
for i in range(60,80):
    for j in range(60,80):
        Grid_B[i][j] = 1 

System = {'B': Grid_B, 'A': Grid_A }
Next_B = Grid_B[:]
Next_A = Grid_A[:]

for t in range(0,steps):
        printtime = t/steps_plot
        if printtime.is_integer() == True: 
            fig, ax = plt.subplots(1,1)
            uminv = [subtract(System['A'][n], System['B'][n]) for n in range(0,gs)]
            ax.imshow(uminv, cmap=cmap, norm=norm)
            ax.set_title('A-B: t = ' + str(t))
            plt.show()
        for i in range(0,gs-1):
            for j in range(0,gs-1): #next is the dynamics
                Next_B[i][j] = System['B'][i][j]+ d_B*LaplaceR(System['B'],i,j)-(k+f)*System['B'][i][j]+pow(System['B'][i][j],2)*System['A'][i][j]
                Next_A[i][j] = System['A'][i][j] + d_A*LaplaceR(System['A'],i,j)+f*(1-System['A'][i][j])-pow(System['B'][i][j],2)*System['A'][i][j]
        System['B'], Next_B = Next_B, System['B'] #swap lists
        System['A'], Next_A = Next_A, System['A'] #swap lists

""" 
#Discrete:...
#PDE:...
# Steady state initial condition:  A = .85, B=-.3
                                     A = .5, B=-.7 
                                     
                                     
#Interesting Initial conditions:
Grid_B = [[ .7 for i in range(0,gs)]
         for i in range(0,int(gs))]
Grid_A = [[0.1 for i in range(0,gs)]
        for i in range(0,int(gs))]
        
    """

   