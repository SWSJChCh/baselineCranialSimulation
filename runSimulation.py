'''
runSimulation.py - Samuel Johnson - 07/11/2025
'''

import imageio.v2 as imageio
import matplotlib.cm as cm
import matplotlib.animation as animation
import matplotlib.patches as patches
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import math
import shutil
import time
import os
import cv2
import copy
import sys
import random
import scipy
from scipy.ndimage import zoom
import datetime
from scipy import signal
from numba import jit, prange
from VEGF import *
from insertCell import *
from moveCell import *
from collisionCell import *
from growthFunction import *

#Filopodial persistence variable
filPersist = 1

#Animation Boolean
animate = True
#Timestep for animation
animStep = 10

#Create output directory
if animate:
    date = datetime.datetime.now().strftime('%H-%M-%S')
    os.makedirs('LeaderFollower' + date)

#Scale factor for mesh spacing (for computational speed-up)
meshScale = 1
#Simulation runtime
finTime = 24
#First cell insertion time
firstInsert = 6
#Timestep (h)
dt = 1 / 60
#List of domain lengths for each timestep (μm)
lengthList = domainLengths(int(finTime + 1))
#Domain length (μm)
Len = int(lengthList[0])
#Length of PDE mesh (for solver on unit-length mesh)
meshLen = int(lengthList[-1] / meshScale)
#Domain width (μm)
Wit = 120
#Width of PDE mesh (for solver on unit-length mesh)
meshWit = int(Wit / meshScale)
#Boundary condition smoothing parameter
bcParam = 0
#Number of leader cells
leadNum = 5
#Proportion of domain for which zero-flux boundary conditions are expressed
spanLen = 0
#Timesteps per attempted cell insertion
insertStep = 1
#Repeats for data averaging
repeats = 1

#VEGF parameters
D = 0.1                        #Diffusion constant
subStep = 1                    #Solver steps per timestep
dx = 1                         #Spacestep (x / μm)
dy = 1                         #Spacestep (y / μm)
c0 = 1.0                       #Initial concentration of reactant
xi = 0.1                       #Sensing parameter

#Cell parameters
cellRad = 7.5                          #Cell radius (μm)
searchRad = 5 * cellRad                #Box size for internalisation (μm)
leadSpeed = 1                          #Speed of leaders (μm / minute)
folSpeed = 1.3 * leadSpeed             #Speed of followers (μm / minute)
chi = 10**-4                           #Logistic production parameter
epsilon = 2 * cellRad                  #Distance for phenotype switch

lenFilo = 3.5 * cellRad                #Length of filopodia (μm)
lenFiloMax = 6 * cellRad               #Maximum length of cell-cell communication (μm)
filoNum = 3                            #Number of filopodia extended per timestep

lmda = 1 * 10**2                       #Internalisation parameter
chi = 10**-4                           #Logistic production parameter
epsilon = 2 * cellRad                  #Distance for phenotype switch

#Repeat simulations
for _ in range(repeats):
    #Initialise time
    t = 0
    #Initialise VEGF Mesh (for solver)
    VEGFMesh = createVEGFArray(meshWit, meshLen, c0, bcParam)
    #Initialise VEGF Array (for cells)
    zoom_factors = (Wit / VEGFMesh.shape[0], Len / VEGFMesh.shape[1])
    VEGFArray = zoom(VEGFMesh, zoom_factors, order=1)
    #List to store cell objects
    cellList = []
    #Data lists for visualisation
    cellMast = []
    filMast = []
    VEGFMast = []
    #Images for movie writer
    ims = []
    #Plot objects
    fig, ax = plt.subplots(1)
    #Counting variable
    counter = 0
    #Boolean for leader insertion
    leaderInsert = False
    #Run main simulation loop
    while (t < finTime):
        #Increase counting variable
        counter += 1
        #Update time
        t += dt
        #Actual domain length (μm)
        Len = int(lengthList[counter])
        #Time derivative of domain length
        lenDot = (lengthList[counter] - lengthList[counter-1]) / dt
        #Initial cells are leaders
        if not leaderInsert:
            #Create initial leader cells (evenly distributed line at LHS)
            initConfiguration(cellList, leadNum, Wit, cellRad, lenFilo)
            leaderInsert = True
        if t > firstInsert:
            #List to track cell data at time t
            cellCopyList = []
            #Insert follower cells at constant time intervals
            if counter % insertStep == 0:
                #Insert follower cell
                cell = followerCell(cellRad, lenFilo)
                insertCell(cell, cellList, Wit, Len)
            #Update chemicals according to PDE
            for _ in range(subStep):
                #Update chemoattractant
                VEGFMesh = updateVEGF(VEGFMesh, D, chi, lmda, cellRad, \
                                      cellList, dt, subStep, searchRad, \
                                      Len, Wit, lenDot, meshScale)
            #Update VEGF Array
            zoom_factors = (Wit / VEGFMesh.shape[0], Len / VEGFMesh.shape[1])
            VEGFArray = zoom(VEGFMesh, zoom_factors, order=1)

            filList = moveCells(VEGFArray, cellList, \
                            filoNum, lenFilo, lenFiloMax, \
                            xi, c0, cellRad, dx, dy, leadSpeed, folSpeed, spanLen, \
                            lengthList[counter-1], lengthList[counter], epsilon, filPersist)

            #List of cell position and VEGF Array
            if counter % animStep == 0 and animate:
                for i in cellList:
                    cellCopyList.append(copy.deepcopy(i))
                cellMast.append(cellCopyList)
                VEGFMast.append(VEGFArray.copy())
                filMast.append(filList.copy())

    #Chain any cells detached due to stochastic effects
    cellList = chainAtEnd(cellList, dx, dy, lenFilo, cellRad)

#Increase space between subplots
fig.tight_layout(pad=2.5)

#Produce .MP4 file of simulation
if animate:
    #Images
    ims = []
    for i in range(len(cellMast)):
        #Clear current axes
        ax.clear()
        #Set axes limits
        ax.set_xlim(0, Len)

        #Add filopodia patches
        for j in filMast[i]:
            rec0 = patches.Rectangle((round(j[0]), round(j[1])), width=0.15, \
                                      height=j[4], angle=j[2], color=j[3], \
                                      alpha=1, zorder=1)
            rec1 = patches.Rectangle((round(j[0]), round(j[1])), width=0.15, \
                                      height=j[4], angle=j[2], color=j[3], \
                                      alpha=1, zorder=1)
            ax.add_patch(rec0)

        #Add cell patches
        for j in cellMast[i]:
            if j.cellType == 'L':
                ax.add_patch(patches.Circle((j.x - 0.5, j.y - 0.5), \
                                cellRad, linewidth = 0, edgecolor = 'g', \
                                facecolor = 'g'))
            elif j.cellType == 'F' and j.attachedTo == 0:
                ax.add_patch(patches.Circle((j.x - 0.5, j.y - 0.5), \
                                cellRad, linewidth = 0, edgecolor = 'r', \
                                facecolor = 'r'))
            elif j.cellType == 'F' and j.attachedTo != 0:
                ax.add_patch(patches.Circle((j.x - 0.5, j.y - 0.5), \
                                cellRad, linewidth = 0, edgecolor = 'r', \
                                facecolor = 'r'))

        #Show VEGF Profile
        im0 = ax.imshow(VEGFMast[i], interpolation = 'none', vmin = 0, \
                           vmax = np.amax(VEGFMast[i]))

        #Title
        ax.set_title('Leader-Follower Simulation [24h]')

        #Colorbar
        cb0 = fig.colorbar(im0, shrink=0.125, aspect=3, ax=ax)

        #Save visualisation to folder
        plt.savefig('LeaderFollower{}/image{}.png'.format(date, i), dpi=400)
        #Remove colorbar for visualisation
        cb0.remove()

    #Produce video from folder
    with imageio.get_writer('LeaderFollower{}.mp4'.\
                            format(date), mode='I', fps=10) as writer:
        for i in range(len(cellMast)):
            filename = 'LeaderFollower{}/image{}.png'.format(date, i)
            image = imageio.imread(filename)
            writer.append_data(image)
            os.remove(filename)

#Delete folder used to make MP4
else:
    os.rmdir('LeaderFollower' + date)
