'''
insertCell.py - Samuel Johnson - 07/11/2025
'''

import math
import numpy as np

'''
Create leader cell
'''
class leaderCell:

    def __init__(self, radius, filLength, cellList):
        self.cellType = "L"                           #Cell type (Leader)
        self.radius = radius                          #Cell radius
        self.filLength = filLength                    #Filopodium length
        self.attachedTo = 0                           #Attachment ID (0 = None)
        #Create new chain for each leader cell
        self.chain = 1
        self.filoAge   = 0                            #Filopodium age variable

'''
Create follower cell
'''
class followerCell:

    def __init__(self, radius, filLength):
        self.cellType = "F"                           #Cell type (Follower)
        self.radius = radius                          #Cell radius
        self.filLength = filLength                    #Filopodium length
        self.attachedTo = 0                           #Attachment ID (0 = None)
        self.chain = 0                                #Cell chain (0 = None)

'''
Initial lattice configuration
'''
def initConfiguration(cellList, leadNum, width, radius, filLength):
    #List of leader objects for initial conditions
    initList = []
    #Create leader cells
    for _ in range(leadNum):
        initList.append(leaderCell(radius, filLength, cellList))
    #Evenly distributed y coordinates of leader cells
    yList = list(np.linspace(2 * radius, width - 2 * \
                 radius, len(initList)))
    #Update coordinates of leader cells and append cells to main cell list
    for i in range(len(initList)):
        initList[i].x = round(initList[i].radius)
        initList[i].y = yList[i]
        #Append created cells to main cell list
        cellList.append(initList[i])

'''
Insert cell into array
'''
def insertCell(cell, cellList, width, length):
    #Cells are inserted at LHS of domain
    xIns = round(cell.radius)
    #Cells are inserted into central migratory corridor
    yIns = np.random.uniform(round(cell.radius), \
           round(width - cell.radius))
    #Insertion Boolean
    insert = True
    #Determine if insertion causes overlap with other cell
    for i in cellList:
        if math.sqrt((i.x - xIns)**2 + (i.y - yIns)**2) < 2 * cell.radius:
            insert = False
    #Append cell to list of all cells
    if insert == True:
        cell.x = xIns
        cell.y = yIns
        cellList.append(cell)
    #Delete cell if overlap detected
    else:
        del cell
