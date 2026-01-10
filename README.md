# Agent-based model of cranial neural crest migration

## Overview
This repository contains Python code implementing the agent-based model
described in Chapter 2 of the PhD thesis:

Mathematical models of single and collective cell migration in the cranial neural crest

The code initialises and runs simulations of the agent-based model and returns 
a visualisation of the results.

### Code Author
- Samuel Johnson

### Date
- 10/01/2025

### Requirements
- contourpy==1.3.1
- cycler==0.12.1
- fonttools==4.57.0
- imageio==2.37.0
- imageio-ffmpeg==0.6.0
- kiwisolver==1.4.8
- llvmlite==0.44.0
- matplotlib==3.10.1
- numba==0.61.2
- numpy==2.2.4
- packaging==24.2
- pillow==11.2.1
- pyparsing==3.2.3
- python-dateutil==2.9.0.post0
- scipy==1.15.2
- six==1.17.0

The required libraries can be installed from the requirements file using pip:

```bash
pip install requirements.txt
```

### Script Descriptions

#### VEGF.py
`VEGF.py` contains a forward-Euler solver used to update the chemoattractant profiles within the growing 
simulation domain. 

#### collisionCell.py 
`collisionCell.py` includes functions for cellular collision detection and the detection of cells with filopodia. 

#### growthFunction.py 
`growthFunction.py` includes functions that fit _in vivo_ data of the domain length to a logistic curve, and returns
a time-resolved list of domain lengths for use in the main simulation. 

#### insertCell.py 
`insertCell.py` includes functions that creates leader and follower cell objects. 

#### moveCell.py 
`moveCell.py` includes functions for cell movement according to leader-follower dynamics.  
 runs the main simulation and outputs a video or a .txt containing simulation data. 

### Execution 
Code is executed using the runSimulation.py: 

```bash
python runSimulation.py 
```
