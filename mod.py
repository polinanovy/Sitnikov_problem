#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import subprocess
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

def ReadParams(input_file):
    """
    Function for reading parameters from the input json-file.

    Parameters
    ----------
    input_file : string
        Name of the input json-file.

    Returns
    -------
    e : float
        Eccentricity.
    v0 : float
        Initial velocity for the falling body.
    atol : float
        Absolute tolerance. The solver keeps the local error estimates less than atol + rtol * abs(y)
    rtol : float
        Relative tolerance.
    method : string
        Method for integration.
    z1 : float
        The height to which the body falls.
    z2 : float
        Initial height for the falling body.
    step : float
        Step for heights from z1 to z2.
    n_rot : integer
        Number of rotations.
    save_dir : string
        Directory for saving the results.

    """
    with open(input_file, 'r') as f:
        s = json.load(f)
        
    e = s["input parameters"]["e"]
    v0 = s["input parameters"]["v0"]
        
    atol = s["integrator"]["atol"]
    rtol = s["integrator"]["rtol"]
    method = s["integrator"]["method"]
    
    z1 = s["phase map"]["z1"]
    z2 = s["phase map"]["z2"]
    step = s["phase map"]["step"]
    n_rot = s["phase map"]["n_rot"]
    save_dir = s["phase map"]["save_dir"]
    
    return e, v0, atol, rtol, method, z1, z2, step, n_rot, save_dir

def CreateDir(save_dir):
    """
    Function for creating directory for saving the results.

    Parameters
    ----------
    save_dir : string
        Directory for saving the results.

    Returns
    -------
    None.

    """
    if not(os.path.exists(save_dir)):
        key = True
        while key:
            ans = input("No path %s. Do you want to create path? (y/n) \n" %save_dir)
            if (ans == "y"):
                subprocess.run("mkdir -p %s" %save_dir, shell=True)  
                key = not(key)
            elif(ans == "n"):
                print("Exit")
                key = not(key)
                return
    return

def ODESystem(E, y, e):
    """
    Function for solving the main ODE system for height z and velocity v.

    Parameters
    ----------
    E : float
        Eccentric anomaly.
    y : list
        List of z and v.
    e : float
        Eccentricity.

    Returns
    -------
    array of float
        Array of values of differentials dz/dE and dv/dE.

    """
    z = y[0]
    v = y[1]
    r = 0.5 * (1 - e * np.cos(E))
    dz = 2 * r * v #dz/dE
    dv = -2 * r * z / ((r**2 + z**2)**(3/2)) #dv/dE 
    return np.array([dz, dv])

def PhaseMap(nrot, h, v, e, atol, rtol, method):
    """
    Function for solving the Poincare map. With function solve_ivp from scipy.integrate we can solve the main ODE system from fynction ODESystem.

    Parameters
    ----------
    nrot : integer
        Number of rotations.
    h : float
        Current height of the falling body.
    v : float
        Initial velocity of the falling body.
    e : float
        Eccentricity
    atol : float
        Absolute tolerance. The solver keeps the local error estimates less than atol + rtol * abs(y)
    rtol : float
        Relative tolerance.
    method : string
        Method for integration.

    Returns
    -------
    z : array of float
        Found from integration height.
    v : array of float
        Found from integration velocity.

    """
    t_eval = np.arange(0, nrot*2*np.pi, 2*np.pi)
    y0 = [h, v]
    sol = solve_ivp(ODESystem,[0, np.max(t_eval)], y0, method=method, rtol=rtol, atol=atol, t_eval=t_eval, args=(e,))
    E = sol.t
    z = sol.y[0]
    dz = sol.y[1] 
    r = 0.5 * (1 - e * np.cos(E))
    v = dz / r / 2
    return z, v

def SavePhaseMap(e, z1, z2, step, zv, save_dir):
    """
    Function for saving resultd of integration to '.dat'-files.

    Parameters
    ----------
    e : float
        Eccentricity.
    z1 : float
        The height to which the body falls.
    z2 : float
        Initial height for the falling body.
    step : float
        Step for heights from z1 to z2.
    zv : list
        List of found height z and velocity v from integration.
    save_dir : string
        Directory for saving the results.

    Returns
    -------
    files : list
        List of names of saved files. Important for the function for plotting Poincare map.

    """
    fname = 'e_' + str(e) + '_' + str(z1) + '_' + str(z2) + '_' + str(step)
    files = []
    for i,o in enumerate(zv):
        z = o[0]
        v = o[1]
        file = fname + '_' + str(i) + '_.dat'
        fop = open(Path(save_dir, file),'w')
        fop.write('# z   v' + '\n')
        for zi, vi in zip(z, v):
            fop.write(str(zi) + ' ' + str(vi) + '\n')
        files.append(file)
        fop.close()
    return files

def Plot(e, z2, files_dir, files):
    """
    Function for plotting Poincare map.

    Parameters
    ----------
    e : float
        Eccentricity.
    z2 : float
        Initial height for the falling body.
    files_dir : string
        Directory for saving the results.
    files : list
        List of names of saved files.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots()
    fig.set_size_inches(25,15)
    font = {'size'   : 20, 'family' : 'sans-serif'}
    mpl.rc('font', **font)
    plt.rc('text', usetex=True)
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$v$')    
    plt.title('e = ' + str(e))
    
    for i in files:
        fop = open(Path(files_dir, i),'r')
        z = []; v = []
        for line in fop.readlines():
            if line.startswith("#"):
                continue
            sl = line.split()
            z.append(float(sl[0]))
            v.append(float(sl[1]))
        ax.plot(z, v, 'o', markersize=0.8)
        fop.close()
        
    plt.xlim(-z2,z2)
    plt.ylim(-2,2)
    plt.savefig('e_' + str(e) + '.png')
    plt.close()
    return