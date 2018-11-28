import os as os
import sys as sys
import numpy as np

from polypy import Read as rd
from polypy import Density as Dens
from polypy import Utils as Ut
from polypy import Write as wr
from polypy import Generic as ge

from scipy import stats
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev
                                              

def distances(r1, r0):
    ''' 
    Simple subtraction 
    Parameters
    ----------
    r1       : numpy object
    r0       : numpy object
    
    Returns
    -------
    distance : Float
    '''
    distance = ( r1 - r0 )
    return distance

def square_distance(distance, n):
    '''
    Calculate the MSD for a series of distances 
    Parameters 
    ----------
    distance : Distance between atomic coordinates     : 2D Numpy object
    n        : 1 = 2D array, 0 = 1D array              : Integer

    Return
    ------
    msd      : squared displacement                    : Numpy object 
    '''
    if n == 1:
        msd = (distance[:,0] ** 2) + (distance[:,1] ** 2) + (distance[:,2] ** 2)
    elif n == 0:
        msd = (distance[0] ** 2) + (distance[1] ** 2) + (distance[2] ** 2)

    return msd

def run_msd(trajectories, lv, timesteps, natoms, start, timestep):
    '''
    MSD calculator - Common to all the various funcitons that do some sort of MSD
    Parameters
    ----------
    trajectories  : atomic coordinates                  : 3D numpy array
    lv            : Lattive Vectors                     : 1D numpy array
    timesteps     : Total Number of Timesteps           : Integer
    natoms        : Total Number of Atoms               : Integer
    start         : Total number of trajectory loops    : Integer
    timestep      : Timestep of the simulation          : Integer
            
    Return
    ------
    msd_data : Dictionary {'msd': msd, 'xmsd': xmsd, 'ymsd': ymsd, 'zmsd': zmsd, 'time': time}    
    pmsd     : MSD arrays for every atom          :  1D numpy array
    '''
    trajectories = np.asarray(trajectories)
    msd = np.array([])
    xmsd = np.array([])
    ymsd = np.array([])
    time = np.array([])
    zmsd = np.array([])
    r0 = trajectories[start-1]
    rOd = trajectories[start-1] 

    for j in range((start), timesteps):
        
        vec = lv[j]
        r1 = trajectories[j]
        distance_new = distances(r1, r0)
        r1.tolist()
        rOd.tolist()    
        
        if distance_new.size > 3:
            n = 1

            for k in range(0, distance_new[:,0].size):
                for i in range(0, 3):
                    cross, r_new = ge.pbc(r1[k,i], rOd[k,i], vec[i])
                    if cross == True:
                        r1[k,i] = r_new
                        distance_new[k,i] = r_new - r0[k,i]

        else:
            n = 0
            r1 = r1.flatten()
            rOd = rOd.flatten()
            r0 = r0.flatten()
            distance_new = distance_new.flatten()
            for i in range(0, 3):
                
                cross, r_new = ge.pbc(r1[i], rOd[i], vec[i])
                if cross == True:
                    r1[i] = r_new
                    distance_new[i] = r_new - r0[i]
        if n == 0:
            distance = distance_new.flatten()
        else:
            distance = distance_new

        r1 = np.asarray(r1)
        rOd = np.asarray(rOd)
        rOd = r1    

        msd_new = square_distance(distance, n)
        msd_new = np.average(msd_new)
        msd = np.append(msd, (msd_new))
        time = np.append(time, ((j - start) * timestep))

        if n == 1:
            xmsd = np.append(xmsd, (np.average((distance[:,0] ** 2))))
            ymsd = np.append(ymsd, (np.average((distance[:,1] ** 2))))
            zmsd = np.append(zmsd, (np.average((distance[:,2] ** 2))))
        elif n == 0:
            xmsd = np.append(xmsd, (np.average((distance[0] ** 2))))
            ymsd = np.append(ymsd, (np.average((distance[1] ** 2))))
            zmsd = np.append(zmsd, (np.average((distance[2] ** 2))))
        
        msd_data = {'msd': msd, 'xmsd': xmsd, 'ymsd': ymsd, 'zmsd': zmsd, 'time': time}

    return msd_data

def check_trajectory(trajectory, xc, lv, timesteps, timestep, ul, ll, runs):
    '''
    Check Trajectory - From an assigned bin determine if any part of a trajectory crosses the bin
    Parameters
    ----------
    trajectory   : Trajectories                    : Numpy array
    xc           : Coordinates for one dimension   : Numpy array
    lv           : Lattice vectors                 : Numpy array
    timesteps    : Total Number of Timesteps       : Integer
    timestep     : Timestep of simulation          : Float
    ul           : Upper Bin limit                 : Float
    ll           : Lower Bin Limit                 : Float
    runs         : Number of trajectory sweeps     : Intger
             
    Return
    ------
    Diffusion Coefficient for a given atom within a given bin
    '''
    ib = False
    count = 0
    trajectory_slice = np.array([])
    dco = np.array([])
    xco = np.array([])
    yco = np.array([])
    zco = np.array([])
    vecs = np.array([])
    
    for i in range(0, xc.size):
        if xc[i] > ll and xc[i] < ul:
            ib = True
            count = count + 1
            trajectory_slice = np.append(trajectory_slice, trajectory[i])
            vecs = np.append(vecs, lv[i])

        elif xc[i] < ll or xc[i] > ul:
            if count > 200 and ib == True:

                trajectory_slice = np.split(trajectory_slice, (trajectory_slice.size / 3))
                vecs = np.reshape(vecs, (count, 3))
                do = np.array([])
                xo = np.array([])
                yo = np.array([])
                zo = np.array([])

                for i in range(0, runs):
                    start = i + 5
                    msd_data, pmsd = ta.run_msd(trajectory_slice, vecs, count, 1, start, timestep)
                    d, xd, yd, zd = ta.msd_stats(msd_data)
                    d, xd, yd, zd = ta.diffusion_coefficient(d, xd, yd, zd)
                    do = np.append(do, d)
                    xo = np.append(xo, xd)
                    yo = np.append(yo, yd)
                    zo = np.append(zo, zd)

                dco = np.append(dco, np.average(do))
                xco = np.append(xco, np.average(xo))
                yco = np.append(yco, np.average(yo))
                zco = np.append(zco, np.average(zo))

                count = 0
                trajectory_slice = np.array([])
                vecs = np.array([])
            else:   
                ib = False
                trajectory_slice = np.array([])
                vecs = np.array([])
                count = 0
    if count > 200 and ib == True:

        do = np.array([])
        xo = np.array([])
        yo = np.array([])
        zo = np.array([])
        trajectory_slice = np.split(trajectory_slice, (trajectory_slice.size / 3))
        vecs = np.reshape(vecs, (count, 3))

        for i in range(0, runs):
            start = i + 5
            msd_data, pmsd = ta.run_msd(trajectory_slice, vecs, count, 1, start, timestep)
            d, xd, yd, zd = ta.msd_stats(msd_data)
            d, xd, yd, zd = ta.diffusion_coefficient(d, xd, yd, zd)
            do = np.append(do, d)
            xo = np.append(xo, xd)
            yo = np.append(yo, yd)
            zo = np.append(zo, zd)
        
        dco = np.append(dco, np.average(do))
        xco = np.append(xco, np.average(xo))
        yco = np.append(yco, np.average(yo))
        zco = np.append(zco, np.average(zo))
        count = 0

    return dco, xco, yco, zco

def msd(data, timestep):
    '''
    msd - Function that runs all of the parts of the MSD calcualtion
    Parameters
    ----------
    data          : Dictionary containing atomic trajectories ['trajectories'], lattice vectors ['lv'], timesteps ['timesteps'] and number                     of atoms['natoms']
    timestep      : simulation timestep - float                                 : Float
    conductivity  : True/False True - calculate conductivity                    : Bool
    temperature   : Temperature of the simulation - needed for conductivity     : Integer
    
    Returns
    -------
    Outputs diffusion info
    '''
     
    if data['timesteps'] == 1:
        print("ERROR: - Only one timestep has been found")
    if data['timesteps'] < 100:
        print("WARNING: Small number of timesteps - Poor statistics likely")
    if len(np.unique(data['label'])) > 1:
        print("ERROR: MSD can only handle one atom type. Exiting...")
        sys.exit(0)
    

    trajectories = np.split(data['trajectories'], data['timesteps'])
    msd_data = run_msd(trajectories, data['lv'], data['timesteps'], data['natoms'], 1, timestep)

    return msd_data

def smooth_msd(data, timestep, runs=None, conductivity=None, temperature=None):
    '''
    smooth_msd - MSD Launcher for a Smoothed MSD calc
    Parameters
    ----------
    data          : Dictionary containing atomic trajectories ['trajectories'], lattice vectors ['lv'], timesteps ['timesteps'] and number                     of atoms['natoms']
    timestep      : simulation timestep - float                                 : Float
    runs          : How many sweeps across the trajectory                       : Integer
    conductivity  : True/False True - calculate conductivity                    : Bool
    temperature   : Temperature of the simulation - needed for conductivity     : Integer
    
    Returns
    -------
    Outputs diffusion info
    '''
    if conductivity is None:
        conductivity = False
    if temperature is None and con == True:
        print("Temperature is needed for conductivity calulcation- exiting....")
        sys.exit(0)
            
    if runs is None:
        runs = 5
        
    dc = np.array([])
    xdc = np.array([])
    ydc = np.array([])
    zdc = np.array([])
    smsd = np.array([])
    sxmsd = np.array([])
    symsd = np.array([])
    szmsd = np.array([])
    stime = np.array([])

    trajectories = np.split(data['trajectories'], data['timesteps'])
    x = np.array([])
    
    for i in range(1, runs):
        start = i * 10
        print("Starting Run", i, "of", runs)
        msd_data, pmsd = run_msd(trajectories, data['lv'], data['timesteps'], data['natoms'], start, timestep)
        d, xd, yd, zd = msd_stats(msd_data)    
        
        dc = np.append(dc, d)
        xdc = np.append(xdc, xd)
        ydc = np.append(ydc, yd)
        zdc = np.append(zdc, zd)
        smsd = np.append(smsd, msd_data['msd'])
        sxmsd = np.append(sxmsd, msd_data['xmsd'])
        symsd = np.append(symsd, msd_data['ymsd'])
        szmsd = np.append(szmsd, msd_data['zmsd'])
        stime = np.append(stime, msd_data['time'])

    d, xd, yd, zd = diffusion_coefficient(dc, xdc, ydc, zdc)
    smsd_data = {'time': stime, 'msd': smsd, 'xmsd': sxmsd, 'ymsd': symsd, 'zmsd': szmsd}
    if conductivity == True:
        
        volume = (np.average(data['lv'][:,0])) *  (np.average(data['lv'][:,1])) * (np.average(data['lv'][:,2]))
        cond = ta.conductivity(data['natoms'], volume, d, temperature)
        wr.diffusion_output(d, xd, yd, zd, cond)
    else:
        wr.diffusion_output(d, xd, yd, zd)
    wr.msd_plot(smsd_data)
       
def plane_msd(data, timestep, runs=None, ul=None, ll=None, direction=None, conductivity=None, temperature=None):
    '''
    PlaneMSD - Calculate an MSD value within a area of a structure 
    Parameters
    ----------
    data          : Dictionary containing atomic trajectories ['trajectories'], lattice vectors ['lv'], timesteps ['timesteps'] and number of atoms['natoms']    
    timestep      : Simulation timestep                        : Float
    runs          : Number of trajectory sweeps                : Integer
    ul            : Upper bin limit                            : Float
    ll            : Lower bin limit                            : Float
    direction     : Direction normal to slices                 : str
    conductivity  : True/False - True for conductivity calc    : Bool
    temperature   : Temperature of simulation                  : Integer
    
    Returns
    -------
    File containing the result of the calc
    '''
    if runs is None:
        runs = 1
    if ul is None:
        sys.exit(0)
    if ll is None:
        sys.exit(0)
    if direction is None:
        direction = "x"
    if conductivity:
        conductivity = True
        c = trajectories
    else:
        con = False
    if temperature is None and conductivity == True: 
        sys.exit(0)
            
    if direction == "x":
        val = 0
        area = [1, 2]
    elif direction == "y":
        val = 1
        area = [0, 2]
    elif direction == "z":
        val = 2
        area = [0, 1]
                      
    xc = np.reshape(data['trajectories'][:,val], ((data['timesteps']), data['natoms']))
    trajectories = np.split(data['trajectories'], data['timesteps'])
    trajectories = np.asarray(trajectories)
    d = np.array([])
    x = np.array([])
    y = np.array([])
    z = np.array([])

    for i in range(0, (data['natoms'])):

        dd, xd, yd, zd = check_trajectory(trajectories[:,i], xc[:,i], data['lv'], data['timesteps'], timestep, ul, ll, runs)
        d = np.append(d, dd)
        x = np.append(x, xd) 
        y = np.append(y, yd) 
        z = np.append(z, zd) 

    nt = d.size
    diffusion = np.average(d)
    xd = np.average(x)
    yd = np.average(y)
    zd = np.average(z)

    if conductivity == True:

        width = ul - ll
        volume = width * (np.average(data['lv'][:,[area[0]]])) * (np.average(data['lv'][:,[area[1]]]))
        plane = ta.one_dimensional_density_sb(c, ul=ul, ll=ll, direction=direction)
        plane = plane / data['timesteps']
        cond = ta.conductivity(plane, volume, diffusion, temperature)
        wr.plane_msd_output(diffusion, xd, yd, zd, ul, ll, nt, cond)

    else:
        wr.plane_msd_output(diffusion, xd, yd, zd, ul, ll, nt)