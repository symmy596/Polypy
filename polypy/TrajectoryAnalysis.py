import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt

import Generic as ge
import Read as rd
import TrajectoryAnalysis as ta
import Write as wr


from scipy import stats
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev


def system_volume(lv, timesteps, timestep=None, output=None):  
    '''
    system volume - Calculate the volume at each timestep and return a volume as a function of time plot
    
    Parameters
    ----------
    first  : numpy
             Lattice vectors
    second : int
             Total number of timesteps
    third  : float
             timestep between records
    forth  : str
             output file name
            
    
    Returns
    -------
    first  : numpy
             volume at each timestep
    second : numpy
             time
    
    '''
    
    if timestep is None:
        print("No timestep provided - Exiting ")
        sys.exit(0)
        
    if output is None:
        filename = "Volume.png"
    
    volume = np.array([])
    time = np.array([])
    for i in range(0, timesteps):
        volume = np.append(volume, (np.prod(lv[i])))
        time = np.append(time, (i * timestep))
        
    wr.line_plot(time, volume, "Timestep", "System Volume (" r'$\AA$' ")", filename)

    return volume, time


def one_dimensional_density_sb(trajectories, ul=None, ll=None, direction=None):
    '''
    one_dimensional_density_sb - will return total number of species that spend a timestep within a bin range
    
    Parameters
    ----------
    first  : numpy 
             Atomic coordinates
    second : float
             Upper bin limit
    third  : float
             lower bin limit
    forth  : str
             Direction normal to slice
    
    Returns
    -------
    first  : Int
             Number of atomic species within bin
    '''
    
    if ul is None:
        print("No Upper Bin Limit provided - Please check and rerun")
    
    if ll is None:
        print("No Lower Bin Limit provided - Please check and rerun")
    
    if direction is None:
        direction = "x"
        print("No direction specified - Using default - x")
        
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2
        
    c = trajectories[:,val]
    l = (c.size)
    x = c.tolist()

    plane = 0
    
    for j in range(0, l):
              
        if c[j] > ll and c[j] < ul:
            plane = plane + 1       

    u = str(ul)
    l = str(ll)
    filename = "1D-Density-" + l + " - " + u             
         
    wr.one_dimensional_density_sb_output(plane, ul, ll, filename)   
    
    return plane
  
    
def one_dimensional_density(trajectories, timesteps, lv, Bin=None, direction=None, output=None):
    
    '''
    1D Atomic Density Analysis
    Last updated : 19/06/2018
    
    Parameters
    ----------
    
    first  : 2D numpy array
             Atomic Coordinates - (Number of Atoms * Number of Timesteps) x 3
    second : integer
             Number of timesteps
    third : 1D numpy array
             Lattice Vectors
    forth  : float
             Bin Value
    filth  : string
             Direction normal to the bin
    sixth  : str
             Output file name
        
    
    Returns
    -------
    
    text file
    matplolib plot
    "Conc.txt" - single column of concentrations in each bin
    "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
    
    '''
    if Bin is None:
        Bin = 0.1
        print("Bin Value not specified - Using Default - 0.1")
        
    if direction is None:
        direction = "x"
        print("No direction specified - Using default - x")
        
    if output is None:
        filename = "1D-Density.png"
      
    
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2

        
    c = trajectories[:,val]
    c = c + (np.average(lv[:,val]) / 2 ) 

    l = (c.size)
    x = c.tolist()
    
    x = np.amax(c) / Bin
    x = round(x, 0)
    x = int(x)
    
   
    bin_array = np.zeros((x))
    
    for j in range(0, l):
        plane = 0
        
        plane = ge.bin_choose(c[j], Bin)
        bin_array[plane] = bin_array[plane] + 1       

    x = np.arange( 0, ( bin_array.size ) )
    x = ((x * Bin)  - ( (np.average(lv[:,val]) / 2 )))
    y = ( bin_array / timesteps)
        
    wr.line_plot(x, y, "XCoordinate (" r'$\AA$' ")", "Number Density", filename)

    
def two_dimensional_density(trajectories, timesteps, lv, box=None, direction=None, output=None, log=None):
    
    '''
    2D Atomic Density Analysis
    Last updated : 19/06/2018
    
    
    Parameters
    ----------
    
    first     : 2D numpy array
                Atomic Coordinates - (Number of Atoms * Number of Timesteps) x 3
    second    : int
                number of timesteps
    third     : 1D numpy array
                Lattice Vectors
    forth     : float
                Box Value
    filth     : string
                Direction normal to the box
    sixth     : str
                output file name
    seventh   : boolean
                True for log plot, False for no log
    
    Returns
    -------
    
    matplolib plot
    "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
    
    '''
    
    if box is None:
        box = 0.1
        print("No box size specified - Using Default - 0.1")
    
    if direction is None:
        direction = "x"
        print("No direction specified - Using default - x")
        
    if output is None:
        filename = "2D-Density.png"
  
    if log == True:
        log = True
    else:
        log = False
        
    if direction == "x":
        val = [1, 2]
    elif direction == "y":
        val = [0, 2]
    elif direction == "z":
        val = [0, 1]

    xc = trajectories[:,val[0]]
    xc = xc + ( np.average(lv[:,[val[0]]]) / 2 )             
    yc = trajectories[:,val[1]]
    yc = yc + ( np.average(lv[:,[val[1]]]) / 2 ) 

    l = (xc.size)
    
    xc = xc.tolist()
    yc = yc.tolist()

    x = np.amax(xc) / box
    y = np.amax(yc) / box

    x = round(x, 0)
    x = int(x)
    
    y = round(y, 0)
    y = int(y)
    
    
    bin_array = np.zeros(((y), (x)))
    
    for j in range(0, l):        
        
        xbox = 0
        ybox = 0
        
        xbox = ge.bin_choose(xc[j], box)
        ybox = ge.bin_choose(yc[j], box)

        bin_array[ybox, xbox] = bin_array[ybox, xbox] + 1       

    bin_array = bin_array / timesteps

    x = np.arange((x))
    y = np.arange((y))
    
    x = ((x * box)) - (np.average(lv[:,[val[0]]]) / 2 )
    y = ((y * box))

    bin_array = bin_array + 0.001

    wr.contour_plot(x, y, bin_array, filename, log)


def conductivity(plane, volume, diff, temperature):
    '''
    conductivity - Calculate the ionic conductivity 
    
    Parameters
    ----------
    first  : int
             Total number of charge carriers
    second : numpy
             lattice vectors
    third  : float
             diffusion coefficient
    forth  : int
             Temperature
              
    Returns
    -------
    first   : float
              Conductivity
    '''

    volume = volume * (10 ** -30)
    diff = diff * (10 ** -9)
    conc = plane / volume
    
    EV = ev ** 2
    constants = kb * temperature
    conductivity = ((diff * conc) * EV) / constants
    
    return conductivity
                                                      
    
def average_position(trajectory, timesteps, natoms, lv):
    
    '''
    Average position calculator
    Parameters
    ----------
    first  : numpy 2D array
             Atomic Coorinates in one dimension - Atoms x Timesteps 
    second : integer
             Total number of timesteps
    third  : integer
             Total number of atoms
    fourth : numpy 1D array
             Lattice vectors
             
    Returns
    -------
    
    numpy 1D array
        average position in the desired dimension for each atom.
    
    '''
    
    average = np.array([])
    for i in range(0, natoms):
        for j in range(1, timesteps):
            cross, x_new = ge.pbc(trajectory[j,i], trajectory[(j-1),i], lv)
            if cross == True:
                trajectory[j,i] = x_new
        average = np.append(average, (np.average(trajectory[:,i])))
        
    return average


def distances(r1, r0):
    ''' 
    
    Simple subtraction 
    
    Parameters
    ----------
    first  : some sort of numpy object
    second : some sort of numpy object
    
    Returns
    -------
    floats
        result of second - first
        
    '''
    
    distance = ( r1 - r0 )
    return distance


def square_distance(distance, n):
    '''
    Calculate the MSD for a series of distances 
    
    Parameters 
    ----------
    first : 2D Numpy object
            Distance between atomic coordinates
    second : integer
             1 = 2D array
             0 = 1D array
    
    Return
    ------
    first : Numpy object 
            squared displacement
    '''
    
    if n == 1:
        msd_new = (distance[:,0] ** 2) + (distance[:,1] ** 2) + (distance[:,2] ** 2)
    elif n == 0:
        msd_new = (distance[0] ** 2) + (distance[1] ** 2) + (distance[2] ** 2)

    return msd_new


def run_msd(trajectories, lv, timesteps, natoms, start, timestep):
    '''
    
    MSD calculator - Common to all the various funcitons that do some sort of MSD
    
    Parameters
    ----------
    first  : 3D numpy array
             atomic coordinates 
    second : 1D numpy array
             Lattive Vectors
    third  : Integer
             Total Number of Timesteps
    forth  : Integer
             Total Number of Atoms
    filth  : Integer
             Total number of trajectory loops
    sixth  : Integer
             Timestep of the simulation
             
    Return
    ------
    first  : 1D numpy array
             MSD values
    second : 1D numpy array
             MSD values for the X direction
    third  : 1D numpy array
             MSD values for the y direction
    fourth : 1D numpy array
             MSD values for the z direction
    filth  : 1D numpy array
             Timesteps
    sixth  : 1D numpy array
             MSD arrays for every atom - needs reshaped
    
    '''
    
    trajectories = np.asarray(trajectories)
    msd = np.array([])
    xmsd = np.array([])
    ymsd = np.array([])
    time = np.array([])
    zmsd = np.array([])
    pmsd = np.array([])
    
    r0 = trajectories[start-1]
    rOd = trajectories[start-1] 
    for j in range((start), timesteps):
        vec = lv[j]
        r1 = trajectories[j]
        
        distance_new = distances(r1, r0)
        
        x = distance_new.size / 3
        x = int(x)
        r1.tolist()
        rOd.tolist()    
        if distance_new.size > 3:
            n = 1
            for k in range(0, x):
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
        pmsd = np.append(pmsd, msd_new)
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
        

    return msd, xmsd, ymsd, zmsd, time, pmsd


def msd_stats(msd, xmsd, ymsd, zmsd, time):
    
    '''
    Linear Regression 
    
    Parameters
    ----------
    first  : 1D numpy array
             MSD values
    second : 1D numpy array
             MSD values for the X direction
    third  : 1D numpy array
             MSD values for the y direction
    fourth : 1D numpy array
             MSD values for the z direction
    filth  : 1D numpy array
             Timesteps
             
    Return
    ------
    first  : float
             Overal gradient
    second : float
             gradient for X
    third  : float
             gradient for y
    fourth : float
             gradient for z
    '''
    
    ddiffusion, dintercept, dr_value, dp_value, dstd_err = stats.linregress(time, msd)
    xdiffusion, xintercept, xr_value, xp_value, xstd_err = stats.linregress(time, xmsd)
    ydiffusion, yintercept, yr_value, yp_value, ystd_err = stats.linregress(time, ymsd)
    zdiffusion, zintercept, zr_value, zp_value, zstd_err = stats.linregress(time, zmsd)
    
    return (ddiffusion, xdiffusion, ydiffusion, zdiffusion)


def diffusion_coefficient(ddiffusion, xdiffusion, ydiffusion, zdiffusion):
    
    '''
    Calculate the diffusion coefficient from the slope of MSD vs Time
    
    Parameters
    ----------
    first  : float
             Overal gradient
    second : float
             gradient for X
    third  : float
             gradient for y
    fourth : float
             gradient for z
    
    
    Return
    ------
    first  : float
             Overal Diffusion coefficient
    second : float
             Diffusion Coefficient for X
    third  : float
             Diffusion Coefficient for Y
    fourth : float
             Diffusion Coefficient for Z
    '''
    
    ddiffusion = ((np.average(ddiffusion)) / 6) * 10
    xdiffusion = ((np.average(xdiffusion)) / 2) * 10
    ydiffusion = ((np.average(ydiffusion)) / 2) * 10
    zdiffusion = ((np.average(zdiffusion)) / 2) * 10
    
    return ddiffusion, xdiffusion, ydiffusion, zdiffusion


def check_trajectory(trajectory, xc, lv, timesteps, timestep, ul, ll, runs):
    '''
    
    Check Trajectory - From an assigned bin determine if any part of a trajectory crosses the bin
    
    Parameters
    ----------
    first   : numpy array
              Trajectories
    second  : 2D numpy array
              Coordinates for one dimension
    third   : numpy array
              lattice vectors
    fourth  : integer
              Total Number of Timesteps
    filth   : float
              Timestep of simulation
    sixth   : float
              Upper Bin limit
    seventh : float 
              Lower Bin Limit
    eight   : integer
              Total number of trajectory loops
             
    Return
    ------
    float
        Diffusion Coefficient for a given atom within a given bin
    
    
    '''
        
    ib = False
    count = 0
    
    trajectory_slice = np.array([])
    diffusionco = np.array([])
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
                
                for i in range(0, runs):
                    start = i + 5

                    msd, xmsd, ymsd, zmsd, time, pmsd = ta.run_msd(trajectory_slice, vecs, count, 1, start, timestep)
                    d, xd, yd, zd = ta.msd_stats(msd, xmsd, ymsd, zmsd, time)
                    d, xd, yd, zd = ta.diffusion_coefficient(d, xd, yd, zd)
                    do = np.append(do, d)
                
                diffusionco = np.append(diffusionco, np.average(do))
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
        trajectory_slice = np.split(trajectory_slice, (trajectory_slice.size / 3))
        vecs = np.reshape(vecs, (count, 3))

        for i in range(0, runs):
            start = i + 5
            msd, xmsd, ymsd, zmsd, time, pmsd = ta.run_msd(trajectory_slice, vecs, count, 1, start, timestep)
            d, xd, yd, zd = ta.msd_stats(msd, xmsd, ymsd, zmsd, time)
            d, xd, yd, zd = ta.diffusion_coefficient(d, xd, yd, zd)
            do = np.append(do, d)
        diffusionco = np.append(diffusionco, np.average(do))
        count = 0

    return diffusionco


def msd(trajectories, lv, timesteps, natoms, timestep, conductivity=None, temperature=None):
    
    '''
    MSD Launcher
    
    Parameters
    ----------
    first  : 2D numpy array
             Atomic Coordinates 
    second : 1D numpy array
             Lattice vectors
    third  : integer
             total number of timesteps
    forth  : integer 
             total number of atoms
    
    Returns
    -------
    None
    '''
    if conductivity is None:
        conductivity = False
    if temperature is None and conductivity == True:
        print("Temperature is needed for conductivity calulcation- exiting....")
        sys.exit(0)
        
    trajectories = np.split(trajectories, timesteps)
    X = np.array([])
    start = 1
    
    msd, xmsd, ymsd, zmsd, time, pmsd = run_msd(trajectories, lv, timesteps, natoms, start, timestep)
    d, xd, yd, zd = msd_stats(msd, xmsd, ymsd, zmsd, time)    
    d, xd, yd, zd = diffusion_coefficient(d, xd, yd, zd)
    
    if conductivity == True:
       
        volume = (np.average(lv[:,0])) *  (np.average(lv[:,1])) * (np.average(lv[:,2]))
        cond = ta.conductivity(natoms, volume, d, temperature)
        wr.diffusion_output(d, xd, yd, zd, cond)

    else:
        wr.diffusion_output(d, xd, yd, zd)
    wr.msd_output(msd, xmsd, ymsd, zmsd, time)
    wr.msd_plot(time, msd, xmsd, ymsd, zmsd)

    
def smooth_msd(trajectories, lv, timesteps, natoms, timestep, runs=None, conductivity=None, temperature=None):
    
    '''
    MSD Launcher for a Smoothed MSD calc
    
    Parameters
    ----------
    first  : 2D numpy array
             Atomic Coordinates 
    second : 1D numpy array
             Lattice vectors
    third  : integer
             Number of MSD sweeps to carry out
    forth  : integer
             total number of timesteps
    filth  : integer 
             total number of atoms
    
    Returns
    -------
    None
    
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

    trajectories = np.split(trajectories, timesteps)
    x = np.array([])
    
    for i in range(1, runs):
        start = i * 10
        print("Starting Run", i, "of", runs)
        msd, xmsd, ymsd, zmsd, time, pmsd = run_msd(trajectories, lv, timesteps, natoms, start, timestep)
        d, xd, yd, zd = msd_stats(msd, xmsd, ymsd, zmsd, time)    
        
        dc = np.append(dc, d)
        xdc = np.append(xdc, xd)
        ydc = np.append(ydc, yd)
        zdc = np.append(zdc, zd)
        smsd = np.append(smsd, msd)
        sxmsd = np.append(sxmsd, xmsd)
        symsd = np.append(symsd, ymsd)
        szmsd = np.append(szmsd, zmsd)
        stime = np.append(stime, time)

    d, xd, yd, zd = diffusion_coefficient(dc, xdc, ydc, zdc)
    
    if conductivity == True:
        
        volume = (np.average(lv[:,0])) *  (np.average(lv[:,1])) * (np.average(lv[:,2]))
        cond = ta.conductivity(natoms, volume, d, temperature)
        wr.diffusion_output(d, xd, yd, zd, cond)
    else:
        wr.diffusion_output(d, xd, yd, zd)
    wr.msd_output(msd, sxmsd, symsd, szmsd, stime)
    wr.msd_plot(stime, smsd, sxmsd, symsd, szmsd)

        
def plane_msd(trajectories, lv, timesteps, natoms, timestep, runs=None, ul=None, ll=None, direction=None, conductivity=None, temperature=None):
    '''
    PlaneMSD - Calculate an MSD value within a area of a structure 
    
    Parameters
    ----------
    first  : 
    
    '''
    if runs is None:
        runs = 1
        
    if ul is None:
        print("No Upper bin limit has been provided - Exiting")
        sys.exit(0)
        
    if ll is None:
        print("No Lower bin limit has been provided - Exiting")
        sys.exit(0)
        
    if direction is None:
        direction = "x"
        
    if conductivity:
        conductivity = True
        c = trajectories
    else:
        con = False
        
    if temperature is None and conductivity == True: 
            print("The temperature of the simulation is needed to convert the calculated diffusion coeff to a conductivity")
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
                      
    xc = np.reshape(trajectories[:,val], ((timesteps), natoms))
    trajectories = np.split(trajectories, timesteps)
    trajectories = np.asarray(trajectories)
    d = np.array([])
    for i in range(0, (natoms)):
        dc = check_trajectory(trajectories[:,i], xc[:,i], lv, timesteps, timestep, ul, ll, runs)
        d = np.append(d, dc)
    nt = d.size
    diffusion = np.average(d)
    
    if conductivity == True:
        width = ul - ll
        volume = width * (np.average(lv[:,[area[0]]])) * (np.average(lv[:,[area[1]]]))
        plane = ta.one_dimensional_density_sb(c, ul=ul, ll=ll, direction=direction)
        plane = plane / timesteps
        cond = ta.conductivity(plane, volume, diffusion, temperature)
        wr.plane_msd_output(diffusion, ul, ll, nt, cond)

    else:
        wr.plane_msd_output(diffusion, ul, ll, nt)
    
    
def pmsd(trajectories, lv, timesteps, natoms, timestep, Bin=None, direction=None):
    '''
    MSD Launcher, will return the diffusion coefficient of every atom in the trajectory. 
    
    Parameters
    ----------
    first  : 2D numpy array
             Atomic Coordinates 
    second : 1D numpy array
             Lattice vectors
    third  : integer
             total number of timesteps
    fourth : integer 
             total number of atoms
    filth  : float
             Bin size for each plane
             
    Returns
    -------
    None
    
    '''
    if Bin is None:
        Bin = 5.0
        
    if direction is None:
        direction = "x"
          
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2
    
    xc = np.reshape(trajectories[:,0], ((timesteps), natoms))
    trajectories = np.split(trajectories, timesteps)

    x = np.array([])
    diffusion = np.array([])

    start = 1
    
    msd, xmsd, ymsd, zmsd, time, pmsd = run_msd(trajectories, lv, timesteps, natoms, start, timestep)

    pmsd = np.reshape(pmsd, ((timesteps - 1), natoms)) 
    for p in range(0, (natoms)):
        slope, intercept, r_value, p_value, std_err = stats.linregress(time, pmsd[:,p])
        diffusion = np.append(diffusion, slope)

    vec = np.average(lv[:,[val]])
    average = average_position(xc, timesteps, natoms, vec)
    ll = 0
    ul = 0
    coef = np.average(diffusion)
    bins = ge.bin_choose(vec, Bin)
    diff = np.array([])
    z = average.size
    slice_total = np.array([])

    for i in range(0, bins):
        ll = (i * Bin) - (vec * 0.5)
        ul = (ll + Bin)
        diff = np.append(diff, ((ul - (Bin / 2))))
       
        bin_diff = np.array([])
        for j in range(0, z):
            if average[j] > ll and average[j] < ul:
                bin_diff = np.append(bin_diff, diffusion[j])
        if bin_diff.size == 0:
            slice_total = np.append(slice_total, 0)
        else:
            slice_total = np.append(slice_total, (np.average(bin_diff))) 
    
    wr.pmsd_average_plot(diff, slice_total, coef, direction)
    wr.pmsd_plot(average, diffusion, direction)