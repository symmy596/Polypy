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


def one_dimensional_density_sb(Coords, NAtoms, NConfigs, UL=None, LL=None, Direction=None):
    '''
    one_dimensional_density_sb - will return total number of species that spend a timestep within a bin range
    
    Parameters
    ----------
    first  : numpy 
             Atomic coordinates
    second : Int
             Total number of atoms
    third  : Int
             Total number of timesteps
    fourth : float
             Upper bin limit
    filth  : float
             lower bin limit
    sixth  : str
             Direction normal to slice
    
    Returns
    -------
    first  : Int
             Number of atomic species within bin
    '''
    
    if UL:
        UL = UL
    else:
        print("No Upper Bin Limit provided - Please check and rerun")
    
    if LL:
        LL = LL
    else:
        print("No Lower Bin Limit provided - Please check and rerun")
    
    if Direction:
        Direction = Direction
    else:
        Direction = "x"
        print("No direction specified - Using default - x")
        
    if Direction == "x":
        Val = 0
    elif Direction == "y":
        Val = 1
    elif Direction == "z":
        Val = 2
    Total = []
                   
    C = Coords[:,Val]
    x = (C.size)
    X = C.tolist()

    Plane = 0
    
    for j in range(0, x):
              
        if C[j] > LL and C[j] < UL:
            Plane = Plane + 1       

    U = str(UL)
    L = str(LL)
    filename = "1D-Density-" + L + " - " + U             
         
    wr.one_dimensional_density_sb_output(Plane, UL, LL, filename)   
    
    return Plane
    
def one_dimensional_density(Coords, NAtoms, NConfigs, Vec, Bin=None, Direction=None, output=None):
    
    '''
    1D Atomic Density Analysis
    Last updated : 19/06/2018
    
    Parameters
    ----------
    
    first  : 2D numpy array
             Atomic Coordinates - (Number of Atoms * Number of Timesteps) x 3
    second : integer
             Number of Atoms
    third  : integer
             Number of timesteps
    fourth : 1D numpy array
             Lattice Vectors
    filth  : float
             Bin Value
    sixth  : string
             Direction normal to the bin
        
    
    Returns
    -------
    
    text file
    matplolib plot
    "Conc.txt" - single column of concentrations in each bin
    "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
    
    '''
    if output:
        filename = output + ".png"
    else:
        filename = "1D-Density.png"
    if Bin:
        Bin = Bin
    else:
        Bin = 0.1
        print("Bin Value not specified - Using Default - 0.1")
    if Direction:
        Direction = Direction
    else:
        Direction = "x"
        print("No direction specified - Using default - x")

        
        
    if Direction == "x":
        Val = 0
    elif Direction == "y":
        Val = 1
    elif Direction == "z":
        Val = 2
    Total = []
                   
    C = Coords[:,Val]
    C = C + (np.average(Vec[:,Val]) / 2 ) 

    x = (C.size)
    X = C.tolist()
    
    X = np.amax(C) / Bin
    X = round(X, 0)
    X = int(X)
    
   
    Bin_array = np.zeros((X))
    
    for j in range(0, x):
        Plane = 0
        
        Plane = ge.bin_choose(C[j], Bin)
        Bin_array[Plane] = Bin_array[Plane] + 1       

    X = np.arange( 0, ( Bin_array.size ) )
    X = ((X * Bin)  - ( (np.average(Vec[:,Val]) / 2 )))
    Y = ( Bin_array / NConfigs)
        
    wr.line_plot(X, Y, "XCoordinate (" r'$\AA$' ")", "Number Density", filename)
    
def two_dimensional_density(Coords, NAtoms, NConfigs, Vec, Box=None, Direction=None, output=None, log=None):
    
    '''
    2D Atomic Density Analysis
    Last updated : 19/06/2018
    
    
    Parameters
    ----------
    
    first   : 2D numpy array
              Atomic Coordinates - (Number of Atoms * Number of Timesteps) x 3
    second  : integer
              Number of Atoms
    third   : integer
              Number of timesteps
    fourth  : 1D numpy array
              Lattice Vectors
    filth   : float
              Box Value
    sixth   : string
              Direction normal to the box
    seventh : boolean
              True for log plot, False for no log
    
    Returns
    -------
    
    matplolib plot
    "Heatmap.png" - Contour plot - xy grid of total of atoms within each box
    
    To Do
    ------------------
    
    - Further Functionalise
    - Rewrite this into a class or something so that you just type "oxygen.two_dimensional_bin( BinSize, Direction)"
    
    '''
    if output:
        filename = output + ".png"
    else:
        filename = "2D-Density.png"
    if Box:
        Box = Box
    else:
        Box = 0.1
        print("No box size specified - Using Default - 0.1")
    if Direction:
        Direction = Direction
    else:
        Direction = "x"
        print("No direction specified - Using default - x")
    if log == True:
        log = True
    else:
        log = False
        
    if Direction == "x":
        Val = [1, 2]
    elif Direction == "y":
        Val = [0, 2]
    elif Direction == "z":
        Val = [0, 1]
    Total = []
    XCoords = Coords[:,Val[0]]
    XCoords = XCoords + ( np.average(Vec[:,[Val[0]]]) / 2 )             
    YCoords = Coords[:,Val[1]]
    YCoords = YCoords + ( np.average(Vec[:,[Val[1]]]) / 2 ) 

    x = (XCoords.size)
    XCoords = XCoords.tolist()

    y = (YCoords.size)
    YCoords = YCoords.tolist()

    X = np.amax(XCoords) / Box
    Y = np.amax(YCoords) / Box

    X = round(X, 0)
    X = int(X)
    
    Y = round(Y, 0)
    Y = int(Y)
    
    
    Bin_array = np.zeros(((Y), (X)))
    
    for j in range(0, x):        
        
        XBox = 0
        YBox = 0
        
        XBox = ge.bin_choose(XCoords[j], Box)
        YBox = ge.bin_choose(YCoords[j], Box)

        Bin_array[YBox, XBox] = Bin_array[YBox, XBox] + 1       

    Bin_array = Bin_array / NConfigs

    X = np.arange((X))
    Y = np.arange((Y))
    
    X = ((X * Box)  - ( np.amax(XCoords) / 2 ))
    Y = ((Y * Box))

    Bin_array = Bin_array + 0.001

    wr.contour_plot(X, Y, Bin_array, filename, log)

def system_volume(lv, NConfigs, timestep, output=None):
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
    
    
    if output:
        filename = output + ".png"
    else:
        filename = "Volume.png"
    
    volume = np.array([])
    time = np.array([])
    for i in range(0, NConfigs):
        Vec = lv[i]
        Vec = np.prod(Vec)
        volume = np.append(volume, Vec)
        time = np.append(time, (i * timestep))
        
    wr.line_plot(time, volume, "Timestep", "System Volume (" r'$\AA$' ")", filename)

    return volume, time

def conductivity(NConfigs, plane, volume, diff, temperature):
    '''
    conductivity - Calculate the ionic conductivity 
    
    Parameters
    ----------
    first   : int
              Total number of timesteps
    second  : int
              Total number of charge carriers
    third   : numpy
              lattice vectors
    seventh : float
              diffusion coefficient
    eight   : int
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
                                                      
    
def average_position(Coord, NConfigs, NAtoms, Vec):
    
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
    
    Average = np.array([])
    for i in range(0, NAtoms):
        for j in range(1, NConfigs):
            Cross, Xnew = ge.pbc(Coord[j,i], Coord[(j-1),i], Vec)
            if Cross == True:
                Coord[j,i] = Xnew
        Average = np.append(Average, (np.average(Coord[:,i])))
        
    return Average

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
    Distance = ( r1 - r0 )
    return Distance

def square_distance(Distance, N):
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
    if N == 1:
        MSD_new = (Distance[:,0] ** 2) + (Distance[:,1] ** 2) + (Distance[:,2] ** 2)
    elif N == 0:
        MSD_new = (Distance[0] ** 2) + (Distance[1] ** 2) + (Distance[2] ** 2)

    return MSD_new

def run_msd(Coords, Start, NConfigs, NAtoms, lv, timestep):
    '''
    
    MSD calculator - Common to all the various funcitons that do some sort of MSD
    
    Parameters
    ----------
    first : 3D numpy array
            atomic coordinates 
    second : 1D numpy array
             Lattive Vectors
    third  : Integer
             Starting position for MSD - Smoothed MSD changes the starting position to increase statistics
    fourth : Integer
             Total number of Timesteps
    filth  : Integer 
             Total number of atoms
             
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
    
    Coords = np.asarray(Coords)
    MSD = np.array([])
    XMSD = np.array([])
    YMSD = np.array([])
    Time = np.array([])
    ZMSD = np.array([])
    PMSD = np.array([])
    
    r0 = Coords[Start-1]
    rOd = Coords[Start-1] 
    for j in range((Start), NConfigs):
        Vec = lv[j]
        r1 = Coords[j]
        
        Distance_new = distances(r1, r0)
        
        x = Distance_new.size / 3
        x = int(x)
        r1.tolist()
        rOd.tolist()    
        if Distance_new.size > 3:
            N = 1
            for k in range(0, x):
                for i in range(0, 3):
                    Cross, r_new = ge.pbc(r1[k,i], rOd[k,i], Vec[i])
                    if Cross == True:
                        r1[k,i] = r_new
                        Distance_new[k,i] = r_new - r0[k,i]

        else:
            N = 0
            r1 = r1.flatten()
            rOd = rOd.flatten()
            r0 = r0.flatten()
            Distance_new = Distance_new.flatten()
            for i in range(0, 3):
                
                Cross, r_new = ge.pbc(r1[i], rOd[i], Vec[i])
                if Cross == True:
                    r1[i] = r_new
                    Distance_new[i] = r_new - r0[i]
        if N == 0:
            Distance = Distance_new.flatten()
        else:
            Distance = Distance_new

        r1 = np.asarray(r1)
        rOd = np.asarray(rOd)
        rOd = r1    

        MSD_new = square_distance(Distance, N)
        PMSD = np.append(PMSD, MSD_new)
        MSD_new = np.average(MSD_new)
        MSD = np.append(MSD, (MSD_new))
        Time = np.append(Time, ((j - Start) * timestep))

        if N == 1:
            XMSD = np.append(XMSD, (np.average((Distance[:,0] ** 2))))
            YMSD = np.append(YMSD, (np.average((Distance[:,1] ** 2))))
            ZMSD = np.append(ZMSD, (np.average((Distance[:,2] ** 2))))
        elif N == 0:
            XMSD = np.append(XMSD, (np.average((Distance[0] ** 2))))
            YMSD = np.append(YMSD, (np.average((Distance[1] ** 2))))
            ZMSD = np.append(ZMSD, (np.average((Distance[2] ** 2))))
        

    return MSD, XMSD, YMSD, ZMSD, Time, PMSD

def msd_stats(MSD, XMSD, YMSD, ZMSD, Time):
    
    '''
    Linear Regression 
    
    Parameters
    ----------
    first : 1D numpy array
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
    
    DDiffusion, DDintercept, DDr_value, DDp_value, DDstd_err = stats.linregress(Time, MSD)
    XDiffusion, Xintercept, Xr_value, Xp_value, Xstd_err = stats.linregress(Time, XMSD)
    YDiffusion, Yintercept, Yr_value, Yp_value, Ystd_err = stats.linregress(Time, YMSD)
    ZDiffusion, Zintercept, Zr_value, Zp_value, Zstd_err = stats.linregress(Time, ZMSD)
    
    return (DDiffusion, XDiffusion, YDiffusion, ZDiffusion)

def diffusion_coefficient(DDiffusion, XDiffusion, YDiffusion, ZDiffusion):
    
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
    
    DDiffusion = ((np.average(DDiffusion)) / 6) * 10
    XDiffusion = ((np.average(XDiffusion)) / 2) * 10
    YDiffusion = ((np.average(YDiffusion)) / 2) * 10
    ZDiffusion = ((np.average(ZDiffusion)) / 2) * 10
    return DDiffusion, XDiffusion, YDiffusion, ZDiffusion

def check_trajectory(NConfigs, XCoords, Coords, UL, LL, Runs, lv, timestep):
    '''
    
    Check Trajectory - From an assigned bin determine if any part of a trajectory crosses the bin
    
    Parameters
    ----------
    first  : integer
             Number of timesteps
    second : 2D numpy array
             Coordinates for one dimension - hard coded for x atm
    third  : 3D array 
             Coordinates
    fourth : float
             Upper limit of bin
    filth  : float
             Lower limit of bin
    sixth  : numpy array
             Lattice vectors
             
    Return
    ------
    float
        Diffusion Coefficient for a given atom within a given bin
    
    
    '''
    
    InBin = False
    Count = 0
    Trajectory = np.array([])
    DiffusionCo = np.array([])
    vecs = np.array([])
    for i in range(0, XCoords.size):

        if XCoords[i] > LL and XCoords[i] < UL:
            InBin = True
            Count = Count + 1
            Trajectory = np.append(Trajectory, Coords[i])
         
            vecs = np.append(vecs, lv[i])

        elif XCoords[i] < LL or XCoords[i] > UL:

            if Count > 200 and InBin == True:

                Trajectory = np.split(Trajectory, (Trajectory.size / 3))
                
                vecs = np.reshape(vecs, (Count, 3))
                DO = np.array([])
                for i in range(0, Runs):
                    Start = i + 5
                    MSD, XMSD, YMSD, ZMSD, Time, PMSD = ta.run_msd(Trajectory, Start, Count, 1, vecs, timestep)
                    D, XD, YD, ZD = ta.msd_stats(MSD, XMSD, YMSD, ZMSD, Time)
                    D, XD, YD, ZD = ta.diffusion_coefficient(D, XD, YD, ZD)
                    DO = np.append(DO, D)
                
                DiffusionCo = np.append(DiffusionCo, np.average(DO))
                Count = 0
                Trajectory = np.array([])
                vecs = np.array([])
            else:   

                InBin = False
                Trajectory = np.array([])
                vecs = np.array([])
                Count = 0
    if Count > 200 and InBin == True:

        DO = np.array([])
        Trajectory = np.split(Trajectory, (Trajectory.size / 3))
        vecs = np.reshape(vecs, (Count, 3))

        for i in range(0, Runs):
            Start = i + 5
            MSD, XMSD, YMSD, ZMSD, Time, PMSD = ta.run_msd(Trajectory, Start, Count, 1, vecs, timestep)
            D, XD, YD, ZD = ta.msd_stats(MSD, XMSD, YMSD, ZMSD, Time)
            D, XD, YD, ZD = ta.diffusion_coefficient(D, XD, YD, ZD)
            DO = np.append(DO, D)
        DiffusionCo = np.append(DiffusionCo, np.average(DO))
        Count = 0

    return DiffusionCo


def msd(Coords, NConfigs, NAtoms, timestep, lv, temperature=None, con=None):
    
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
    if con == True:
        con = True
    else:
        con = False
    if temperature:
        temperature = temperature
    else:
        temperature
        
    Coords = np.split(Coords, NConfigs)
    X = np.array([])
    Start = 1

    MSD, XMSD, YMSD, ZMSD, Time, PMSD = run_msd(Coords, Start, NConfigs, NAtoms, lv, timestep)
    DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo = msd_stats(MSD, XMSD, YMSD, ZMSD, Time)    
    DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo = diffusion_coefficient(DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo)
    
    if con == True:
       
        volume = (np.average(lv[:,0])) *  (np.average(lv[:,1])) * (np.average(lv[:,2]))
        conductivity = ta.conductivity(NConfigs, NAtoms, volume, DiffusionCo, temperature)
        wr.diffusion_output(DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo, conductivity)

    else:
        wr.diffusion_output(DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo)
        wr.msd_output(MSD, XMSD, YMSD, ZMSD, Time)
        wr.msd_plot(Time, MSD, XMSD, YMSD, ZMSD)

def smooth_msd(Coords, NConfigs, NAtoms, lv, timestep, Runs=None, conductivity=None, temperature=None):
    
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
    if conductivity:
        conductivity = True
    else:
        conductivity = False
    if temperature:
        temperature = temperature
    else:
        if conductivity == True:
            print("Temperature not provided - Temperature is needed for conductivity calculation")
            sys.exit(0)
            
    if Runs:
        Runs = Runs
    else:
        Runs = 5
        
    DiffusionCo = np.array([])
    XDiffusionCo = np.array([])
    YDiffusionCo = np.array([])
    ZDiffusionCo = np.array([])
    SMSD = np.array([])
    SXMSD = np.array([])
    SYMSD = np.array([])
    SZMSD = np.array([])
    STime = np.array([])

    Coords = np.split(Coords, NConfigs)
    X = np.array([])
    
    for i in range(1, Runs):
        Start = i * 10
        print("Starting Run", i, "of", Runs)
        MSD, XMSD, YMSD, ZMSD, Time, PMSD = run_msd(Coords, Start, NConfigs, NAtoms, lv, timestep)
        DDiffusion, XDiffusion, YDiffusion, ZDiffusion = msd_stats(MSD, XMSD, YMSD, ZMSD, Time)    
        
        DiffusionCo = np.append(DiffusionCo, DDiffusion)
        XDiffusionCo = np.append(XDiffusionCo, XDiffusion)
        YDiffusionCo = np.append(YDiffusionCo, YDiffusion)
        ZDiffusionCo = np.append(ZDiffusionCo, ZDiffusion)
        SMSD = np.append(SMSD, MSD)
        SXMSD = np.append(SXMSD, XMSD)
        SYMSD = np.append(SYMSD, YMSD)
        SZMSD = np.append(SZMSD, ZMSD)
        STime = np.append(STime, Time)

    DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo = diffusion_coefficient(DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo)
    
    if conductivity == True:
        volume = (np.average(lv[:,0])) *  (np.average(lv[:,1])) * (np.average(lv[:,2]))
        conductivity = ta.conductivity(NConfigs, NAtoms, volume, DiffusionCo, temperature)
        wr.diffusion_output(DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo, conductivity)
    else:
        wr.diffusion_output(DiffusionCo, XDiffusionCo, YDiffusionCo, ZDiffusionCo)
        wr.msd_output(SMSD, SXMSD, SYMSD, SZMSD, STime)
        wr.msd_plot(STime, SMSD, SXMSD, SYMSD, SZMSD)


def plane_msd(Coords, NConfigs, NAtoms, UL, LL, Direction, Runs, lv, timestep, con=None, temperature=None):
    '''
    PlaneMSD - Calculate an MSD value within a area of a structure 
    
    Parameters
    ----------
    first  : 
    
    '''
    if con:
        Con = True
        C = Coords
    else:
        Con = False
        
    if temperature: 
        temperature = temperature
    else:
        if con == True:
            print("The temperature of the simulation is needed to convert the calculated diffusion coeff to a conductivity")
        else:
            Temperature = 0
            
            
    if Direction == "x":
        Val = 0
        area = [1, 2]
    elif Direction == "y":
        Val = 1
        area = [0, 2]
    elif Direction == "z":
        Val = 2
        area = [0, 1]
                      
    XCoords = np.reshape(Coords[:,Val], ((NConfigs), NAtoms))
    Coords = np.split(Coords, NConfigs)
    Coords = np.asarray(Coords)
    Diffusion = np.array([])
    for i in range(0, (NAtoms)):
        DiffusionCo = check_trajectory(NConfigs, XCoords[:,i], Coords[:,i], UL, LL, Runs, lv, timestep)
        Diffusion = np.append(Diffusion, DiffusionCo)
    nt = Diffusion.size
    Diffusion = np.average(Diffusion)
    
    if con == True:
        width = UL - LL
        volume = width * (np.average(lv[:,[area[0]]])) * (np.average(lv[:,[area[1]]]))
        plane = ta.one_dimensional_density_sb(C, NAtoms, NConfigs, UL=UL, LL=LL, Direction="x")
        plane = plane / NConfigs
        conductivity = ta.conductivity(NConfigs, plane, volume, Diffusion, temperature)
        wr.plane_msd_output(Diffusion, UL, LL, nt, conductivity)

    else:
        wr.plane_msd_output(Diffusion, UL, LL)
    
    
def pmsd(Coords, lv, NConfigs, NAtoms, timestep, Bin=None, Direction=None):
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
    if Bin:
        Bin = Bin
    else:
        Bin = 5.0
        
    if Direction:
        Direction = Direction
    else:
        Direction = "x"
        
          
    if Direction == "x":
        Val = 0
    elif Direction == "y":
        Val = 1
    elif Direction == "z":
        Val = 2
    
    XCoords = np.reshape(Coords[:,0], ((NConfigs), NAtoms))
    Coords = np.split(Coords, NConfigs)

    X = np.array([])
    Diffusion = np.array([])

    Start = 1
    
    MSD, XMSD, YMSD, ZMSD, Time, PMSD = run_msd(Coords, Start, NConfigs, NAtoms, lv, timestep)

    PMSD = np.reshape(PMSD, ((NConfigs - 1), NAtoms)) 
    for p in range(0, (NAtoms)):
        DDiffusion, DDintercept, DDr_value, DDp_value, DDstd_err = stats.linregress(Time, PMSD[:,p])
        Diffusion = np.append(Diffusion, DDiffusion)

    Vec = np.average(lv[:,[Val]])
    Average = average_position(XCoords, NConfigs, NAtoms, Vec)
    LL = 0
    UL = 0
    Coef = np.average(Diffusion)
    Bins = ge.bin_choose(Vec, Bin)
    Diff = np.array([])
    Z = Average.size
    BinDiffTot = np.array([])

    for i in range(0, Bins):
        LL = (i * Bin) - (Vec * 0.5)
        UL = (LL + Bin)
        Diff = np.append(Diff, ((UL - (Bin / 2))))
       
        BinDiff = np.array([])
        for j in range(0, Z):
            if Average[j] > LL and Average[j] < UL:
                BinDiff = np.append(BinDiff, Diffusion[j])
        if BinDiff.size == 0:
            BinDiffTot = np.append(BinDiffTot, 0)
        else:
            BinDiffTot = np.append(BinDiffTot, (np.average(BinDiff))) 
    
    wr.pmsd_average_plot(Diff, BinDiffTot, Coef, Direction)
    wr.pmsd_plot(Average, Diffusion, Direction)