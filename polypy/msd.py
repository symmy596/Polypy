import sys as sys
import numpy as np
from polypy import utils as ut
from scipy.constants import codata

kb = codata.value('Boltzmann constant')
ev = codata.value('electron volt')
ev = -ev


def square_distance(distance, n):
    '''Calculate the MSD for a series of distances

    Parameters
    ----------
    distance : array like
        Distance between atomic coordinates
    n : integer
        1 = 2D array, 0 = 1D array

    Returns
    -------
    msd : array like
        squared displacement
    '''
    if n == 1:
        msd = (distance[:, 0] ** 2) + (
               distance[:, 1] ** 2) + (
               distance[:, 2] ** 2)
    elif n == 0:
        msd = (distance[0] ** 2) + (
               distance[1] ** 2) + (
               distance[2] ** 2)
    return msd


def run_msd(trajectories, lv, timesteps, natoms, start, timestep):
    '''MSD calculator for a trajectory or series of trajectories

    Parameters
    ----------
    trajectories : array like
        atomic coordinates
    lv : array like
        Lattive Vectors
    timesteps : int
        Total Number of Timesteps
    natoms : int
        Total Number of Atoms
    start : int
        Total number of trajectory loops
    timestep : int
        Timestep of the simulation

    Returns
    -------
    msd_data : dictionary
        Dictionary containing 3D msd, 1D msd in the x, y, z directions
        and the time.
    '''
    trajectories = np.asarray(trajectories)
    msd_data = {'msd': np.array([]),
                'xmsd': np.array([]),
                'ymsd': np.array([]),
                'zmsd': np.array([]),
                'time': np.array([])}
    r0 = trajectories[start-1]
    rOd = trajectories[start-1]

    for j in range((start), timesteps):

        vec = lv[j]
        r1 = trajectories[j]
        distance_new = r1 - r0
        r1.tolist()
        rOd.tolist()

        if distance_new.size > 3:
            n = 1
            for k in range(0, distance_new[:, 0].size):
                for i in range(0, 3):
                    cross, r_new = ut.pbc(r1[k, i], rOd[k, i], vec[i])
                    if cross is True:
                        r1[k, i] = r_new
                        distance_new[k, i] = r_new - r0[k, i]
        else:
            n = 0
            r1 = r1.flatten()
            rOd = rOd.flatten()
            r0 = r0.flatten()
            distance_new = distance_new.flatten()
            for i in range(0, 3):

                cross, r_new = ut.pbc(r1[i], rOd[i], vec[i])
                if cross is True:
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
        msd_data['msd'] = np.append(msd_data['msd'], (msd_new))
        msd_data['time'] = np.append(msd_data['time'], ((j - start) * timestep))

        if n == 1:
            msd_data['xmsd'] = np.append(msd_data['xmsd'],
                               (np.average((distance[:, 0] ** 2))))
            msd_data['ymsd'] = np.append(msd_data['ymsd'],
                               (np.average((distance[:, 1] ** 2))))
            msd_data['zmsd'] = np.append(msd_data['zmsd'],
                               (np.average((distance[:, 2] ** 2))))
        elif n == 0:
            msd_data['xmsd'] = np.append(msd_data['xmsd'],
                               (np.average((distance[0] ** 2))))
            msd_data['ymsd'] = np.append(msd_data['ymsd'],
                               (np.average((distance[1] ** 2))))
            msd_data['zmsd'] = np.append(msd_data['zmsd'],
                               (np.average((distance[2] ** 2))))

    return msd_data


def check_trajectory(trajectory, xc, lv, timesteps, timestep, ul, ll, runs):
    '''Given an upper and lower limit of a 1D slice,
       determine if any part of a trajectory crosses
       a given 1D bin.
    Parameters
    ----------
    trajectory : array like
        Trajectories
    xc : array like
        Coordinates for one dimension
    lv : array like
        Lattice vectors
    timesteps : int
        Total Number of Timesteps
    timestep : float
        Timestep of simulation
    ul : float
        Upper Bin limit
    ll : float
        Lower Bin Limit
    runs : float
        Number of trajectory sweeps
    Return
    ------
    dco : float
        Diffusion Coefficient for a given atom within a given bin
    '''
    ib = False
    count = 0
    trajectory_slice = np.array([])
    dco = np.array([])
    vecs = np.array([])
    for i in range(0, xc.size):
        if xc[i] > ll and xc[i] < ul:
            ib = True
            count = count + 1
            trajectory_slice = np.append(trajectory_slice, trajectory[i])
            vecs = np.append(vecs, lv[i])

        elif xc[i] < ll or xc[i] > ul:
            if count > 200 and ib is True:

                trajectory_slice = np.split(trajectory_slice,
                                            (trajectory_slice.size / 3))
                vecs = np.reshape(vecs, (count, 3))
                do = np.array([])

                for i in range(0, runs):
                    start = i + 5
                    msd_data = run_msd(trajectory_slice,
                                       vecs,
                                       count,
                                       1,
                                       start,
                                       timestep)
                    d = ut.linear_regression(msd_data['time'],
                                             msd_data['msd'])[0]
                    d = ut.three_d_diffusion_coefficient(d)
                    do = np.append(do, d)

                dco = np.append(dco, np.average(do))
                count = 0
                trajectory_slice = np.array([])
                vecs = np.array([])
            else:
                ib = False
                trajectory_slice = np.array([])
                vecs = np.array([])
                count = 0
    if count > 200 and ib is True:

        do = np.array([])
        trajectory_slice = np.split(trajectory_slice,
                                    (trajectory_slice.size / 3))
        vecs = np.reshape(vecs, (count, 3))

        for i in range(0, runs):
            start = i + 5
            msd_data = run_msd(trajectory_slice,
                               vecs,
                               count,
                               1,
                               start,
                               timestep)
            d = ut.linear_regression(msd_data['time'],
                                     msd_data['msd'])[0]
            d = ut.three_d_diffusion_coefficient(d)
            do = np.append(do, d)

        dco = np.append(dco, np.average(do))
        count = 0

    return dco


def msd(data, timestep):
    '''Function that runs all of the parts of the MSD calcualtion.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : float
        Simulation timestep.

    Returns
    -------
    msd_data : dictionary
        Dictionary containing 3D msd, 1D msd in the x, y, z directions
        and the time.
    '''
    if data['timesteps'] == 1:
        print("ERROR: - Only one timestep has been found")
    if data['timesteps'] < 100:
        print("WARNING: Small number of timesteps - Poor statistics likely")
    if len(np.unique(data['label'])) > 1:
        print("ERROR: MSD can only handle one atom type. Exiting")
        sys.exit(0)

    trajectories = np.split(data['trajectories'], data['timesteps'])
    msd_data = run_msd(trajectories, data['lv'],
                       data['timesteps'],
                       data['natoms'],
                       1,
                       timestep)
    return msd_data


def smooth_msd(data, timestep, runs=5):
    '''MSD Launcher for a Smoothed MSD calc.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : float
        simulation timestep
    runs : int (optional)
        How many sweeps across the trajectory

    Returns
    -------
    msd_data : dictionary
        Dictionary containing 3D msd, 1D msd in the x, y, z directions
        and the time.
    '''
    smsd_data = {'time': np.array([]),
                 'msd':  np.array([]),
                 'xmsd': np.array([]),
                 'ymsd': np.array([]),
                 'zmsd': np.array([])}

    trajectories = np.split(data['trajectories'], data['timesteps'])

    for i in range(1, runs):
        start = i * 2
        msd_data = run_msd(trajectories, data['lv'],
                           data['timesteps'],
                           data['natoms'],
                           start,
                           timestep)
        smsd_data['msd'] = np.append(smsd_data['msd'], msd_data['msd'])
        smsd_data['xmsd'] = np.append(smsd_data['xmsd'], msd_data['xmsd'])
        smsd_data['ymsd'] = np.append(smsd_data['ymsd'], msd_data['ymsd'])
        smsd_data['zmsd'] = np.append(smsd_data['zmsd'], msd_data['zmsd'])
        smsd_data['time'] = np.append(smsd_data['time'], msd_data['time'])

    msd_data = {'time': ut.smooth_msd_data(smsd_data['time'],
                                           smsd_data['msd'])[0],
                 'msd':  ut.smooth_msd_data(smsd_data['time'],
                                            smsd_data['msd'])[1],
                 'xmsd': ut.smooth_msd_data(smsd_data['time'],
                                            smsd_data['xmsd'])[1],
                 'ymsd': ut.smooth_msd_data(smsd_data['time'],
                                            smsd_data['ymsd'])[1],
                 'zmsd': ut.smooth_msd_data(smsd_data['time'],
                                            smsd_data['ymsd'])[1]}
    return msd_data


def plane_msd(data, timestep, ul, ll, runs=1,
              direction="x"):
    '''Calculate an MSD value within a area of a structure.

    Parameters
    ----------
    data : dictionary
        Dictionary containing atom labels, trajectories,
        lattice vectors, total number of timesteps and atoms.
    timestep : float
        Simulation timestep.
    ul : float
        Upper bin limit.
    ll : float
        Lower bin limit.
    runs : int (optional)
        Number of trajectory sweeps.
    direction : str (optional)
        Direction normal to slices.

    Returns
    -------
    diffusion : float
        Diffusion coefficient in region.
    '''
    if direction == "x":
        val = 0
    elif direction == "y":
        val = 1
    elif direction == "z":
        val = 2

    xc = np.reshape(data['trajectories'][:, val], ((data['timesteps']),
                    data['natoms']))
    trajectories = np.split(data['trajectories'], data['timesteps'])
    trajectories = np.asarray(trajectories)
    d = np.array([])

    for i in range(0, (data['natoms'])):

        dd = check_trajectory(trajectories[:, i], xc[:, i], data['lv'],
                              data['timesteps'], timestep, ul, ll, runs)
        d = np.append(d, dd)
    diffusion = np.average(d)
    return diffusion
