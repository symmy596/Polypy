import os as os
import sys as sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math as mt
import seaborn as sns
sns.set(style="white")
sns.set_style("ticks")


def msd_plot(msd_data, set_style="default", palette="tab10", figsize=None, output=None):
    '''
    MSDPlot - Plot MSD 
    Parameters 
    ----------
    msd_data  : Dictionary {'msd': msd, 'xmsd': xmsd, 'ymsd': ymsd, 'zmsd': zmsd, 'time': time}

    Return
    ------
    matplotlib plot - png
    '''
    sns.palette=palette
    plt.style.use(set_style)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    ax.set_ylim(ymin=0, ymax=np.amax(msd_data['msd']))
    ax.set_xlim(xmin=0, xmax=np.amax(msd_data['time']))
    ax.plot(msd_data['time'], msd_data['msd'], label="MSD")
    ax.plot(msd_data['time'], msd_data['xmsd'], label="XMSD")
    ax.plot(msd_data['time'], msd_data['ymsd'], label="YMSD")
    ax.plot(msd_data['time'], msd_data['zmsd'], label="ZMSD")
    ax.tick_params(labelsize=12)
    ax.set_xlabel("Time (ps)", fontsize=15)
    ax.set_ylabel("MSD ($\AA$)", fontsize=15)
    plt.legend()
    if output:
        plt.savefig(output, dpi=600)
    plt.show()
    plt.close()

def volume_plot(x, y, xlab="Timestep (ps)", ylab="System Volume ($\AA$)",
                output=None, set_style="default", palette="tab10",
                figsize=None):
    '''Plots the system volume vs timestep.
    
    Parameters
    ----------
    x : array like
        Timesteps
    y : array like
        Volume
    xlab : str
        X label
    ylab : str
        Y label
    output : str
        Output filename
    set_style : str
        Plot style
    palette : str
        Color palette
    figsize : tuple (optional)
        Size of plot

    Return
    ------
    matplotlib plot
    
    '''
    sns.palette=palette
    plt.style.use(set_style)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    ax.plot(x, y)
    ax.set_xlabel(xlab, fontsize=13)
    ax.set_ylabel(ylab, fontsize=13)
    ax.tick_params(labelsize=12)
    if output:
        plt.savefig(output, dpi=600)
    plt.tight_layout()
    plt.show()
    plt.close()

def line_plot(X, Y, XLab, YLab, output):
    '''
    LinePlot - Simple line plot
    
    Parameters
    ----------
    first  : numpy object
             X axis values
    second : numpy object
             Y axis values
    third  : str
             X label
    fourth : str
             Y label
             
    Return
    ------
    matplotlib plot
    
    '''
        
    plt.plot(X, Y, color="crimson")
    plt.xlabel(XLab, fontsize=13)
    plt.ylabel(YLab, fontsize=13)
    plt.tight_layout()
    plt.tick_params(labelsize=12)
    plt.savefig(output, dpi=600)
    plt.show()
    plt.close()

def contour_plot(X, Y, Z, output):
    '''
    CountourPlot - Contour plotting tool
    
    Parameters 
    ----------
    first  : numpy 
             X axis
    second : numpy 
             Y axis
    third  : 2D numpy array
             grid
             
    Return
    ------
    matplotlib plot
    
    '''
    plt.contourf(X, Y, Z, cmap="gray")
    plt.xlabel("X Coordinate (" r'$\AA$' ")", fontsize=15)
    plt.ylabel("Y Coordinate (" r'$\AA$' ")", fontsize=15)
    plt.tick_params(labelsize=12)
    plt.savefig(output, dpi=600)
    plt.show()
    plt.close()
    
def combined_density_plot(X, Y, Y2, Z, output):
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.contourf(X, Y, Z, cmap="gray")
    ax1.set_xlabel("X Coordinate (" r'$\AA$' ")", fontsize=15)
    ax1.set_ylabel("Y Coordinate (" r'$\AA$' ")", fontsize=15)
    ax1.set_xlim([np.amin(X), np.amax(X)]) 
    ax1.tick_params(labelsize=12)
    ax2.plot(X, Y2, color="white")
    ax2.set_ylabel("Number Density", fontsize=15)
    ax2.tick_params(labelsize=12)
    plt.savefig("Combined_Density.png", dpi=600)
    plt.show()