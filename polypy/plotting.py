"""
Plotting functions included with `polypy`.
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import numpy as np
import matplotlib.pyplot as plt
from polypy import fig_params

def line_plot(x, y, xlab, ylab, 
              figsize=(4, 3)):
    """
    Simple line plotting function. Designed to be generic and used in several different applications.

    Args:
        x (:py:attr:`array like`): x axis points.
        y (:py:attr:`array like`): y axis points.
        xlab (:py:attr:`str`): x axis label.
        ylab (:py:attr:`str`): y axis label.
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    ax = plt.subplots(figsize=figsize)[1]
    ax.plot(x, y)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.tick_params()
    plt.tight_layout()
    return ax

def msd_plot(msd_data, 
            show_all_dimensions=True,
            figsize=(4, 3)):
    """
    Plotting function for the mean squared displacements (MSD).

    Args:
        msd_data ():py:class:`polypy.msd.MSDContainer`): MSD data.
        show_all_dimensions (:py:attr:`bool`): Display all MSD data or the total MSD. Default is :py:attr:`bool`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    ax = plt.subplots(figsize=figsize)[1]
    ax.set_ylim(ymin=0, ymax=np.amax(msd_data.msd))
    ax.set_xlim(xmin=0, xmax=np.amax(msd_data.time))
    ax.plot(msd_data.time, msd_data.msd, label="XYZMSD")
    if show_all_dimensions:
        ax.plot(msd_data.time, msd_data.xymsd, label="XYMSD")
        ax.plot(msd_data.time, msd_data.xzmsd, label="XZMSD")
        ax.plot(msd_data.time, msd_data.yzmsd, label="YZMSD")
        ax.plot(msd_data.time, msd_data.xmsd, label="XMSD")
        ax.plot(msd_data.time, msd_data.ymsd, label="YMSD")
        ax.plot(msd_data.time, msd_data.zmsd, label="ZMSD")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("MSD ($\AA$)")
    ax.legend(frameon=True, edgecolor="black", loc=2)
    return ax


def volume_plot(x, y, 
                xlab="Timestep (ps)", 
                ylab="System Volume ($\AA$)",
                figsize=(4, 3)):
    """
    Gathers the data and creates a line plot for the system volume as a function of simulation timesteps

    Args:
        x (:py:attr:`array like`): x axis points - simulation timesteps
        y (:py:attr:`array like`): y axis points - Volume
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"Timestep (ps)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"System Volume ($\AA$)"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    line_plot(x, y, xlab, ylab, figsize)


def electric_field_plot(x, y,
                        xlab="X Coordinate ($\AA$)",
                        ylab="Electric Field (V)",
                        figsize=(4, 3)):
    """
    Gathers the data and creates a line plot for the electric field in one dimension.

    Args:
        x (:py:attr:`array like`): x axis points - position in simulation cell
        y (:py:attr:`array like`): y axis points - electric field
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Electric Field (V)"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    line_plot(x, y, xlab, ylab, figsize)


def electrostatic_potential_plot(x, y,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Electrostatic Potential (V)",
                                 figsize=(4, 3)):
    """
    Gathers the data and creates a line plot for the electrostatic potential in one dimension.

    Args:
        x (:py:attr:`array like`): x axis points - position in simulation cell
        y (:py:attr:`array like`): y axis points - electrostatic potential
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Electrostatic Potential (V)"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    line_plot(x, y, xlab, ylab, figsize)


def one_dimensional_charge_density_plot(x, y,
                                        xlab="X Coordinate ($\AA$)",
                                        ylab="Charge Density",
                                        figsize=(4, 3)):
    """
    Gathers the data and creates a line plot for the charge density in one dimension.

    Args:
        x (:py:attr:`array like`): x axis points - position in simulation cell
        y (:py:attr:`array like`): y axis points - charge density
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Charge Density"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    line_plot(x, y, xlab, ylab, figsize)


def one_dimensional_density_plot(x, y, data_labels,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Particle Density",
                                 figsize=(4, 3)):
    """
    Plots the number density of all given species in one dimension. 

    Args:
        x (:py:attr:`list`): x axis points - list of numpy arrays containing x axis coordinates.
        y (:py:attr:`list`): y axis points - list of numpy arrays containing y axis coordinates.
        data_labels (:py:attr:`list`): List of labels for legend.
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Particle Density"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    ax = plt.subplots(figsize=figsize)[1]
    for i in range(len(x)):
        ax.plot(x[i], y[i], label=data_labels[i])
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.tick_params()
    ax.legend(frameon=True, edgecolor="black")
    plt.tight_layout()
    return ax

def two_dimensional_charge_density_plot(x, y, z, 
                                        xlab="X Coordinate ($\AA$)",
                                        ylab="Y Coordinate ($\AA$)",
                                        palette="viridis",
                                        figsize=None,
                                        colorbar=True):
    """
    Plots the charge density in two dimensions. 

    Args:
        x (:py:attr:`array like`): x axis points -  x axis coordinates.
        y (:py:attr:`array like`): y axis points -  y axis coordinates.
        z (:py:attr:`array like`): z axis points -  2D array of points.
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Y Coordinate ($\AA$)"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.
        colorbar (:py:class:`bool`): Include the colorbar or not.


    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    fig, ax = plt.subplots(figsize=figsize)
    CM = ax.contourf(x, y, z, cmap=palette)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.tick_params()
    if colorbar:
        cbar = fig.colorbar(CM)
        cbar.set_label('Charge Density', labelpad=-40, y=1.07, rotation=0)
    return fig, ax

def two_dimensional_density_plot(x, y, z,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Y Coordinate ($\AA$)",
                                 palette="viridis",
                                 figsize=None, 
                                 colorbar=True):
    """
    Plots the distribution of an atom species in two dimensions. 

    Args:
        x (:py:attr:`array like`): x axis points -  x axis coordinates.
        y (:py:attr:`array like`): y axis points -  y axis coordinates.
        z (:py:attr:`array like`): z axis points -  2D array of points.
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Y Coordinate ($\AA$)"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.
        colorbar (:py:class:`bool`): Include the colorbar or not.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    fig, ax = plt.subplots(figsize=figsize)
    CM = ax.contourf(x, y, z, cmap=palette)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.tick_params()
    if colorbar:
        cbar = fig.colorbar(CM)
        cbar.set_label('Particle Density', labelpad=-40, y=1.1, rotation=0)
    plt.tight_layout()
    return fig, ax

def combined_density_plot(x, y, z,
                          xlab="X Coordinate ($\AA$)",
                          ylab="Y Coordinate ($\AA$)",
                          y2_lab="Particle Density",
                          palette="viridis",
                          figsize=None):
    """
    Plots the distribution of an atom species in two dimensions. 

    Args:
        x (:py:attr:`array like`): x axis points -  x axis coordinates.
        y (:py:attr:`array like`): y axis points -  y axis coordinates.
        z (:py:attr:`array like`): z axis points -  2D array of points.
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Y Coordinate ($\AA$)"`
        y2_lab (:py:attr:`str`): second y axis label. Default is :py:attr:`"Particle Density"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    y2 = np.sum(z, axis=1)
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.contourf(x, y, z, cmap=palette)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax1.set_xlim([np.amin(x), np.amax(x)])
    ax1.tick_params()
    ax2.plot(x, y2)
    ax2.set_ylabel(y2_lab)
    ax2.tick_params()
    plt.tight_layout()
    return fig, ax1
