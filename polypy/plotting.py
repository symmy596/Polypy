"""
Plotting functions included with `polypy`.
"""

# Copyright (c) Adam R. Symington
# Distributed under the terms of the MIT License
# author: Adam R. Symington

import numpy as np
import matplotlib.pyplot as plt
from polypy import fig_params
from matplotlib.gridspec import GridSpec
from matplotlib import ticker

def line_plot(x, y, xlab, ylab, 
              figsize=(10, 6)):
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
            figsize=(10, 6)):
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
                figsize=(10, 6)):
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
    ax = line_plot(x, y, xlab, ylab, figsize)
    return ax


def electric_field_plot(x, y,
                        xlab="X Coordinate ($\AA$)",
                        ylab="Electric Field (V)",
                        figsize=(10, 6)):
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
    ax = line_plot(x, y, xlab, ylab, figsize)
    return ax


def electrostatic_potential_plot(x, y,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Electrostatic Potential (V)",
                                 figsize=(10, 6)):
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
    ax = line_plot(x, y, xlab, ylab, figsize)
    return ax


def one_dimensional_charge_density_plot(x, y,
                                        xlab="X Coordinate ($\AA$)",
                                        ylab="Charge Density",
                                        figsize=(10, 6)):
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
    ax = line_plot(x, y, xlab, ylab, figsize)
    return ax


def one_dimensional_density_plot(x, y, data_labels,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Particle Density",
                                 figsize=(10, 6)):
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
                                        figsize=(10, 6),
                                        colorbar=True, log=False):
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
        (:py:class:`matplotlib.Fig`): Figure object
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.

    """
    fig, ax = plt.subplots(figsize=figsize)
    if log:
        CM = ax.contourf(x, y, z, cmap=palette, locator=ticker.LogLocator())
    else:
        CM = ax.contourf(x, y, z, cmap=palette)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.tick_params()
    if colorbar:
        cbar = fig.colorbar(CM)
        cbar.set_label('Charge Density', labelpad=-40, y=1.1, rotation=0)
    return fig, ax

def two_dimensional_density_plot(x, y, z,
                                 xlab="X Coordinate ($\AA$)",
                                 ylab="Y Coordinate ($\AA$)",
                                 palette="viridis",
                                 figsize=(10, 6), 
                                 colorbar=True, log=False):
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
        log (:py:class:`bool`): Log the z data or not? This can sometimes be useful but obviously one needs to be careful
        when drawing conclusions from the data. 

    Returns:
        (:py:class:`matplotlib.Fig`): Figure object
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    fig, ax = plt.subplots(figsize=figsize)
    if log:
        CM = ax.contourf(x, y, z, cmap=palette, locator=ticker.LogLocator())
    else:
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
                          y2_lab="Number Density",
                          palette="viridis", linecolor="black",
                          figsize=(10, 6), log=False):
    """
    Plots the distribution of an atom species in two dimensions and overlays the one 
    dimensional density on top. Think of it as a combination of the two_dimensional_density_plot
    and one_dimensional_density_plot functions

    Args:
        x (:py:attr:`array like`): x axis points -  x axis coordinates.
        y (:py:attr:`array like`): y axis points -  y axis coordinates.
        z (:py:attr:`array like`): z axis points -  2D array of points.
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Y Coordinate ($\AA$)"`
        y2_lab (:py:attr:`str`): second y axis label. Default is :py:attr:`"Particle Density"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.
        log (:py:class:`bool`): Log the z data or not? This can sometimes be useful but obviously one needs to be careful
        when drawing conclusions from the data. 

    Returns:
        (:py:class:`matplotlib.Fig`): Figure object
        (:py:attr:`list`): List of axes objects.
    """
    y2 = np.sum(z, axis=0)
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = GridSpec(5, 2, figure=fig)
    gs.update(wspace=0.025, hspace=0.05)
    ax2 = fig.add_subplot(gs[0,:])
    ax1 = fig.add_subplot(gs[1:, :])
    ax = [ax1, ax2]
    if log:
        ax1.contourf(x, y, z, cmap=palette, locator=ticker.LogLocator())
    else:
        ax1.contourf(x, y, z, cmap=palette)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax1.set_xlim([np.amin(x), np.amax(x)])
    ax1.tick_params()
    ax2.plot(x, y2, color=linecolor)
    ax2.set_xlim([np.amin(x), np.amax(x)])
    ax2.axis('off')
    plt.tight_layout()
    return fig, ax

def two_dimensional_density_plot_multiple_species(x_list, y_list, z_list, palette_list,
                          xlab="X Coordinate ($\AA$)",
                          ylab="Y Coordinate ($\AA$)",
                          y2_lab="Number Density",
                          figsize=(10, 6), log=False):
    """
    Plots the distribution of a list of atom species in two dimensions. Returns
    heatmaps for each species stacking on top of one another. This is limited
    to four species.

    Args:
        x_list (:py:attr:`array like`): x axis points -  x axis coordinates.
        y_list (:py:attr:`array like`): y axis points -  y axis coordinates.
        z_list (:py:attr:`array like`): z axis points -  2D array of points.
        palette_list (:py:attr:`array like`): Color palletes for each atom species. 
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Y Coordinate ($\AA$)"`
        y2_lab (:py:attr:`str`): second y axis label. Default is :py:attr:`"Particle Density"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.
        log (:py:class:`bool`): Log the z data or not? This can sometimes be useful but obviously one needs to be careful
        when drawing conclusions from the data. 

    Returns:
        (:py:class:`matplotlib.Fig`): Figure object
        (:py:class:`matplotlib.axes.Axes`): The axes with new plots.
    """
    fig, ax1 = plt.subplots(figsize=figsize)
    alphas = [1.0, 0.7, 0.5, 0.3]
    if log:
        for i in range(len(x_list)):
            ax1.contourf(x_list[i], y_list[i], z_list[i], cmap=palette_list[i], locator=ticker.LogLocator())

    else:
        for i in range(len(x_list)):
            ax1.contourf(x_list[i], y_list[i], z_list[i], cmap=palette_list[i], alpha=alphas[i])
     
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax1.tick_params()
    plt.tight_layout()
    return fig, ax1

def combined_density_plot_multiple_species(x_list, y_list, z_list, palette_list, label_list, color_list,
                          xlab="X Coordinate ($\AA$)",
                          ylab="Y Coordinate ($\AA$)",
                          figsize=(10, 6), log=False):
    """
    Plots the distribution of a list of atom species in two dimensions. Returns
    heatmaps for each species stacking on top of one another. It also plots the
    same density in one dimension on top of the heatmaps. 
    
    Args:
        x (:py:attr:`list`): x axis points -  x axis coordinates.
        y (:py:attr:`list`): y axis points -  y axis coordinates.
        z (:py:attr:`list`): z axis points -  2D array of points.
        palette_list (:py:attr:`list`): Color palletes for each atom species. 
        label_list (:py:attr:`list`): List of species labels.
        color_list (:py:attr:`list`): List of colors for one dimensional plot.
        xlab (:py:attr:`str`): x axis label. Default is :py:attr:`"X Coordinate ($\AA$)"`
        ylab (:py:attr:`str`): y axis label. Default is :py:attr:`"Y Coordinate ($\AA$)"`
        fig_size (:py:class:`tuple`): Horizontal and veritcal size for figure (in inches). Default is :py:attr:`(10, 6)`.

    Returns:
        (:py:class:`matplotlib.Fig`): Figure object
        (:py:attr:`list`): List of axes objects.
    """
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = GridSpec(5, 2, figure=fig)
    gs.update(wspace=0.025, hspace=0.05)
    ax2 = fig.add_subplot(gs[0,:])
    ax1 = fig.add_subplot(gs[1:, :])
    ax = [ax1, ax2]
    alphas = [1.0, 0.7, 0.5, 0.3]
    if log:
        for i in range(len(x_list)):
            ax1.contourf(x_list[i], y_list[i], z_list[i], cmap=palette_list[i], locator=ticker.LogLocator())

    else:
        for i in range(len(x_list)):
            ax1.contourf(x_list[i], y_list[i], z_list[i], cmap=palette_list[i], alpha=alphas[i])
     
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax1.set_xlim([np.amin(x_list[0]), np.amax(x_list[0])])
    ax1.tick_params()
    for i in range(len(x_list)):
        ax2.plot(x_list[i], np.sum(z_list[i], axis=0), label=label_list[i], color=color_list[i])
    ax2.axis('off')
    ax2.set_ylim(np.amin(z_list[0]), np.amax(np.sum(z_list[0], axis=0)) * 1.4)
    ax2.set_xlim(np.amin(x_list[0]), np.amax(x_list[0]))

    ax2.legend(loc=2, ncol=len(label_list), frameon=False, fontsize=12)
    plt.tight_layout()
    return fig, ax