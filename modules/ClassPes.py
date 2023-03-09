#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 06:12:04 2019

@author: zhtfeng
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from Config.ConfigInput import *


class PES:

    def __init__(self, filename, xrange, yrange):

        self.filename = filename
        self.rawdata = pd.read_csv(filename)
        self.rawdata.values
        self.minpoint = np.nanmin(self.rawdata)
        self.data = (self.rawdata-self.minpoint)*627.5
        self.data = np.array(self.data)
        self.xmin, self.xmax = xrange
        self.ymin, self.ymax = yrange
        self.row, self.col = np.shape(self.data)

    def plot_pes(self, stepsize=None):

        if stepsize == None:

            pass

        else:

            X = np.linspace(self.xmin, self.xmax, num=self.col)
            Y = np.linspace(self.ymin, self.ymax, num=self.row)

        X, Y = np.meshgrid(X, Y)
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(12,12),dpi=240)
        # fig, ax = plt.subplots(figsize=(12,12),dpi=240,)
        # plt.rc('xtick', size=50)
        # surf = ax.contourf(X, Y, self.data, cmap=cm.coolwarm,linewidth=1)
        masked_array = np.ma.array (self.data, mask=np.isnan(self.data))
        # cmap=cm.coolwarm
        # cmap.set_bad('white',1.)
        surf = ax.plot_surface(X, Y, self.data,cmap=cm.jet, alpha=1.,vmin=np.nanmin(self.data), vmax=np.nanmax(self.data),linewidth=1,antialiased=False)
        ax.set_title(surface_title)
        ax.set_xlim(min(self.xmin,self.xmax),max(self.xmin,self.xmax))
        ax.set_ylim(min(self.ymin,self.ymax),max(self.ymin,self.ymax))
        ax.set_zlim(np.nanmin(self.data),np.nanmax(self.data))
        cbar = fig.colorbar(surf,ax=ax,extend='both')
        cbar.minorticks_on()
        ax.view_init(20,60)

#        ax.annotate('acyl-TSS', color='r',
#                    xy=(2.35, 2.49), xycoords='data',
#                    xytext=(.55, .65), textcoords='axes fraction',
#                    arrowprops=dict(facecolor='r', shrink=0.06),
#                    horizontalalignment='right', verticalalignment='top')
#        ax.annotate('acyl-VTS', color='b',
#                    xy=(2.43, 2.44), xycoords='data',
#                    xytext=(.59, .58), textcoords='axes fraction',
#                    arrowprops=dict(facecolor='b', shrink=0.06),
#                    horizontalalignment='right', verticalalignment='top')
#        ax.annotate('iPr-TSS',color='crimson',
#                    xy=(2.51, 2.25), xycoords='data',
#                    xytext=(.95, .50), textcoords='axes fraction',
#                    arrowprops=dict(facecolor='crimson', shrink=0.06),
#                    horizontalalignment='right', verticalalignment='top')
#        ax.annotate('iPr-VTS',color='navy',
#                    xy=(2.51, 2.26), xycoords='data',
#                    xytext=(.95, .56), textcoords='axes fraction',
#                    arrowprops=dict(facecolor='navy', shrink=0.06),
#                    horizontalalignment='right', verticalalignment='top')

        # ax.text(2.33,2.51,'Al',fontsize=15)
        # ax.text(2.52,2.28,'Et', fontsize=15)
        return fig, ax

    def plot_gradient(self, stepsize=None):

        if stepsize == None:

            pass

        else:

            X = np.linspace(self.xmin, self.xmax, num=self.col)
            Y = np.linspace(self.ymin, self.ymax, num=self.row)

        X, Y = np.meshgrid(X, Y)
        grad = np.gradient(self.data)
        fig, ax = plt.subplots()
        M = np.hypot(grad[0], grad[1])
        Q = ax.quiver(X, Y, -grad[1], grad[0], M,
                      pivot='mid', units='inches')
        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        ax.set_title(surface_title)
        ax.set_xlim(min(self.xmin, self.xmax), max(self.xmin, self.xmax))
        ax.set_ylim(min(self.ymin, self.ymax), max(self.ymin, self.ymax))
        ax.text(2.33, 2.51, 'Al', fontsize=15)
        ax.text(2.52, 2.28, 'Et', fontsize=15)

        return fig, ax
