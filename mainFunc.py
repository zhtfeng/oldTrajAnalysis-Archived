#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 20:24:40 2019
@author: zhtfeng

"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import Process.ClassPes as pes
import os
import scipy.signal as signal
import Config.ClassConfig as conf
from Config.ConfigInput import *


def analyze_all_traj(allTraj,directTraj):

    nameList = [i.end_point_name for i in allTraj]
    names = []
    for i in nameList:

        if i not in names:
            names.append(i)

    for each_set in names:

        if each_set not in names:
            trajtype_set.append(each_set)

    if directTraj == 'Specified':
        write_to_output('Total {} Trajectories were found in folder {}'.format(len(nameList), os.getcwd()))
        print('Total {} Trajectories were found in folder {}'.format(len(nameList), os.getcwd()))
    else:
        write_to_output('Total {} Trajectories were found in folder {}'.format(len(nameList), directTraj))
        print('Total {} Trajectories were found in folder {}'.format(len(nameList), directTraj))

    for eachtype in names:

        number_of_traj = nameList.count(eachtype)

        if len(eachtype) == 0:

            pass
        elif len(eachtype) == 1:

            eachtype = list(eachtype)
            write_to_output(['There are total ', str(number_of_traj), str(eachtype[0]), ' to ', str(eachtype[0]),
                             'Trajectories, ',
                             str(float(number_of_traj * 100 / len(allTraj))), ' % of total'])
            print('There are total', number_of_traj, eachtype[0], ' to ', eachtype[0], 'Trajectories, ', \
                  float(number_of_traj * 100 / len(allTraj)), '% of total')

        else:

            eachtype = list(eachtype)

            write_to_output(
                ['There are total ', str(number_of_traj), str(eachtype[0]), ' to ', str(eachtype[1]), ' Trajectories, ', \
                 str(float(number_of_traj * 100 / len(allTraj))), ' % of total'])
            print('There are total ', number_of_traj, eachtype[0], 'to', eachtype[1], 'Trajectories, ', \
                  float(number_of_traj * 100 / len(allTraj)), '% of total')


def specify_traj_for_analyze(allTrajArr, mode='All', totalNum=False, random=False,whichtraj='all'):
    lucky_trajList = []

    if mode == 'All':

        lucky_trajList = allTrajArr

    else:

        for traj in allTrajArr:

            if mode == 'All': lucky_trajList = allTrajArr
            if mode == 'Defined':

                if len(traj.end_point_name) != 1:
                    lucky_trajList.append(traj)

            if mode not in ['All', 'Defined']:

                if mode in traj.end_point_name: lucky_trajList.append(traj)

    if totalNum > len(lucky_trajList):
        totalNum = len(lucky_trajList)
        write_to_output(['Number of trajectory specified exceeds the total number, use all instead \n'])

    if type(totalNum) == int:
        lucky_trajList = lucky_trajList[:totalNum+1]

    if isinstance(whichtraj,int):  # Specify which trajectory or which set of trajectories to plot. controlled by whichtraj variable

        if whichtraj > len(lucky_trajList):
            print('trajectory ID exceeds total traj number, use the last one instead')
            return [lucky_trajList[-1]]
        else: return [lucky_trajList[whichtraj]]
    elif isinstance(whichtraj, list):

        output_traj_list=[]
        for n in whichtraj:
            if n > len(lucky_trajList):
                print('trajectory ID exceeds total traj number, skipped')
            else:
                output_traj_list.append(lucky_trajList[n])
        if len(output_traj_list) == 0:

            print('No traj specified for next step, use all instead')
            return allTrajArr
        else:
            return output_traj_list


    else:

        return lucky_trajList





def para_dim_check(trajArr, dim):
    if len(dim) == 2:
        paraForPlot = np.array([traj.distance(dim) for traj in trajArr])
        paraType = ' Bond Length '

    elif len(dim) == 3:
        paraForPlot = np.array([traj.angle(dim) for traj in trajArr])
        paraType = ' Bond Angle '

    elif len(dim) == 4:
        paraForPlot = np.array([traj.dihedral(dim) for traj in trajArr])
        paraType = ' Dihedral Angle '

    return paraForPlot, paraType


def plot_geom_vs_time(trajArr, geom_para, newplot=True, plot_type='Time'):
    geom_para = np.array(geom_para)


    if newplot is True:

        fig, ax = plt.subplots(figsize=(12,8),dpi=160)

    else:

        fig, ax = newplot

    if plot_type == 'Time':

        for traj_para in para_dim_check(trajArr, geom_para)[0]:
            ax.plot(np.arange(len(traj_para)), traj_para)

        ax.set_xlabel('Time (fs)')
        ax.set_ylabel('C-H Bond Length (Å)')
        ax.set_title(geom_vs_time_title)
        write_to_output(
            ['A timing and ', para_dim_check(trajArr, geom_para)[1], str(geom_para), ' diagram has been plotted'])
        plt.axhline(1.67,color='black')
        plt.show()


    return fig, ax


def plot_geometrical_changes(trajArr, geom_para, newPlot=None, zaxis=False):
    geom_x, geom_y = geom_para
    geomxArr, geomType_x = para_dim_check(trajArr, geom_x)
    geomyArr, geomType_y = para_dim_check(trajArr, geom_y)

    if isinstance(zaxis, list):
        geomz = zaxis
        geomzArr, geomType_z = para_dim_check(trajArr, geomz)

    if newPlot is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = newPlot

    for i in range(len(geomxArr)):
        if isinstance(zaxis, list):
            ax = fig.gca(projection='3d')
            ax.plot(geomxArr[i], geomyArr[i], geomzArr[i])

        else:
            ax.plot(geomxArr[i], geomyArr[i])

    ax.set_title(geom_vs_geom_title)
    ax.set_xlabel(x_axis)
    ax.set_ylabel(y_axis)
    ax.plot(2.35, 2.49, 'or', color="red")
    ax.plot(2.43, 2.44, 'or', color="b")
    ax.plot(2.51, 2.25,'or', color="crimson")
    ax.plot(2.51, 2.26, 'or', color="navy")
    plt.savefig('pes.png')

    write_to_output(
        ['A plot of ', str(geomType_x), str(geom_x), ' and ', str(geomType_y), str(geom_y), ' has been plotted. '])
    plt.show()
    return fig, ax


def plot_timing_histogram(trajArr, mode='Together', box_width=5, names=[]):
    if mode == 'Together':

        timeArr = np.array([traj.time for traj in trajArr])
        fig, ax = plt.subplots()
        ax.set_ylabel('Probability density')
        ax.set_xlabel('Time/fs')
        ax.hist(timeArr, box_width, density=True, histtype='bar', facecolor='dodgerblue', rwidth=0.6)
        plt.show()
        ax.set_title(timing_histogram_title)

    write_to_output(['A Timing Histogram has been plotted'])
    return fig, ax


def plot_surface(filename, stepsize, surface, trajArr=None, xyrange=[[0, 1], [0, 1]], geom=None):
    xrange, yrange = xyrange
    if trajArr is None:

        if stepsize is not None:
            if surface == 'PES':
                pes.PES(filename, xrange, yrange).plot_pes(stepsize)
            if surface == 'Gradient':
                pes.PES(filename, xrange, yrange).plot_gradient(stepsize)
        plt.show()

    else:

        if stepsize is not None:

            if surface == 'PES':

                plotPES = pes.PES(filename, xrange, yrange).plot_pes(stepsize)
                plot_geometrical_changes(trajArr, geom, newPlot=plotPES)

            elif surface == 'Gradient':

                plotGradient = pes.PES(filename, xrange, yrange).plot_gradient(stepsize)
                plot_geometrical_changes(trajArr, geom, newPlot=plotGradient)

            plt.show()


def plot_energy_vs_time(trajArr):
    fig, ax = plt.subplots()

    for traj in trajArr:
        ax.plot(range(traj.time), traj.energy)

    plt.show()
    return fig, ax


def count_number_of_peaks(trajArr):
    num_of_peakArr = np.array([])

    for traj in trajArr:
        num_of_peakArr = np.append(num_of_peakArr,
                                   len(signal.argrelextrema(np.array((traj.distance([6, 8]))), np.greater)[0]))

    write_to_output(['The average number of peaks of the selected trajectories is ', str(np.average(num_of_peakArr))])
    print('The average number of peaks of the selected trajectories is ', str(np.average(num_of_peakArr)))

    return num_of_peakArr


def write_to_output(inputline):
    filename = 'out.log'
    output = open(filename, 'a')
    inputline = list(inputline)
    inputline.append('\n')
    output.writelines(inputline)

def linear_kinetic_check(trajArr,atoms,cutoff,mass):

    timeArr =np.array([])
    for traj in trajArr:

        timeArr = np.append(timeArr,traj.truncate(atoms[0],atoms[1],cutoff)**2*mass)
        print(traj.truncate(atoms[0],atoms[1],cutoff)**2*mass)

    if len(timeArr) != 0 : print('Average kinetic Energy is : ', np.mean(timeArr))

def momentum_decomposition(trajArr,atomList):

    pass

def bond_Energy_Time_Plot(trajNum):

    pass

def kinetic_energy_direction_plot(trajArr,nameList,atom1,atom2,plane1,plane2,t=1):

    fig, ax = plt.subplots(figsize=(8,7))
    plt.rc('legend', fontsize=14)
    plt.rc('figure', titlesize=30)
    grp1,grp2 = nameList
    aArr = []
    bArr = []
    cArr = []
    dArr = []
    eArr = []
    fArr = []

    for traj in trajArr:

        if grp1 in traj.end_point_name:

            initalA = traj.velTheta(atom1, plane1[0], plane1[1], plane1[2])[t]
            initalB = traj.velTheta(atom2, plane2[0], plane2[1], plane2[2])[t]
            aArr.append(initalA)
            bArr.append(initalB)

        elif grp2 in traj.end_point_name:

            initalA = traj.velTheta(atom1, plane1[0], plane1[1], plane1[2])[t]
            initalB = traj.velTheta(atom2, plane2[0], plane2[1], plane2[2])[t]
            cArr.append(initalA)
            dArr.append(initalB)

        else:

            initalA = traj.velTheta(atom1, plane1[0], plane1[1], plane1[2])[t]
            initalB = traj.velTheta(atom2, plane2[0], plane2[1], plane2[2])[t]
            eArr.append(initalA)
            fArr.append(initalB)
    ax.scatter(aArr, bArr, color='orangered',label = str(grp1+' Migration'))
    ax.scatter(cArr, dArr, color='limegreen', label = str(grp2+' Migration'))
    ax.scatter(eArr, fArr, color='indigo', label = 'Not reactive')
    ax.set_xlim(0, 370)
    ax.set_ylim(0, 370)
    plt.rc('font', size=12)
    ax.set_xlabel(str('Angle of '+ 'i-Pr'+ ' Momentum (degree)'))
    ax.set_ylabel(str('Angle of '+ 'Acyl'+ ' Momentum (degree)'))
    ax.set_title('Momentum Angle Plot of ' + '$R^1 = i-Pr, R^2 = Acyl$')
    ax.xaxis.label.set_fontsize(12)
    ax.yaxis.label.set_fontsize(12)
    ax.legend(loc='upper right')
    plt.show()

def kinetic_energy_Planedirection_plot(trajArr,nameList,atom1,atom2,plane1,plane2,t=1):

    fig, ax = plt.subplots()
    grp1,grp2 = nameList
    aArr = []
    bArr = []
    cArr = []
    dArr = []
    eArr = []
    fArr = []

    for traj in trajArr:

        if grp1 in traj.end_point_name:

            initalA = traj.velThetaPlane(atom1, plane1[0], plane1[1], plane1[2])[t]
            initalB = traj.velThetaPlane(atom2, plane2[0], plane2[1], plane2[2])[t]

            aArr.append(initalA)
            bArr.append(initalB)

        elif grp2 in traj.end_point_name:

            initalA = traj.velThetaPlane(atom1, plane1[0], plane1[1], plane1[2])[t]
            initalB = traj.velThetaPlane(atom2, plane2[0], plane2[1], plane2[2])[t]
            cArr.append(initalA)
            dArr.append(initalB)

        else:

            initalA = traj.velThetaPlane(atom1, plane1[0], plane1[1], plane1[2])[t]
            initalB = traj.velThetaPlane(atom2, plane2[0], plane2[1], plane2[2])[t]
            eArr.append(initalA)
            fArr.append(initalB)

    ax.scatter(aArr, bArr, color='orangered',label = str(grp1+' Migration'))
    ax.scatter(cArr, dArr, color='limegreen', label = str(grp2+' Migration'))
    ax.scatter(eArr, fArr, color='indigo', label = 'Not reactive')

    ax.set_xlabel(str('Angle of '+ grp1+ ' Momentum (degree)'))
    ax.set_ylabel(str('Angle of '+ grp2+ ' Momentum (degree)'))

    ax.legend()
    plt.show()

def kinetic_energy_magnitude_plot(trajArr,nameList,atom1,atom2,vecA,vecB,t=1):

    fig, ax = plt.subplots()
    grp1,grp2 = nameList
    aArr = []
    bArr = []
    cArr = []
    dArr = []
    eArr = []
    fArr = []
    kin1 = []
    kin2 = []
    for traj in trajArr:

        if grp1 in traj.end_point_name:

            initalA = 0
            initalB = 0
            for x in atom1:
                initalA += traj.atom_velPrj(x,vecA[0],vecA[1])[t]
            for y in atom2:
                initalB += traj.atom_velPrj(y,vecB[0],vecB[1])[t]
            kin1.append(abs(initalA))
            kin2.append(abs(initalB))
            aArr.append(initalA)
            bArr.append(initalB)

        elif grp2 in traj.end_point_name:

            initalA = 0
            initalB = 0
            for x in atom1:
                initalA += traj.atom_velPrj(x,vecA[0],vecA[1])[t]
            for y in atom2:
                initalB += traj.atom_velPrj(y,vecB[0],vecB[1])[t]
            kin1.append(abs(initalA))
            kin2.append(abs(initalB))
            cArr.append(initalA)
            dArr.append(initalB)

        else:

            initalA = 0
            initalB = 0
            for x in atom1:
                initalA += traj.atom_velPrj(x,vecA[0],vecA[1])[t]
            for y in atom2:
                initalB += traj.atom_velPrj(y,vecB[0],vecB[1])[t]
            kin1.append(abs(initalA))
            kin2.append(abs(initalB))
            eArr.append(initalA)
            fArr.append(initalB)


    ax.scatter(aArr, bArr, color='tomato', label=str(grp1 + ' Migration'))
    ax.scatter(cArr, dArr, color='forestgreen', label=str(grp2 + ' Migration'))
    ax.scatter(eArr, fArr, color='mediumblue', label='Not reactive')
    ax.set_xlabel(str('Velocity of ' + grp1 + ' (m/s )' + str(np.mean(kin1))))
    ax.set_ylabel(str('Velocity of ' + grp2 + ' (m/s )' + str(np.mean(kin2))))
    ax.legend()
    plt.show()

def avgSmooth(array,step): # a small tool to smooth oscillating curve by averaging every nth step

    rawlength = len(array)
    remain = rawlength % int(step)
    array = array[:-remain]
    newlength = len(array)
    sectionNum = int(newlength/step)

    avgArr = np.mean(array.reshape(sectionNum,step),axis=1)
    return avgArr

def exludeArr(total,excludeArr):

    outputArr = []
    for i in range(total):

        if i+1 not in excludeArr:

            outputArr.append(i+1)

    return outputArr






# def functional_group_KE_analysis(groups):
#
#    me_list = []
#    for i,group in enumerate(groups):
#
#        me_list.append(traj.func_grp_ke([26,25,27,4])[0])
#
