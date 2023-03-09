# -*- coding: utf-8 -*-

import os
import Process.ClassTraj as t
import Process.TrajProcessFunc as func
from Analysis.mainFunc import *
from Config.ConfigInput import *
import numpy as np
import matplotlib as mpl
from Process.ClassFreq import *
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

mpl.rc('font',family='Times New Roman')
plt.rcParams.update({'font.size': 12})
output = open('out.log', 'w')
confDict.print_config()
plotDict.print_config()
tot_trajList = []
last_peak_array = np.array([])
which_dirArr = np.array([])
filenameList = []
# directTraj = '/Users/zhtfeng/Desktop/Migratory Aptitude/Uphill Dynamics/Secondary/cppac/'
# directTraj = 'D:\\Projects\\Nopsane Project\\MD\\waterleaving-downhill\\
directTraj = 'D:\\Projects\\ShmidtReaction\\Uphill\\'
# directTraj = 'D:\\Projects\\Migratory Aptitude\\Uphill Dynamics\\Secondary\\vih-solvent\\'
# directTraj = 'D:\\Port\\'
# =============================================================================
# Read the def.txt file of the definition of the structures
if confDict.config['Analyze']: anaRules = func.read_ana_file(directTraj)
# =============================================================================

# =============================================================================
# Reading trajectories at current folder and catagorize it if necessary

if confDict.config['Traj']:

    index = 0
    trajname_list = []
    direction = confDict.config['Direction']

    if directTraj is not 'Specified':

        directory_info_list = list(os.walk(directTraj))

    else:

        directory_info_list = list(os.walk(os.getcwd()))

    for i in range(0, len(directory_info_list[0][1]) + 1):

        for filename in directory_info_list[i][2]:

            if filename[0:4] == 'traj' and len(filename) > 4:

                if direction == "Both":
                    dyn_coord = func.get_raw_dyn_downhill(directory_info_list[i][0] + '/' + filename)
                    traj = t.Downhill_trajectory(dyn_coord[0], dyn_coord[1], filename, index, dyn_coord[2],
                                                 dyn_coord[3])
                if direction == "Only":

                    # filenameList.append(filename)
                    # which_dirArr=np.append(which_dirArr,int(directory_info_list[i][0][-1]))
                    dyn_coord = func.get_raw_dyn_uphill(directory_info_list[i][0] + '/' + filename)
                    traj = t.Uphill_trajectory(dyn_coord[0], dyn_coord[1], filename, index, dyn_coord[2], True)
                    # last_peak = signal.argrelextrema(np.array((traj.distance([6, 8]))), np.greater)[0][-1]
                    # last_peak_array = np.append(last_peak_array,last_peak)


                if traj.complete:

                    if confDict.config['Analyze']: func.ana_traj_file(anaRules[0], traj)
                    if confDict.config['Analyze']: func.assign_traj(traj, anaRules[1], anaRules[2])

                tot_trajList.append(traj)
                index += 1

# =============================================================================


# =============================================================================
#  pull out the trajectory of interest
if confDict.config['Analyze']:

    analyze_all_traj(tot_trajList, directTraj)
    refined_traj = specify_traj_for_analyze(tot_trajList, mode=confDict.config['TrajOfInterest'],
                                            totalNum=confDict.config['Number'],
                                            random=confDict.config['Random'])

    if confDict.config['Random']:
        Indexpoll = np.arange(len(refined_traj))
        np.random.shuffle(Indexpoll)
        refined_traj = [refined_traj[Indexpoll[0]]]


else:
    refined_traj = tot_trajList
    write_to_output(['No analyze file found, all trajectories were analyzed instead'])

# =============================================================================
# =============================================================================
# This Section deals will plot the surfaces

if confDict.config['Geometry_vs_time']: plot_geom_vs_time(refined_traj, geom_para1)
if confDict.config['Geom_vs_Geom_plot']: plot_geometrical_changes(refined_traj, [geom_para1, geom_para2],
                                                                  zaxis=geom_para3)
if confDict.config['Timing_histogram']: plot_timing_histogram(refined_traj, mode='Together', box_width=5)
if confDict.config['Plot_surface']:
    plot_surface(pesfilename, stepsize, plotType, None, xyrange=[xrange, yrange])
if confDict.config['Trajecory_on_surface']: plot_surface(pesfilename, stepsize, plotType, refined_traj,
                                                         xyrange=[xrange, yrange], geom=[geomx, geomy])
# if confDict.config['Group_kinetic_energy_analysis']:
if confDict.config['Plot_energy_vs_time']: plot_energy_vs_time(refined_traj)
if confDict.config['Count_number_of_peaks']: count_number_of_peaks(refined_traj)
if confDict.config['Linear_kinetic_energy']: linear_kinetic_check(refined_traj, atoms=[3, 5], cutoff=1.5, mass=45)
# if confDict.config['3D_Bond_Energy_Time_Plot']: pass
if confDict.config['Customize']:


    fig, ax = plt.subplots(figsize=(12,8),dpi=400)
    for traj in refined_traj:
        # ax.plot(range(traj.time),traj.distance([5,7]))
        ax.plot(range(traj.time), traj.distance([5, 7]))
    ax.set_xlabel('Time (fs)',fontsize=32,fontname='Times New Roman')
    ax.set_ylabel('Bond Length (Å)',fontsize=32,fontname='Times New Roman')
    plt.savefig('D:\\Projects\\Migratory Aptitude\\waterNopsane22')
    plt.show()

# =============================================================================

write_to_output(['Normal Termination'])
output.close()