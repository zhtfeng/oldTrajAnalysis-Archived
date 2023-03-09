#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 07:45:39 2019

@author: zhtfeng
"""
import re
import numpy as np
def is_number_check(num):
    pattern = re.compile(r'^[-+]?[-0-9]\d*\.\d*|[-+]?\.?[0-9]\d*$')
    result = pattern.match(num)
    if result: return  True
    else: return False

def process(coord):
    ### we next do a sanity check on the second line ###

    initial_line = re.split(r' +',coord[0]) # initial line contains several important info about the point
    atom_list = []

    cartesian = np.zeros((len(coord) - 1, 3))



    for i in range(1,len(coord)):# extract out the coordinate of each atom and put it into a list
        
        line = re.split(r'[\s]',coord[i])
        atom_list.append(line[0])
        line = [x for x in line if x!= '']

        cartesian[i-1,0] = float(line[1])
        cartesian[i-1,1] = float(line[2])
        cartesian[i-1,2] = float(line[3])          
            
    return cartesian,initial_line,atom_list

#### this funtion works on the entire trajectory and will give you the list of cartesian coordinates
#### and edit the order of the points so that it makes a continuous trajectory of downhill dynamics
    
def get_raw_dyn_downhill(filename):


    read_file = open(str(filename))
    data = read_file.readlines()
    atm_num = int(data[0])
    temp = []
    traj_coord = []
    point_order = []
    total_coord_list,total_energy_list = [],[]
    completeness = True

    for line_num,data_line in enumerate(data):
        
        if line_num % (atm_num+2) == 0:
            
            traj_coord.append(temp)
            temp = []
            
        else:
            
            temp.append(data_line) 
            
    traj_coord.pop(0)
    cartesian_list = []
    energy_list = []
    
    for each in traj_coord:
        
        pro = process(each)
        atom_list = pro[2]
        cartesian_list.append(pro[0])
        energy_list.append(float(pro[1][0]))
        point_order.append(int(pro[1][-3]))
        
### now we edit the trajectory to make it in order####
    if 1 in point_order[1:]:
        
        split_point = point_order[1:].index(1)
        back_coord = cartesian_list[split_point+1:]
        back_energy = energy_list[split_point+1:]
        back_coord.reverse()
        back_energy.reverse()
        forward_coord = cartesian_list[0:split_point]
        forward_energy = energy_list[0:split_point]
        total_coord_list = back_coord + forward_coord
        total_energy_list = back_energy + forward_energy
        read_file.close()
        
    else:
        
        completeness = False

    return total_coord_list,total_energy_list,atom_list,completeness

#### similar function for uphill dynamics ####
    
def get_raw_dyn_uphill(filename):
    
    read_file = open(filename)
    data = read_file.readlines()
    atm_num = int(data[0])
    temp = []
    traj_coord = []

    #### check the third line to be coordinates or comments
    line = re.split(r'[\s]', data[2])
    line = [x for x in line if x != '']
    for each in line[1:]:
         if is_number_check(each) == False:
             startingLine = 3
             break
         else:
             startingLine = 2
    for line_num,data_line in enumerate(data):
        
        if line_num % (atm_num+startingLine) == 0:

            if startingLine == 3 and len(temp) > 1:
                temp.pop(1)
            traj_coord.append(temp)

            temp = []
            
        else:

            temp.append(data_line)

    traj_coord.pop(0)
    cartesian_list = []
    energy_list = []
    
    for each in traj_coord:
        
        pro = process(each)
        cartesian_list.append(pro[0])
        pro_1 = [x for x in pro[1] if x!= '']
        energy_list.append(float(pro_1[0]))
        
    atom_list = pro[2]   # occurence of errors if QCEIMS traj files are used !!!
    read_file.close()
    
    return cartesian_list,energy_list,atom_list

def read_ana_file(directTraj = 'Specified'):
    
    global struct_list
    if directTraj == 'Specified': ana_file = open('def.txt')
    else: ana_file = open(directTraj+'def.txt')
    ana = ana_file.readlines()
    print('The analyze file has defined ',ana.count('next\n')+1,' different outcomes')
    sub_command = [] # criteria of a single structure 
    sub_command_lines = [] # collection of structures
    struct_list,type_list = [],[]
    
    for line in ana:   
        
        if line[0:4] == 'next' or line[0:3] == 'end':  
            
            sub_command_lines.append(sub_command)
            sub_command = []
            
        else:     
            
            sub_command.append(line) 
        
    for sub_index,sub in enumerate(sub_command_lines):
        
        struct_list.append('Unknown Structure'+str(sub_index))
        type_list.append(None)
        
        for line in sub:
            
            line = line.split(' ')
            if line[0] == 'print':
                
                structure_name = line[1]
                struct_list[sub_index] = structure_name[:-1]
                
            if line[0] == 'type':
                
                type_list[sub_index] = line[-1]
                
    ana_file.close()
               
    return sub_command_lines,struct_list,type_list

def ana_traj_file(rules,traj):
   
    traj.initialize_attribute(full = len(rules))  #### Now we apply the criteria to the trajctories ####
    
    for sub_index,sub in enumerate(rules):

        for line in sub:
                     
            line = line.split(' ')

            if line[0] == 'bond':
               
                geom_parameter = traj.dist(int(line[1]),int(line[2]),-1)
                starting_geom_parameter = traj.dist(int(line[1]),int(line[2]),0)
               
            elif line[0] == 'angle':
               
                geom_parameter = traj.ang(int(line[1]),int(line[2]),int(line[3]),-1)
                starting_geom_parameter = traj.ang(int(line[1]),int(line[2]),int(line[3]),0)
                
            elif line[0] == 'dihedral':
                
                geom_parameter = traj.dih(int(line[1]),int(line[2]),int(line[3]),int(line[4]),-1)
                starting_geom_parameter = traj.dih(int(line[1]),int(line[2]),int(line[3]),int(line[4]),0)
                
            elif line[0] == 'time':
                
                geom_parameter = float(line[-1])
           
            if line[-2] == '<':
    
                if geom_parameter > float(line[-1]):
                       
                    traj.attribute[sub_index] = 0
                    
                if starting_geom_parameter > float(line[-1]):
                    
                    traj.starting_attribute[sub_index] = 0
                    
                       
            elif line[-2] == '>':
                
                if geom_parameter < float(line[-1]):
                       
                    traj.attribute[sub_index] = 0
                
                if starting_geom_parameter > float(line[-1]):
                    
                    traj.starting_attribute[sub_index] = 0
            
        
def assign_traj(traj,struc_list,type_list): #### put names on to the trajs, name obtained from print
    
    
    if 1 not in traj.attribute:
        
        traj.end_point_name.add('Unknown')
        
    if 1 not in traj.starting_attribute:
        
        traj.end_point_name.add('Unknown')
        
    for i,j in enumerate(traj.attribute):

        if j == 1:

            traj.end_point_name.add(struc_list[i])
            traj.end_point_type.add(type_list[i])
            
    for i,j in enumerate(traj.starting_attribute):
                
        if j == 1:
             
            traj.end_point_name.add(struc_list[i])
            traj.end_point_type.add(type_list[i])

        traj.end_point_name = set(traj.end_point_name)          
            
            



    
