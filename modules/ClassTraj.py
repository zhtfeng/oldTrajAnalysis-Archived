#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 07:47:57 2019

@author: zhtfeng
"""
######## the class traj contains all the atrributes of a trajctory,
######## and with functions we can pull out "any" property needed

import numpy as np

massDict = {'H': 1.0082, 'He': 4.0026, 'Li': 6.997, 'Be': 9.0122, 'B': 10.821, 'C': 12.012, 'N': 14.008, 'O': 16.0,
            'F': 18.998,
            'Ne': 20.180, 'Na': 22.99, 'Mg': 24.307, 'Al': 26.982, 'Si': 28.086, 'P': 30.974, 'S': 32.076, 'Cl': 35.457,
            'Ar': 39.963,
            'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
            'Fe': 55.845, 'Co': 58.933,
            'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'As': 74.922, 'Se': 78.971, 'Br': 79.907}


class Trajectory:

    def __init__(self, coord_list, energy_list, name, index, atom_list, completeness):

        self.coord = np.array(coord_list)
        self.energy = energy_list
        self.index = index
        self.name = name
        self.time = len(self.energy)
        self.attribute = [1]
        self.end_point_name = set()
        self.atom = atom_list
        self.starting_attribute = [1]
        self.plotting_status = False
        self.complete = completeness
        self.recrossing = False
        self.end_point_type = set()


    def distance(self, bond):

        a, b = bond
        distArr = np.array([])

        for point in range(self.time):

            dist = self.dist(a, b, point)
            distArr = np.append(distArr, dist)

        return distArr

    def angle(self, bondangle):

        a, b, c = bondangle
        angle_list = []

        for point in range(self.time):
            ang = self.ang(a, b, c, point)
            angle_list.append(ang)

        return angle_list

    def dihedral(self, dihedralangle):

        a, b, c, d = dihedralangle
        dihedralArr = np.array([])

        for point in range(self.time):
            dih = self.dih(a, b, c, d, point)
            dihedralArr = np.append(dihedralArr, dih)

        return dihedralArr

    def initialize_attribute(self, full):

        self.attribute = self.attribute * full
        self.starting_attribute = self.starting_attribute * full

        return self.attribute, self.starting_attribute

    def dist(self, a, b, point):

        a -= 1
        b -= 1
        cart = self.coord[int(point)]
        dist = np.sqrt((cart[a, 0] - cart[b, 0]) ** 2 + (cart[a, 1] - cart[b, 1]) ** 2 +
                       (cart[a, 2] - cart[b, 2]) ** 2)

        return dist

    def ang(self, a, b, c, point):

        a -= 1
        b -= 1
        c -= 1
        cart = self.coord[int(point)]
        vec_a = np.array([(cart[a, 0] - cart[b, 0]), (cart[a, 1] - cart[b, 1]), (cart[a, 2] - cart[b, 2])])
        vec_b = np.array([(cart[c, 0] - cart[b, 0]), (cart[c, 1] - cart[b, 1]), (cart[c, 2] - cart[b, 2])])
        cos_theta = vec_a.dot(vec_b) / (np.sqrt(vec_a.dot(vec_a) * vec_b.dot(vec_b)))
        theta = np.arccos(cos_theta)
        theta = 180 * theta / np.pi

        return theta

    def dih(self, a, b, c, d, point):

        a -= 1
        b -= 1
        c -= 1
        d -= 1
        cart = self.coord[int(point)]
        vec_a_1 = np.array([(cart[a, 0] - cart[b, 0]), (cart[a, 1] - cart[b, 1]), (cart[a, 2] - cart[b, 2])])
        vec_a_2 = np.array([(cart[c, 0] - cart[b, 0]), (cart[c, 1] - cart[b, 1]), (cart[c, 2] - cart[b, 2])])
        vec_b_1 = np.array([(cart[b, 0] - cart[c, 0]), (cart[b, 1] - cart[c, 1]), (cart[b, 2] - cart[c, 2])])
        vec_b_2 = np.array([(cart[d, 0] - cart[c, 0]), (cart[d, 1] - cart[c, 1]), (cart[d, 2] - cart[c, 2])])
        norm_a = np.cross(vec_a_1, vec_a_2)
        norm_b = np.cross(vec_b_1, vec_b_2)
        cos_theta = norm_a.dot(norm_b) / (np.sqrt(norm_a.dot(norm_a) * norm_b.dot(norm_b)))
        theta = np.arccos(cos_theta)
        theta = 180 * theta / np.pi

        return theta

    def generate_mass_list(self):

        mass_list = np.zeros((len(self.atom), 1))
        for i, each_atom in enumerate(self.atom):

            mass_list[i] = massDict[each_atom]

        return mass_list

    def velocity(self,n):

        n -= 1
        xposArr = np.array([pt[n,0] for pt in self.coord])
        yposArr = np.array([pt[n,1] for pt in self.coord])
        zposArr = np.array([pt[n,2] for pt in self.coord])

        xvel = 1e5*np.gradient(xposArr)
        yvel = 1e5*np.gradient(yposArr)
        zvel = 1e5*np.gradient(zposArr)

        return np.column_stack((xvel,yvel,zvel)) # veloctiy in m/s


    def funcgrp_keTot(self, func_list):

        '''
        Iuput: Atom Number you are interested in

        Output: Array of total kinetic energy of the functional group along the trajecotry

        '''


        massList = self.generate_mass_list()
        keArr = np.zeros((self.time))

        for atomNum in func_list: # Loop over all the atoms mentioned

            xvel,yvel,zvel = self.velocity(atomNum)

            mass = massList[atomNum-1]
            keArr += 0.5*mass*(xvel**2+yvel**2+zvel**2)
        print(self.end_point_name)

        return keArr

    def funcgrp_velTot(self, func_list):



        velArr = np.zeros((self.time))

        for atomNum in func_list:

            xvel,yvel,zvel = self.velocity(atomNum-1)
            velArr += np.sqrt(xvel**2+yvel**2+zvel**2)

        return velArr

    def atom_velPrj(self,atomNum,vecptA,vecptB,velAbs = True):

        xvel, yvel, zvel = self.velocity(atomNum - 1)
        vecptA -= 1
        vecptB -= 1
        if velAbs: vecPrjArr = np.zeros((self.time,1))  # the velocity magnitude array is only a 1D array
        else: vecPrjArr = np.zeros((self.time,3))

        for pt in range(self.time):

            tempVelArr = np.array([xvel[pt],yvel[pt],zvel[pt]])
            tempbasisArr = self.coord[pt,vecptB] - self.coord[pt,vecptA] # generate the vector for the bond that we want to project on
            vectorPrjcted = tempVelArr*(np.dot(tempbasisArr,tempVelArr)/np.dot(tempbasisArr,tempbasisArr))
            if velAbs: vectorPrjcted = np.dot(tempVelArr, tempbasisArr)/np.linalg.norm(tempbasisArr)
            vecPrjArr[pt] = vectorPrjcted*1e5

        return vecPrjArr
    def center_of_mass(self):

        massList = self.generate_mass_list()
        center_of_massArr = np.zeros((self.time,3))
        totalMass = np.sum(self.generate_mass_list())

        for i in range(self.time):

            massWeightedCoord = self.coord[i] * massList
            center_of_massArr[i] = np.sum(massWeightedCoord,axis=0)/totalMass

        return center_of_massArr

    def center_of_massPartial(self,n):

        massList = self.generate_mass_list()
        center_of_massArr = np.zeros((self.time,3))
        totalMass = np.sum([massList[i-1] for i in n])

        for i in range(self.time):

                massWeightedCoord = self.coord[i] * massList

                for atom in n:
                    center_of_massArr[i] += massWeightedCoord[atom-1]/totalMass


        return center_of_massArr

    def transitionalVel(self):

        certerOfMassArr = self.center_of_mass()

        xcenterPosArr = np.array([center[0] for center in certerOfMassArr])
        ycenterPosArr = np.array([center[1] for center in certerOfMassArr])
        zcenterPosArr = np.array([center[2] for center in certerOfMassArr])   # parse the coordinates into 3 different directions


        xtransVel = 1e5*np.gradient(xcenterPosArr)
        ytransVel = 1e5*np.gradient(ycenterPosArr)
        ztransVel = 1e5*np.gradient(zcenterPosArr) # first-order derivatives of displacement is velocity

        return xtransVel,ytransVel,ztransVel # unit of m/s

    def transitionalVelPartial(self,n):

        certerOfMassArr = self.center_of_massPartial(n)
        xcenterPosArr = np.array([center[0] for center in certerOfMassArr])
        ycenterPosArr = np.array([center[1] for center in certerOfMassArr])
        zcenterPosArr = np.array([center[2] for center in certerOfMassArr])   # parse the coordinates into 3 different directions

        xtransVel = 1e5*np.gradient(xcenterPosArr)
        ytransVel = 1e5*np.gradient(ycenterPosArr)
        ztransVel = 1e5*np.gradient(zcenterPosArr) # first-order derivatives of displacement is velocity


        return xtransVel,ytransVel,ztransVel # unit of m/s

    def transKE(self):

        transVel = self.transitionalVel()
        totalMass = np.sum(self.generate_mass_list())

        return 0.5*totalMass*(transVel[0]**2 + transVel[1]**2 + transVel[2] **2)/4.184e6

    def transKEPartial(self,n):

        massList = self.generate_mass_list()
        transVel = self.transitionalVelPartial(n)
        totalMass = np.sum([massList[i-1] for i in n])

        return 0.5*totalMass*(transVel[0]**2 + transVel[1]**2 + transVel[2] **2)/4.184e6

    def transKEFrag(self,n):

        massList = self.generate_mass_list()
        transVel = self.transitionalVel()
        totalMass = np.sum([massList[i - 1] for i in n])


        return 0.5*totalMass*(transVel[0]**2 + transVel[1]**2 + transVel[2] **2)/4.184e6

    def GeneratePrincipalInertia(self,t):

        def GenerateInertiaTensor(t,coordP):

            inertiaTensor = np.zeros((3,3),dtype='float64')
            massArr = self.generate_mass_list()

            for i in range(len(self.atom)):


                inertiaTensor[0,0] += massArr[i]*(coordP[i,1]**2 + coordP[i,2]**2)
                inertiaTensor[0,1] += -massArr[i]*coordP[i,1]*coordP[i,0]
                inertiaTensor[0,2] += -massArr[i]*coordP[i,2]*coordP[i,0]
                inertiaTensor[1,0] += -massArr[i]*coordP[i,0]*coordP[i,1]
                inertiaTensor[1,1] += massArr[i]*(coordP[i,2]**2 + coordP[i,0]**2)
                inertiaTensor[1,2] += -massArr[i]*coordP[i,2]*coordP[i,1]
                inertiaTensor[2,0] += -massArr[i]*coordP[i,0]*coordP[i,2]
                inertiaTensor[2,1] += -massArr[i]*coordP[i,1]*coordP[i,2]
                inertiaTensor[2,2] += massArr[i]*(coordP[i,1]**2 + coordP[i,0]**2)

            return inertiaTensor/1e23 # the unit is kg*m^2/mol

        #### reorient the coordinates and velocities to the center_of_mass ####
        reorientCoord = np.zeros((len(self.atom),3))
        reorientVel = np.zeros((len(self.atom),3))
        centerOfMass = self.center_of_mass()
        centerofMassVel = np.column_stack(
            (self.transitionalVel()[0], self.transitionalVel()[1], self.transitionalVel()[2]))


        for atom in range(len(self.atom)): # Loop over all atom at a time point and reorient

            reorientCoord[atom] = self.coord[t,atom] - centerOfMass[t] #subtract the coord of center of mass at that time
            reorientVel[atom] = self.velocity(atom+1)[t] - centerofMassVel[t]

        #### Generate Inertia Tensor of given time T####

        inertia = GenerateInertiaTensor(t, reorientCoord)
        eigval, rotationMatx = np.linalg.eig(inertia)
        rotateCoord = np.zeros_like(reorientCoord)
        rotateVel = np.zeros_like(reorientCoord)

        for atom in range(len(self.atom)):
            rotateCoord[atom] = 1e-10*np.dot(rotationMatx,reorientCoord[atom])
            rotateVel[atom] = np.dot(rotationMatx,reorientVel[atom])

        return eigval,rotateCoord,rotateVel # all units in SI

    def GeneratePrincipalInertiaPartial(self,t,n):

        def GenerateInertiaTensor(t,coordP,n):

            inertiaTensor = np.zeros((3,3))
            massArr = self.generate_mass_list()

            for counter,atomIndex in enumerate(n):

                atomIndex -= 1
                inertiaTensor[0,0] += massArr[atomIndex]*(coordP[counter,1]**2 + coordP[counter,2]**2)
                inertiaTensor[0,1] += -massArr[atomIndex]*coordP[counter,1]*coordP[counter,0]
                inertiaTensor[0,2] += -massArr[atomIndex]*coordP[counter,2]*coordP[counter,0]
                inertiaTensor[1,0] += -massArr[atomIndex]*coordP[counter,0]*coordP[counter,1]
                inertiaTensor[1,1] += massArr[atomIndex]*(coordP[counter,2]**2 + coordP[counter,0]**2)
                inertiaTensor[1,2] += -massArr[atomIndex]*coordP[counter,2]*coordP[counter,1]
                inertiaTensor[2,0] += -massArr[atomIndex]*coordP[counter,0]*coordP[counter,2]
                inertiaTensor[2,1] += -massArr[atomIndex]*coordP[counter,1]*coordP[counter,2]
                inertiaTensor[2,2] += massArr[atomIndex]*(coordP[counter,1]**2 + coordP[counter,0]**2)

            return inertiaTensor/1e23 # the unit is kg*m^2/mol

        #### reorient the coordinates and velocities to the center_of_mass ####
        reorientCoord = np.zeros((len(n),3))
        reorientVel = np.zeros((len(n),3))
        centerOfMass = self.center_of_massPartial(n)
        centerofMassVel = np.column_stack(
            (self.transitionalVelPartial(n)[0], self.transitionalVelPartial(n)[1], self.transitionalVelPartial(n)[2]))


        for counter,atomIndex in enumerate(n): # Loop over all atoms at a time point and reorient

            atomIndex -= 1
            reorientCoord[counter] = self.coord[t,atomIndex] - centerOfMass[t]
            #subtract the coord of center of mass at that time
            reorientVel[counter] = self.velocity(atomIndex+1)[t] - centerofMassVel[t]
        #### Generate Inertia Tensor of given time T####

        inertia = GenerateInertiaTensor(t, reorientCoord,n)
        eigval, rotationMatx = np.linalg.eig(inertia)
        rotateCoord = np.zeros_like(reorientCoord)
        rotateVel = np.zeros_like(reorientCoord)

        for atom in range(len(n)):

            rotateCoord[atom] = 1e-10*np.matmul(rotationMatx,reorientCoord[atom].T)
            rotateVel[atom] = np.matmul(rotationMatx,reorientVel[atom].T)
        return eigval,rotateCoord,rotateVel # all units in SI

    def GenerateTotalAngularMomentum(self,t):

        _,rotateCoord,rotateVel = self.GeneratePrincipalInertia(t)
        massArr = self.generate_mass_list()
        Lx = np.zeros((3))
        Ly = np.zeros((3))
        Lz = np.zeros((3))

        for i in range(len(self.atom)):
            rx,ry,rz = rotateCoord[i]
            vx,vy,vz = rotateVel[i]
            Lx += 1e-3*massArr[i] * np.cross(np.array([0, ry, rz]), np.array([0, vy, vz]))
            Ly += 1e-3*massArr[i] * np.cross(np.array([rx, 0, rz]), np.array([vx ,0, vz]))
            Lz += 1e-3*massArr[i] * np.cross(np.array([rx, ry, 0]), np.array([vx, vy, 0]))

        return Lx, Ly, Lz

    def GenerateTotalAngularMomentumPartial(self, t,n):

        _, rotateCoord, rotateVel = self.GeneratePrincipalInertiaPartial(t,n)
        massArr = self.generate_mass_list()
        Lx = np.zeros((3))
        Ly = np.zeros((3))
        Lz = np.zeros((3))

        for counter,atomIndex in enumerate(n):
            atomIndex -= 1
            rx, ry, rz = rotateCoord[counter]
            vx, vy, vz = rotateVel[counter]
            Lx += 1e-3 * massArr[atomIndex] * np.cross(np.array([0, ry, rz]), np.array([0, vy, vz]))
            Ly += 1e-3 * massArr[atomIndex] * np.cross(np.array([rx, 0, rz]), np.array([vx, 0, vz]))
            Lz += 1e-3 * massArr[atomIndex] * np.cross(np.array([rx, ry, 0]), np.array([vx, vy, 0]))


        return Lx, Ly, Lz

    def GenerateAngularMomentumFrag(self,t,n):

        _, rotateCoord, rotateVel = self.GeneratePrincipalInertia(t)
        massArr = self.generate_mass_list()
        Lx = np.zeros((3))
        Ly = np.zeros((3))
        Lz = np.zeros((3))

        for counter,atomIndex in enumerate(n):
            atomIndex -= 1
            rx, ry, rz = rotateCoord[counter]
            vx, vy, vz = rotateVel[counter]
            Lx += 1e-3 * massArr[atomIndex] * np.cross(np.array([0, ry, rz]), np.array([0, vy, vz]))
            Ly += 1e-3 * massArr[atomIndex] * np.cross(np.array([rx, 0, rz]), np.array([vx, 0, vz]))
            Lz += 1e-3 * massArr[atomIndex] * np.cross(np.array([rx, ry, 0]), np.array([vx, vy, 0]))


        return Lx, Ly, Lz


    def RotationalKE(self,t):

        Lx,Ly,Lz = self.GenerateTotalAngularMomentum(t)
        principalMoments, _, _ = self.GeneratePrincipalInertia(t)
        angVelX = Lx[0]/principalMoments[0]
        angVelY = Ly[1]/principalMoments[1]
        angVelZ = Lz[2]/principalMoments[2]

        Erotation = 0.5*(principalMoments[0]*angVelX**2 + principalMoments[1]*angVelY**2 + principalMoments[2]*angVelZ**2)/4184 # unit in kcal/mol

        return Erotation,angVelX, angVelY,angVelZ

    def RotationalKEPartial(self,t,n):

        Lx,Ly,Lz = self.GenerateTotalAngularMomentumPartial(t,n)
        principalMoments, _, _ = self.GeneratePrincipalInertiaPartial(t,n)
        angVelX, angVelY, angVelZ = 0.0,0.0,0.0

        if principalMoments[0] > 1e-48: angVelX = Lx[0]/principalMoments[0]
        if principalMoments[1] > 1e-48: angVelY = Ly[1]/principalMoments[1]
        if principalMoments[2] > 1e-48: angVelZ = Lz[2]/principalMoments[2]

        Erotation = 0.5*(principalMoments[0]*angVelX**2 + principalMoments[1]*angVelY**2 + principalMoments[2]*angVelZ**2)/4184 # unit in kcal/mol

        return Erotation,angVelX, angVelY,angVelZ

    def RotationalKEFrag(self,t,n):

        Lx,Ly,Lz = self.GenerateAngularMomentumFrag(t,n)
        principalMoments, _, _ = self.GeneratePrincipalInertia(t)
        angVelX = Lx[0]/principalMoments[0]
        angVelY = Ly[1]/principalMoments[1]
        angVelZ = Lz[2]/principalMoments[2]

        Erotation = 0.5 * (principalMoments[0] * angVelX ** 2 + principalMoments[1] * angVelY ** 2 + principalMoments[
            2] * angVelZ ** 2) / 4184  # unit in kcal/mol

        return Erotation,angVelX, angVelY,angVelZ



    def VibKE(self,t,vel=False):

        _, rotateCoord, rotateVel = self.GeneratePrincipalInertia(t)
        _, angVelX, angVelY,angVelZ = self.RotationalKE(t)
        velRot = np.cross(np.array([angVelX,angVelY,angVelZ]),rotateCoord)

        velVib = rotateVel - velRot

        if vel:
            return velVib

        else:
            massArr = self.generate_mass_list()
            Evib = 0
            for i in range(len(self.atom)):

                Evib += 0.5e-3 * massArr[i] * np.sum(velVib[i]**2)

            return Evib/4184

    def VibKEPartial(self,t,n):

        _, rotateCoord, rotateVel = self.GeneratePrincipalInertiaPartial(t,n)
        _, angVelX, angVelY,angVelZ = self.RotationalKEPartial(t,n)

        velRot = np.cross(np.array([angVelX,angVelY,angVelZ]),rotateCoord)
        velVib = rotateVel - velRot
        massArr = self.generate_mass_list()
        Evib = 0
        for counter,atomIndex in enumerate(n):

            atomIndex -= 1
            Evib += 0.5e-3 * massArr[atomIndex] * np.sum(velVib[counter]**2)

        return Evib/4184

    def VibKEFrag(self,t,n):

        _, rotateCoord, rotateVel = self.GeneratePrincipalInertia(t)
        _, angVelX, angVelY,angVelZ = self.RotationalKE(t)
        velRot = np.cross(np.array([angVelX,angVelY,angVelZ]),rotateCoord)

        velVib = rotateVel - velRot
        massArr = self.generate_mass_list()
        Evib = 0
        for atom in n:

            atom -= 1
            Evib += 0.5e-3 * massArr[atom] * np.sum(velVib[atom]**2)

        return Evib/4184



    def truncate(self, a, b, trunc):

        time = 0
        for pt in range(self.time):

            if self.dist(a, b, pt) > trunc:
                time += 1

        return time



class Uphill_trajectory(Trajectory):

    def traj_filter(self, traj_plot_rule, total_traj_amount, traj_plot_amount):

        if traj_plot_rule == 'All':

            self.plotting_status = True

        elif traj_plot_rule == 'Defined':

            if len(self.end_point_name) != 1:
                self.plotting_status = True

        elif type(traj_plot_rule) == str:

            if traj_plot_rule in self.end_point_name:
                self.plotting_status = True

        elif traj_plot_rule == 'Random':

            randomNums = np.arange(total_traj_amount)
            np.random.shuffle(randomNums)
            if self.index in randomNums[:int(traj_plot_amount)]: self.plotting_status = True


class Downhill_trajectory(Trajectory):

    def traj_filter(self, traj_plot_rule, total_traj, traj_plot_amount):

        if traj_plot_rule == 'All':

            self.plotting_status = True

        elif traj_plot_rule == 'Defined':

            if len(self.end_point_name) != 1:
                self.plotting_status = True

        elif traj_plot_rule == 'Not_recrossing':

            if len(self.end_point_name) == 2 and 'Unknown' not in self.end_point_name:
                self.plotting_status = True

        elif type(traj_plot_rule) == str:

            if traj_plot_rule in self.end_point_name:
                self.plotting_status = True

    def get_endpoint_type(self):

        for i in list(self.structname):
            self.end_point_type.add(i.type)

            print(self.end_point_type)

    def recrossing_check(self):

        pass
