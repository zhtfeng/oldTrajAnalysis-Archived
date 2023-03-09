from Process.ClassTraj import *
import numpy as np
import re

class Freq:

    def __init__(self,filename):

        self.freqArr, self.redMass, self.displaceMatx,self.forceArr,self.atomNum = self.read_freq_file(filename)

    def read_freq_file(self, filename):

        # Input filename

        # Output: freq - list of frequencies of all the modes
        # reduced mass - list of reduced mass of each mode
        # freq_matx - a [num_of_mode,num_of_atom,3] 3D array that contains the cartesian info
        # of all the vibration modes
        # force_list - list of all the force constants

        line_num_list = []
        freq_list, reduced_mass_list, force_list = [], [], []
        freq_file = open(filename)
        freq_lines = freq_file.readlines()
        matx_disorder = []

        for line_num, each_line in enumerate(freq_lines):

            re.split(r'[\s]', each_line)
            if 'Frequencies' in each_line:
                line_num_list.append(line_num)

        atom_num = int((line_num_list[1] - line_num_list[0] - 7) / 3)
        line_num_list = line_num_list[:(len(line_num_list) - atom_num + 2)]
        freq_matx = np.zeros((3 * atom_num - 6, atom_num, 3))
        raw_matx = np.zeros((3 * atom_num, 3 * atom_num - 6))
        which_freq = 0

        for i, each_line_num in enumerate(line_num_list):

            freqline = re.findall("-?\d+\.*\d*", freq_lines[each_line_num])
            reduced_mass_line = re.findall("-?\d+\.*\d*", freq_lines[each_line_num + 1])
            force_line = re.findall("-?\d+\.*\d*", freq_lines[each_line_num + 2])
            freq_list += list(map(float, freqline))
            reduced_mass_list += list(map(float, reduced_mass_line))
            force_list += list(map(float, force_line))

            for atom_linenum, matx_lines in enumerate(range(each_line_num + 5, each_line_num + 5 + 3 * atom_num)):

                matx_disorder = list(map(float, re.findall("-?\d+\.*\d*", freq_lines[matx_lines])[3:]))

                for item in range(len(matx_disorder)):
                    raw_matx[atom_linenum, item + which_freq] = matx_disorder[item]

            which_freq += len(matx_disorder)

        row_num, col_num = np.shape(raw_matx)

        freq_list = np.array(freq_list)
        for cols in range(col_num):

            for rows in range(row_num):

                if rows % 3 == 0:

                    freq_matx[cols, int(rows / 3), 0] = raw_matx[rows, cols]

                elif rows % 3 == 1:

                    freq_matx[cols, int((rows - 1) / 3), 1] = raw_matx[rows, cols]

                elif rows % 3 == 2:

                    freq_matx[cols, int((rows - 2) / 3), 2] = raw_matx[rows, cols]

        freq_modified = np.zeros_like(freq_list)
        for i, each_freq in enumerate(freq_list):

            if each_freq < 0:

                freq_modified[i] = 2

            else:

                freq_modified[i] = each_freq

        return freq_modified, np.array(reduced_mass_list), freq_matx, np.array(force_list), atom_num


class IVR:

    def __init__(self,freq,traj):

        self.freq = freq
        self.traj = traj

    def vibDecomposition(self,t,group="all"):

        if group == "all": group = range(self.freq.atomNum+1)
        velocityTotal = self.traj.VibKE(t,vel=True)
        vibEnergyArr = np.zeros_like(self.freq.freqArr)

        massArr = self.traj.generate_mass_list()
        for counter,Vibmode in enumerate(self.freq.freqArr):

            velofMode = 0
            for atom in range(self.freq.atomNum):

                if atom+1 in group:

                    veldotDisp = np.dot(velocityTotal[atom],self.freq.displaceMatx[counter,atom]) * self.freq.displaceMatx[counter,atom]
                    velofMode += (1.195e-7 * massArr[atom]) * np.sum(veldotDisp** 2)
                    velocityTotal[atom] -= veldotDisp


            vibEnergyArr[counter] = velofMode


        return vibEnergyArr







