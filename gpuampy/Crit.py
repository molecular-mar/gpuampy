#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tools for working with GPUAM Crit Calculations.
"""

class Crit:
    """
    Class for recovering and storing data from a GPUAM Crit calculation.

    """

    def __init__(self, filename):
        self.filename = filename
        self.dataframes = [None,None,None,None]
        self.last_common_line = 51
        self.line_to_skip = 0
        self.n_cp_lines = [26,34,26,26]
        self.n_cp = [0,0,0,0]
        self.nuclei = None
        self.nuclei_columns = ['id','type','coord','density','laplacian']
        self.cp_columns = [
                      ['id','coord','density','normgrad','laplacian',
                       'rho/lap','G','K','virial','lagran_rho',
                       'G/rho','elf','lol','eigenvalues','shannon_rho',
                       'diseq_rho','complex_rho','shannon_sigma',
                       'diseq_sigma','complex_sigma'],
                      ['id','atom1_id','atom1_type','atom2_id',
                       'atom2_type','d_attract','d_att1','d_att2',
                       'angle','coord','density','normgrad','laplacian',
                       'rho/lap','G','K','virial','lagran_rho',
                       'G/rho','elf','lol','eigenvalues','ellip',
                       'eta','shannon_rho','diseq_rho','complex_rho',
                       'shannon_sigma','diseq_sigma','complex_sigma']
                     ]

    def get_basic_data(self,data_lines):
        self.n_primitives = int(data_lines[34].split(':')[1])
        self.n_orbitals = int(data_lines[35].split(':')[1])
        self.n_nuclei = int(data_lines[36].split(':')[1])
        self.n_electrons = int(data_lines[37].split(':')[1].split()[0])
        self.n_cp[0] = int(data_lines[38].split(':')[1])
        self.n_cp[1] = int(data_lines[39].split(':')[1])
        self.n_cp[2] = int(data_lines[40].split(':')[1])
        self.n_cp[3] = int(data_lines[41].split(':')[1])
        self.n_totalcp = int(data_lines[42].split(':')[1])
        self.poincare_hopf = int(data_lines[45].split('=')[1])

    def set_initial_data(self):
        self.n_dif_cp = ((self.n_cp[0] != 0) + (self.n_cp[1] != 0) 
          + (self.n_cp[2] != 0) + (self.n_cp[3] != 0))
        # 2 lines for the R,B,C CP header. 1 separator
        # line for each kind of CP (n_dif_cp). 2 more
        # lines from the first CP header.
        self.line_to_skip = ( self.last_common_line + self.n_nuclei 
                    + self.n_totalcp + 2 + self.n_dif_cp + 2 )

    def set_nuclei_data(self,data_lines):
        n_lines = self.n_nuclei
        #data_lines[line_to_skip + 5 * n_bcp_lines]
        nuclei_data = [[] for i in range(self.n_nuclei)]
        for i_nuclei in range(n_lines):
            line_data = data_lines[51+i_nuclei].split()
            nuclei_data[i_nuclei].append(int(line_data[0]))
            nuclei_data[i_nuclei].append(line_data[1])
            nuclei_data[i_nuclei].append(list(map(float,line_data[2:5])))
            nuclei_data[i_nuclei].append(float(line_data[5]))
            nuclei_data[i_nuclei].append(float(line_data[6]))
        self.nuclei = pd.DataFrame(nuclei_data,columns=self.nuclei_columns)


    def read_block(self,id,data_lines):
        n_lines = self.n_cp_lines[id]
        #data_lines[line_to_skip + 5 * n_bcp_lines]
        cp_data = [[] for i in range(self.n_cp[id])]
        for i_cp in range(self.n_cp[id]):
            for n_line in range(n_lines):
                if n_line == 0: continue
                current_line = self.line_to_skip + n_line + (i_cp * n_lines)
                line_data = data_lines[current_line].split()
                if not line_data: continue
                elif line_data[0] == 'Critical':
                    cp_data[i_cp].append(int(line_data[2]))
                elif line_data[0] == 'Between':
                    cp_data[i_cp].append(int(line_data[3]))
                    cp_data[i_cp].append(line_data[4])
                    cp_data[i_cp].append(int(line_data[6]))
                    cp_data[i_cp].append(line_data[7])
                elif line_data[0] == 'Coordinates':
                    cp_data[i_cp].append(list(map(float,line_data[-3:])))
                elif line_data[0] == 'Eigenvalues':
                    cp_data[i_cp].append(list(map(float,line_data[-3:])))
                else:
                    cp_data[i_cp].append(float(line_data[-1]))
        
        if id == 1:
            cp_columns = self.cp_columns[1]
        else:
            cp_columns = self.cp_columns[0]
        
        self.dataframes[id] = pd.DataFrame(cp_data,columns=cp_columns)
        self.line_to_skip += self.n_cp[id] * (self.n_cp_lines[id]) + 2
         

    def read_data(self):
        try:
            inp = open(self.filename,'r')
            data = inp.read()
            inp.close()
            data_lines = data.strip().split('\n')
            pass
        except FileNotFoundError:
            print(f"File not found: {self.filename}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")

        self.get_basic_data(data_lines)
        self.set_initial_data()
        self.set_nuclei_data(data_lines)

        for cp_id in range(4):
            if self.n_cp[cp_id] == 0: continue
            else: self.read_block(cp_id,data_lines)

    def get_nuclei(self):
        return self.nuclei
        
    def get_data(self):
        return self.dataframes

    def get_nna(self):
        return self.dataframes[0]

    def get_bcp(self):
        return self.dataframes[1]

    def get_rcp(self):
        return self.dataframes[2]

    def get_ccp(self):
        return self.dataframes[3]

# Example usage for displaying critical points data:
# import gpuampy
# gpuam_data = gpuampy.Crit("moleculeCrit.log")
# gpuam_data.read_data()
# gpuam_data.get_bcp()

