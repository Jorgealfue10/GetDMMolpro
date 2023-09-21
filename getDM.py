import torch as th
import numpy as np
import csv
import os

dir="./"

dimensions = {
    '7300.2': 6,
    '7100.2': 5,
    '7400.2': 1,
    '7600.2': 2
}

coordinates = {
    "X": 0,
    "Y": 1,
    "Z": 2
}

nstates = 14

dists = np.loadtxt(dir+"dists.dat",delimiter=None,dtype=float)

def extract_matrix_elements(lines):
    """
    Extracts matrix elements from a list of lines. 

    :param lines: A list of lines containing matrix elements.
    :type lines: list

    :return: A tuple containing the matrix elements, file1, and file2.
    :rtype: tuple
    """
    # Gets the file names containing the wf
    file1 = lines[4].split()[5]
    file2 = lines[5].split()[5]
    # Each filename has a defined number of states
    dim_ket = dimensions[file1]
    dim_bra = dimensions[file2]
    dm = np.zeros((dim_bra, dim_ket, 3))

    # Extracts the matrix elements
    for line in lines[6:]:
        if '!MRCI expec' in line or '!MRCI trans' in line:
            #print(line)
            splited = line.split()
            if "L" in splited[2]:
                continue
            # Extracts index of the matrix
            dims = splited[2]
            dim1, dim2, coor = int(dims[1]) - 1, int(dims[9]) - 1, coordinates[dims[7]]
            
            # Extracts the matrix value
            dm[dim1, dim2, coor] = float(splited[3])

    return dm,file1,file2

dm_tot_comp = np.zeros((100,nstates, nstates,3))
dm_tot = np.zeros((100,nstates, nstates))

for file in os.listdir(dir):
    print(file)
    if file.endswith(".out"):
        with open(dir+file) as f:
            lines = f.readlines()

        k = int(file[9:12])-1
        istart, istop = [], []
        i,j = 0, 0
        check = False

        nstates = 14

        for iline, line in enumerate(lines):
            if "Transition moment calculation" in line:
                istart.append(iline)
                check = True
            if '**********************************************************************************************************************************' in line and check:
                istop.append(iline)
                check = False

        for i in range(len(istart)):
            dm,file2,file1 = extract_matrix_elements(lines[istart[i]:istop[i]])
            #print(dm)
            #print(i,istart[i],istop[i],file1,file2)
            if file1 == "7300.2":
                if file2 == file1:
                    dm_tot_comp[k,:6,:6,:] = dm
                elif file2 == "7100.2":
                    dm_tot_comp[k,:6,6:11,:] = dm
            elif file1 == "7100.2":
                if file2 == file1:
                    dm_tot_comp[k,6:11,6:11,:] = dm
                elif file2 == "7400.2":
                    dm_tot_comp[k,6:11,:6,:] = dm
            elif file1 == "7400.2":
                if file2 == file1:
                    dm_tot_comp[k,11,11,:] = dm
                elif file2 == "7600.2":
                    dm_tot_comp[k,11,12:,:] = dm
            elif file1 == "7600.2":
                if file2 == file1:
                    dm_tot_comp[k,12:,12:,:] = dm
                elif file2 == "7400.2":
                    dm_tot_comp[k,12:,11:11,:] = dm

        for i in range(nstates):
            for j in range(nstates):
                dm_tot[k,i,j] = np.sqrt((dm_tot_comp[k,i,j,0])**2 + (dm_tot_comp[k,i,j,1])**2 + (dm_tot_comp[k,i,j,2])**2)
        

with open("DM-GS-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,0,0,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-2Delta-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,1,1,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-12Sigma+-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,2,2,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-22Sigma+-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,3,3,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-2Pi-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,4,4,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-22Pi-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,5,5,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-GS-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,6,6,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-2Delta-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,7,7,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-12Sigma-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,8,8,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-2Pi-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,9,9,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-22Pi-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,10,10,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-14Pi-Ap.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,11,11,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-aSigma-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,12,12,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])

with open("DM-14Pi-App.dat",'w') as f:
    data = np.column_stack((dists,dm_tot_comp[:,13,13,2]))
    np.savetxt(f,data,delimiter=' ',fmt=['%12.8f','%16.12f'])