import numpy as np
import os


nstates = 34

dir="./"

dists = np.loadtxt(dir+"dists.dat",delimiter=None,dtype=float)

def get_DM(lines,nstates):
    """
    Given a list of lines and the number of states,
    calculates and returns a density matrix (dm) as a numpy array.
    
    Args:
    - lines: list of strings representing the input data.
    - nstates: integer representing the number of states.
    
    Returns:
    - dm: numpy array representing the density matrix.
    """
    n_el = 8
    num_lines = int(nstates/8)
    last_line = int(nstates%8)
    dm = np.zeros((nstates,nstates))

    n=0
    for line in lines:
        splited = line.split()
        #print(splited)
        if "VALUE:" in splited and n != num_lines:
            for i in range(n_el):
                nel = n*8+i 
                dm[nel,nel] = float(splited[i+1])
                #print(dm[nel,nel],n,nel)
            n += 1
        elif "VALUE:" in splited and n == num_lines:
            for i in range(last_line):
                nel = n*8+i 
                dm[nel,nel] = float(splited[i+1])
                #print(dm[nel,nel],n,nel)

    return dm

dm_tot = np.zeros((100,nstates, nstates))

for file in os.listdir(dir):
    #print(file)
    if file.startswith("mrci") and file.endswith(".out"):
        #print(file)
        with open(dir+file) as f:
            lines = f.readlines()

        k = int(file[9:12])-1

        for iline, line in enumerate(lines):
            if "Expectation values <i|DMZ|i>" in line:
                #istart.append(iline)
                istart = iline
                check = True
            if 'Transition matrix elements <i|DMZ|1>' in line and check:
                #istop.append(iline)
                istop = iline
                check = False

        dm_tot[k,:,:] = get_DM(lines[istart:istop],nstates)

dir = "./DM_bystate/"


for i in range(nstates):
    for j in range(nstates):
#        dm_notzeros = np.argwhere(dm_tot[:,i,j] != 0)
        #dm_notzeros = np.ma.masked_equal(dm_tot[:,i,j],0)
        if i == j:
            dm_not = np.ma.masked_equal(dm_tot[:,i,j],0)
            dm_notzeros = dm_not.compressed()
        if i != j:
            # dm_not = np.ma.masked_equal(dm_tot[:,i,j],0)
            # dm_notzeros = dm_not.compressed()
            # dm_notzeros = np.zeros((len(dists)))
            dm_notzeros = dm_tot[:,i,j]
        if i+1 < 10:
            if j+1 < 10:
                np.savetxt(dir+"dm_tot_0"+str(i+1)+"_0"+str(j+1)+".dat",dm_notzeros[:],fmt='%8f')
            elif j+1 >= 10:
                np.savetxt(dir+"dm_tot_0"+str(i+1)+"_"+str(j+1)+".dat",dm_notzeros[:],fmt='%8f')
        else:
            if j+1 < 10:
                np.savetxt(dir+"dm_tot_"+str(i+1)+"_0"+str(j+1)+".dat",dm_notzeros[:],fmt='%8f')
            elif j+1 >= 10:
                np.savetxt(dir+"dm_tot_"+str(i+1)+"_"+str(j+1)+".dat",dm_notzeros[:],fmt='%8f')