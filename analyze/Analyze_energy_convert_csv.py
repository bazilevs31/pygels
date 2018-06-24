import sys, string
import re
import numpy as np
import pandas as pd
import glob
import os.path


def get_datname(_filename):
    """reads profile file and analyzes which data file it belonged to"""
    searchObj = re.search( r'(.*)_sc(.*?)_press.txt', _filename)
    if searchObj:
        found = searchObj.group(2)
    print found
    dataname = './data/sc'+found+'.dat'
    return dataname

def get_profile_filename(_filename):
    """reads profile file and analyzes which data file it belonged to"""
    searchObj = re.search( r'(.*?)_sc(.*?)_coulsim_(.*?).txt', _filename)
    if searchObj:
        found1 = searchObj.group(1)
        found2 = searchObj.group(2)
        found3 = searchObj.group(3)
    print found1, found2
    profiletype = found1
    searchfilename=re.search(r'(.*?)_ljeq1(.*)',found2)
    dataname = 'sc'+searchfilename.group(1)
    simtype = found3
    return profiletype, dataname, simtype

def get_info_dat(inFileName):
    """analyzes a .dat file in the ./data/....dat
    directory to get number of atoms and dimensions
    """
    inFile = open(inFileName, "r")
    lines = inFile.readlines()
    inFile.close()

    verbose = False

    for line in lines:
        if "atoms" in line:
            list_line = line.split(" ")
            N_atoms = int(list_line[0])
            if verbose:
                print "N_atoms", N_atoms
        elif "xlo xhi" in line:
            list_line = line.split(" ")
            xlo = float(list_line[0])
            xhi = float(list_line[1])
            if verbose:
                print xlo, xhi
        elif "ylo yhi" in line:
            list_line = line.split(" ")
            ylo = float(list_line[0])
            yhi = float(list_line[1])
            print ylo, yhi
        elif "zlo zhi" in line:
            list_line = line.split(" ")
            zlo = float(list_line[0])
            zhi = float(list_line[1])
            if verbose:
                print zlo, zhi
    return N_atoms, abs(xlo - xhi), abs(ylo - yhi), abs(zlo - zhi)


def analyze_energy_profile(inFileName):
    """analyzes pressure profile and saves the output data onto
    csv file
    with press<name of atomtype><name pressure type>_<gel name>.csv
    """

    # inFileName = "energyalltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt"
    inFile = open(inFileName, "r")
    lines = inFile.readlines()
    inFile.close()


    # sc1844_nm100_l50_fl01_fr05_t11_ljeq1_coul_eq_coulsim.dat




    # find slab thickness (delta):
    for line in lines:
        if line[0] != '#': # ignore comments
            words = string.split(line)
            if len(words) == 2 or len(words) == 3:
                nBins = int(words[1])

    # datname = get_datname(inFileName)
    # N_atoms, lx, ly, lz = get_info_dat(datname)
    # slabVolume = lx * ly * lz / (nBins - 1)
    # delta = lz / float(nBins - 1)

    x = np.zeros(nBins, dtype=[('Id', np.float32), ('Coord1', np.float32), ('Ncount', np.float32),  ('NumDensity', np.float32), ('Energy', np.float32)])


    flag_reading = True
    traj_pd = []
    count = 0


    # calc & output P tensor components (file) and P_L-P_N (screen):
    # outFile = open('zP_xx_yy_zz', 'w')
    for line in lines:
        flag_reading = True
        if line[0] != '#': # ignore comments
            words = string.split(line)
            if flag_reading and (len(words) != 3):
                # print "reading Energy data, count", count
                x[count]['Id'] = int(words[0])
                x[count]['Coord1'] = float(words[1])
                x[count]['NumDensity'] = float(words[3])
                nCount = float(words[2])
                x[count]['Ncount'] = nCount
                x[count]['Energy'] = float(words[4]) * nCount
                df = pd.DataFrame.from_records(x)
                # print df
                count += 1
            if len(words) == 3:
                # flag_reading = False
                # print "len is 3"
                if count > 0:
                    traj_pd.append(df)
                    # traj_pandas = []
                    # flag_reading = False
                    frame = int(words[0])
                    count = 0

    profiletype, outputfile, simtype = get_profile_filename(inFileName)
    Combine_PD = pd.concat(traj_pd)
    # Combine_PD.to_csv(profiletype+"_"+outputfile + '.csv', encoding='utf-8', index=False)
    Combine_PD.to_csv(profiletype+"_"+outputfile + "_" + simtype+ '.csv', encoding='utf-8', index=False)
    return None


def main():
    # directory = 'profiles'
    for fname in glob.glob('./energy*.txt'):
        print fname
        # analyze_energy_profile("energyalltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt")
        if os.path.isfile(fname):
            print "file exists"
            analyze_energy_profile(fname)


if __name__ == '__main__':
    main()
