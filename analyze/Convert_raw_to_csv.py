#!/usr/bin/env python

import numpy as np

import pandas as pd
from gridData import Grid
import glob
import os.path

import Analyze_file_info
import matplotlib.pyplot as plt
plt.style.use('presentation')


def find_dx_file(dirname, filename):
    print "looking for matches for the file", filename

    _, gelname, _  = Analyze_file_info.get_profile_filename(filename)
    for fname in glob.glob(dirname+'/*.dx'):
        print "looking for dx files"
        # analyze_Pressure_profile("energyalltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt")
        if os.path.isfile(fname) and (gelname in fname):
            print "file exists", "i found it gelnmae, fname, filename ", gelname, fname, filename
            return fname, gelname

dirname = '/Users/bazilevs/Work/cedar/17-11-29_results'
# dirname = '/Users/bazilevs/Work/cedar/18-02-05'
for fname in glob.glob(dirname+'/pressure_profiles_all/pressalltot*.txt'):
# for fname in glob.glob(dirname+'/profiles/pressalltot*.txt'):
        print fname
        # analyze_Pressure_profile("energyalltot_sc1844_nm100_l50_fl02_fr05_t03_ljeq1_coul_eq_coulsim_press.txt")
        if os.path.isfile(fname):
            print "file exists now analyzing the pressure"
            print fname
            print "looking for dx files"
            dxfile, gelname = find_dx_file(dirname + '/dx_all', os.path.basename(fname))

            df_p = Analyze_file_info.analyze_pressure_profile(fname)
            df_phi = Analyze_file_info.analyze_dx(dxfile)

            df_means = df_p.groupby(df_p['Id']).mean()
            df_stds = df_p.groupby(df_p['Id']).std()
            df_stds.add_suffix('_Error')
            df_p_full = pd.concat([df_means, df_stds], axis=1)

            print df_p.head()
            print df_phi.head()
            # result = pd.concat([df_p, df_means, df_stds], axis=1)


            _ = fname.replace('alltot','cipair')
            ci_name = _.replace('total','coul')

            df_ci = Analyze_file_info.analyze_pressure_profile(ci_name)
            df_ci_means = df_ci.groupby(df_ci['Id']).mean()
            df_ci_stds = df_ci.groupby(df_ci['Id']).std()

            result = pd.concat([df_phi, df_means,df_stds.add_suffix('_error'),df_ci_means.add_suffix('_ci'),df_ci_stds.add_suffix('_ci_error')], axis=1)

            result.to_csv('./'+gelname+'.csv')
