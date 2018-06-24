import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('presentation')
import argparse
import pandas as pd

# parser = argparse.ArgumentParser(description='program that looks at data file created and analyzes its ions distribution ')
# parser.add_argument('-f', '--filename', default='mol')
# args = parser.parse_args()

args = {}
args['filename'] = 'sc1844_nm100_l50_fl01_fr1_test'

filename = (args['filename']).rsplit('.', 1)[0]
input_data = filename + '.data'
# input_data = 'sc1844_nm100_l50_fl01_fr1' + '.data'
# input_dcd = 'sc1844_nm100_l50_fl01_fr1_t03_ljeq1_coul_eq_coulsim_press' + '.dcd'
# u = MDAnalysis.Universe(input_data, input_dcd)
# u.trajectory[1]

u = MDAnalysis.Universe(input_data)

u.universe.atoms.pack_into_box()
all_atoms = u.select_atoms("all")

x_all = u.atoms.positions[:,2]
N_all = len(x_all)

ci_indices = np.where(u.atoms.charges==1.)
gel_all_indices = np.where(u.atoms.charges<=0.)
gel_neut_indices = np.where(u.atoms.charges==0.)
gel_ions_indices = np.where(u.atoms.charges==-1.)


x_countions = x_all[ci_indices]
N_cions = len(x_countions)

x_gel_all = x_all[gel_all_indices]
N_gel_all = len(x_gel_all)
x_gel_neut = x_all[gel_neut_indices]
N_gel_neut = len(x_gel_neut)
x_gel_ions = x_all[gel_ions_indices]
N_gel_ions = len(x_gel_ions)


density_histograms = {}
Nbins = 10
lx = u.trajectory.ts.dimensions[0]
ly = u.trajectory.ts.dimensions[1]
lz = u.trajectory.ts.dimensions[2]
A = lx*ly
dv = A*lz/float(Nbins)

hist_array_ci, bin_all = np.histogram(x_countions, Nbins)
hist_array_g_ions, bin_all = np.histogram(x_gel_ions, Nbins)
hist_array_g_all, bin_all = np.histogram(x_gel_all, Nbins)
hist_array_g_neut, bin_all = np.histogram(x_gel_neut, Nbins)
hist_array_all, bin_all = np.histogram(x_all, Nbins)
bin_plot = 0.5*(bin_all[1:] + bin_all[:-1])


# density_histograms['x/lx'] = bin_plot.flatten()/bin_plot.max()
# density_histograms['counter_ions'] = hist_array_ci.flatten()/float(N_cions)
# density_histograms['gel_neut'] = hist_array_g_neut.flatten()/float(N_gel_neut)
# density_histograms['gel_ions'] = hist_array_g_ions.flatten()/float(N_gel_ions)
# density_histograms['gel_all'] = hist_array_g_all.flatten()/float(N_gel_all)
# density_histograms['all_atoms'] = hist_array_all.flatten()/float(N_all)


# for fname in density_histograms.keys():
#     if fname is 'x/lx':
#         continue
#     print fname
#     plt.plot(density_histograms['x/lx'], density_histograms[fname], label=fname,markersize=20)
# plt.title(filename)
# plt.xlabel('bin in x')
# plt.ylabel('density')
# plt.grid('on')
# plt.legend(loc='best')
# plt.savefig('init_dens_' + filename+'.pdf')



plt.clf()
numbers_histograms = {}
numbers_histograms['x/lx'] = bin_plot.flatten()/bin_plot.max()
numbers_histograms['counter ions'] = hist_array_ci.flatten()
# numbers_histograms['gel_neut'] = hist_array_g_neut.flatten()
# numbers_histograms['gel_ions'] = hist_array_g_ions.flatten()
numbers_histograms['backbone ions'] = hist_array_g_all.flatten()
# numbers_histograms['all_atoms'] = hist_array_all.flatten()


for fname,marker in zip(numbers_histograms.keys(), ['ro--','go:']):
    if fname is 'x/lx':
        continue
    print fname
    plt.plot(numbers_histograms['x/lx'][:-1], numbers_histograms[fname][:-1]/dv, marker, label=fname,lw=2.5)

# plt.grid('on')
plt.legend(loc='best', framealpha=0.0)
# plt.title(filename)
plt.xlabel('$z/l_z$')
plt.ylabel('$c(z/l_z), 1/\sigma^3$')
plt.savefig('init_numb_' + filename+'.pdf')
plt.savefig('init_numb_' + filename+'.png')


# df = pd.DataFrame({"name1" : a, "name2" : b})

