import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import Analyze_file_info
import pandas as pd
from cycler import cycler
def plot_axes(x, y, ax, **kwargs):
    """plot two arrays with a certain style in the current axes"""
    return None

def get_temp_arrays_from_file(**kwargs):
    """plots arrays from a parameter file"""
    return None

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.

    Args:
        n (TYPE): Description
        name (str, optional): Description

    Returns:
        TYPE: Description
    '''
    return plt.cm.get_cmap(name, n)



# for myfile in fileslist:
    # get arrays using get_temp_arrays_from_file
    # finish everything

    # folder = '/Users/bazilevs/Work/cedar/17-12-08/rerun_profiles/'
folder = '/Users/bazilevs/Work/cedar/18-06-12_energy_density_curves/'

# entypelist = [ 'gelpair', 'allpair', 'cipair']
entypelist = ['cipair']
# for simtype in ['coul', 'total', 'lj']:
# for simtype in ['coul', 'lj']:
#     for entype in entypelist:
#         print( simtype)
#         plot_energy_for_type(folder, myenergytype=entype, simtype=simtype)
simtype = 'cipair'
temp_name_array = ['01','02','03','04','05','06','07','08','09','10','11']
temp_color_array = get_cmap(len(temp_name_array))
gelname = 'sc1844_nm100_l50_fl01_fr1'
# params_dict = vars(args)
myenergytype = 'coul'
params = {
'folder' : folder,
'energytype' : myenergytype,
'gelname' : gelname,
'simtype' : simtype,
'scale_power' : 3,
}
scale_factor =  10.**params['scale_power']

# get_cmap
# frho = plt.figure(1)
# fz = plt.figure(2)
# # fz.clf()
# frho.clf()
# axrho = frho.add_subplot(111)

# axz1 = fz.add_subplot(211)
# axz2 = fz.add_subplot(212)

# if simtype is 'coul':
#     axz = axz1
# else:
#     axz = axz2
cc = (cycler(color=list('rgb')) *
        cycler(linestyle=['-', '--', '-.']))

for _tname in temp_name_array:
    params['tempname'] = _tname
    mystring = "{0[folder]}/{0[gelname]}_t{0[tempname]}.csv".format(params)
    df = pd.read_csv(mystring)
    plt.plot(df['Coord1'], df['NumDensity_ci'])

plt.show()
    # params_dict['temp_name_array'] = temp_name_array
    # params_dict['ax_dict'] = ax_dict
    # params_dict['fig_dict'] = fig_dict
    # for myfile in params_dict['myfiles_list']:
    #     c = files_color_dict[myfile]
    #     params_dict['myfile'] = myfile
    #     temp_array, phi3_array_pure, conc_divide_array, f_divide_array, q1_q2_param_array, press_array,work_array, pzz_array, fillstyle = get_arrays_from_file(**params_dict)
