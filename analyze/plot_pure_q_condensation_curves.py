import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
plt.style.use('presentation-new')
# plt.style.use('rsc-style')
# plt.style.use('paper-doublefig')
# plt.style.use('seaborn')
import argparse
import Analyze_file_info
import termcolor
import mpltex


def set_size(fig):
    fig.set_size_inches(6, 3)
    plt.tight_layout()

def get_xi_manning(temp, f):
    """
    analyzes the xi from temp, f
    """
    lb = 1./temp
    nu = 0.5
    b = 1.03
    xi = lb/(b*f**(-nu))
    return xi

# @mpltex.rsc_decorator
def main(**params_dict):
    fig = plt.figure(1)
    fig.clf()
    # fig.set_size_inches(5, 5)
    fig.set_size_inches(6.24,7)
    axleft = fig.add_subplot(311)

    axleft.set(xlim=[0, 1.2],ylim=[-0.1,.6], ylabel=r'$f_{\rm free, 1}$')
    # axleft.set(xlim=[0, 1.2],ylim=[-0.1,.6], ylabel=r'$f_{ free, 1}$')

# $f_{\rm free,i}$
    axright = fig.add_subplot(312)
    # axright.set(xlim=[0, 1.2],ylim=[-0.1,1.3], xlabel='$T \ [\\varepsilon]$', ylabel=r'$f_{\rm free, 2}$')
    axright.set(xlim=[0, 1.2],ylim=[-0.1,1.3], ylabel=r'$f_{\rm free, 2}$')

    # fratio = plt.figure(3)
    # fratio.clf()
    # axratio = fratio.add_subplot(111)
    axratio = fig.add_subplot(313)
    axratio.set(xlim=[0, 1.2], xlabel=r'$T \  [\varepsilon]$', ylabel=r'$\dfrac{f_{\rm free, 2}}{f_{\rm free, 1}}$')


    files_list = ['sc1844_nm100_l50_fl01_fr02','sc1844_nm100_l50_fl01_fr03','sc1844_nm100_l50_fl01_fr05','sc1844_nm100_l50_fl01_fr1','sc1844_nm100_l50_fl025_fr05']
    file_color_list = ['r','m','g','b','c']

        # files_list = ['sc1844_nm100_l50_fl01_fr02']
        # file_color_list = ['r','m','g','b','c']
    files_color_dict = dict(zip(files_list, file_color_list))


    for myfile in params_dict['myfiles_list']:
        c = files_color_dict[myfile]
    # for c, myfile in zip(file_mark_list,files):
        # csvfile = 'COND_q_'+myfile+'_v2'+'.csv'
        # csvfile = 'COND_v1_'+filename+'_interpolated'+'.csv'
        params_dict['myfile'] = myfile

        if params_dict['cond_suffix'] is '':
            csvfile = '{0[cond_prefix]}_{0[myfile]}.csv'.format(params_dict)
            csvfile_2 = '{0[cond_prefix]}_{0[myfile]}_new.csv'.format(params_dict)
        else:
            csvfile = '{0[cond_prefix]}_{0[myfile]}_{0[cond_suffix]}.csv'.format(params_dict)
        print(csvfile)

        if os.path.exists(csvfile):
            print("csvfile exists")
            df4 = pd.read_csv(csvfile)
            df2 = pd.read_csv(csvfile_2)
            df4[df4.isnull()] = df2
            df4 = df4.interpolate()

            if 'fl01_fr02' in myfile:
                loc_ind_array = np.arange(len(df4['T']))[1:]
            else:
                loc_ind_array = np.arange(len(df4['T']))


            fl, fr, _, _ = Analyze_file_info.get_profile_params(myfile+'_t01')
            xi_left_array = get_xi_manning(df4['T'][loc_ind_array], fl)
            ind = np.where(xi_left_array>=1.)
            indices = ind[0]
            xi_left_array_full = fl*np.ones_like(xi_left_array)
            xi_left_array = np.array(xi_left_array)
            xi_left_array_full[indices] = fl/xi_left_array[indices]

            xi_right_array = get_xi_manning(df4['T'][loc_ind_array], fr)
            ind = np.where(xi_right_array>=1.)
            indices = ind[0]
            xi_right_array_full = fr*np.ones_like(xi_right_array)
            xi_right_array = np.array(xi_right_array)
            xi_right_array_full[indices] = fr/xi_right_array[indices]
            # xi_right_array = get_xi_manning(df4['T'][loc_ind_array], fr)

            plt.rcParams["font.family"] = "serif"
            # plt.rcParams["mathtext.fontset"] = "dejavuserif"
            plt.rcParams["mathtext.fontset"] = "stix"
            # plt.rcParams["font.serif"] = "Times New Roman"
            plt.rcParams["font.serif"] = "STIXGeneral"

            axleft.plot(df4['T'][loc_ind_array],(1-df4['q1'][loc_ind_array])*fl,c+'o--',fillstyle='none',mec=c,lw=.8,alpha=0.6,label='$f_l$={fl}'.format(fl=fl))
            axleft.plot(df4['T'][loc_ind_array],xi_left_array_full,c+'o-',fillstyle='full',markersize=5.,lw=.5,label='$f_l$={fl} Manning'.format(fl=fl))
            axright.plot(df4['T'][loc_ind_array],(1-df4['q2'][loc_ind_array])*fr,c+'o--',fillstyle='none',mec=c,lw=.8,alpha=0.6,label='$f_r$={fr}'.format(fr=fr))
            axright.plot(df4['T'][loc_ind_array],xi_right_array_full,c+'o-',fillstyle='full',markersize=5.,lw=.5,label='$f_l$={fr} Manning'.format(fr=fr))
            # axright.plot(df4['T'][loc_ind_array],xi_right_array,c+'o--',fillstyle='full',lw=.8,label='$f_l$={fl} Manning'.format(fl=fl))
            # axratio.plot(df4['T'][loc_ind_array],(1-xi_left_array)/(1-xi_right_array),c+'o--',fillstyle='full',lw=.8,label='$f_l$={fl} Manning'.format(fl=fl))
            axratio.plot(df4['T'][loc_ind_array],(xi_left_array_full)/(xi_right_array_full),c+'o-',fillstyle='full',lw=.8,markersize=5.,label='$f_l$={fl} Manning'.format(fl=fl))
            # axright.plot(df4['T'][loc_ind_array],df4['q2'][loc_ind_array],c+'o--',fillstyle='none',mec=c,lw=.8,label='$f_r$={fr}'.format(fr=fr))
            axratio.plot(df4['T'][loc_ind_array],(fl/fr)*(1-df4['q1'][loc_ind_array])/(1-df4['q2'][loc_ind_array]),c+'o--',fillstyle='none',alpha=.6,mec=c,lw=.8,label='$f_l$={fl}, $f_r$={fr}'.format(fl=fl,fr=fr))

            df4['param'] = (1-df4['q1'])/(1-df4['q2'])
            df4.to_csv('{0[cond_prefix]}_{0[myfile]}_merged.csv'.format(params_dict))


    # line, = ax.plot([1, 2, 3], label='Inline label')

    # box = axright.get_position()
    # axright.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # # Put a legend to the right of the current axrightis
    # axright.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # axleft.legend(frameon=False,loc='best', ncol=2)
    # # # fleft.show()

    # axright.legend(frameon=False,loc='best', ncol=2)
    # axratio.legend(frameon=False, loc='best', ncol=2)
    # fright.show()
    keyword = 'q2_mann'
    # fleft.savefig('left_cond_{keyword}.png'.format(keyword=keyword))
    # fright.savefig('right_cond_{keyword}.png'.format(keyword=keyword))
    # set_size(fig)
    # set_size(fratio)
    # fratio.savefig('ratio_cond_{keyword}.png'.format(keyword=keyword))


    width_single_column = 3.26  # 8.3 cm
    width_double_column = 6.73  # 17.1 cm
    GOLDEN_RATIO = 1/1.61803398875
    # # Default ratio for a single plot figure
    # # I prefer a little higher than goden ratio, from 0.618 to about 0.68
    height_width_ratio = GOLDEN_RATIO * 1.1  # = height / width

    _width = width_single_column
    _height = width_single_column * height_width_ratio
    # plt.figsize: (_width _height))
    plt.rcParams["figure.figsize"] = [_width,_height]


    # draw vertical lines for all of the axes
    xcoords = [0.5, 0.7]
    for a in fig.axes:
        for xc in xcoords:
            a.axvline(x=xc,color='c', linestyle='--',lw=0.3)
    a.text(.2, 1., 'Regime 3', fontsize=15, color='c')
    a.text(.5, 1., 'Regime 2', fontsize=15, color='c')
    a.text(.8, 1., 'Regime 1', fontsize=15, color='c')
    # plt.axvline([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    # fratio.savefig('./ratio_cond_{keyword}.pdf'.format(keyword=keyword))
    # fratio.savefig('./ratio_cond_{keyword}.pdf'.format(keyword=keyword))
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    # plt.tight_layout()
    # fig.savefig('all_{keyword}.png'.format(keyword=keyword))
    fig.savefig('./all_{keyword}.pdf'.format(keyword=keyword))
    return None
# keyword = '_v2_pure'
# for ext in ['.png', '.pdf']:
#     ff.savefig('phi_lnf1f2' + keyword + ext)
#     fp.savefig('phi_deltap' + keyword + ext)
#     fc.savefig('phi_lnc1c2' + keyword + ext)
#     fpT.savefig('press_T' + keyword + ext)
#     fcond.savefig('Condensation_index' + keyword + ext)

# fleft.legend()
# fleft



if __name__ == '__main__':
    # main()


    def is_valid_file(parser, arg):
        """
        Check if arg is a valid file

        Args:
            parser (TYPE): Description
            arg (TYPE): Description

        Returns:
            TYPE: Description
        """
        arg = os.path.abspath(arg)
        if not os.path.exists(arg):
            parser.error("The file %s doesn't exist " % arg)
        else:
            return arg


    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # group = parser.add_mutually_exclusive_group()
    # group.add_argument('--lines', action='store_true', help='Plot the data with the lines style')

    parser.add_argument('-a', action='append', dest='myfiles_list',
                        type=str,
                        default=[],
                        help='Add values to analyze values to a list type (default: %(default)s)',
                        )

    parser.add_argument("-cs",
                        "--cond_suffix",
                        dest="cond_suffix",
                        default='interpolated',
                        type=str,
                        help="cond_suffix for the output csv file with the q1,q2 information (default: %(default)s)")


    parser.add_argument("-cp",
                        "--cond_prefix",
                        dest="cond_prefix",
                        default='COND_q1_method',
                        type=str,
                        help="prefix for the output csv file with the q1,q2 information (default: %(default)s)")

    args = parser.parse_args()
    print(args)
    print(termcolor.colored('parameters have been red', 'green'))


    main(**vars(args))
