__author__ = 'hyojungkim'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import font_manager,rc

def read_strain(s):

    grid = []
    atom=[]
    for line in s.splitlines()[:]: ## skip first 2 lines of xyz file which are comments
        if line != '':
            entries = line.split()
            if len(entries)==5:
                ## atom=[distance from core(2D),x,y,z,max_disp between neighbors, max_angle between neighbors]
                atom.append([float(entries[0]),float(entries[1]),float(entries[2]), float(entries[3]),float(entries[4])])
            else: grid.append([float(entries[0]),float(entries[1]),float(entries[2])])
    print len(atom), 3*len(atom), len(grid)
    return atom, grid

def plot_strain_exx(atom, strain, type, Rcut):
    ### strain matrix
    ### [exx,exy,exz],
    ### [eyx,eyy,eyz],
    ### [ezx,ezy,ezz]
    dist_2d, x,y,z,max_dist=zip(*atom)
    print 'max,dist:' , np.max(max_dist)
  #  dist_2d=[np.linalg.norm([x[i],y[i]]) for i in range(0,len(atom))]

    exx=[strain[3*i][0] for i in range(0,len(atom))]
    eyy=[strain[3*i+1][1] for i in range(0,len(atom))]
    ezz=[strain[3*i+2][2] for i in range(0,len(atom))]
    exy=[strain[3*i][1] for i in range(0,len(atom))]
    eyx=[strain[3*i+1][0] for i in range(0,len(atom))]
    exz=[strain[3*i][2] for i in range(0,len(atom))]
    ezx=[strain[3*i+2][0] for i in range(0,len(atom))]
    eyz=[strain[3*i+1][2] for i in range(0,len(atom))]
    ezy=[strain[3*i+2][1] for i in range(0,len(atom))]

    avg_exy=[(strain[3*i][1]+strain[3*i+1][0])/2.0 for i in range(0,len(atom))]
    avg_exz=[(strain[3*i][2]+strain[3*i+2][0])/2.0 for i in range(0,len(atom))]
    avg_eyz=[(strain[3*i+1][2]+strain[3*i+2][1])/2.0 for i in range(0,len(atom))]

    avg_exx_eyy=[ (exx[i]+eyy[i])/2.0 for i in range(0,len(atom))]
    diff_exx_eyy=[ (exx[i]-eyy[i])/2.0 for i in range(0,len(atom))]
    shear_strain_error=[np.sqrt((exx[i]-eyy[i])**2+(exy[i]+eyx[i])**2)/2.0 for i in range(0,len(atom))]
    if type==1: ## plot exx, eyy, ezz
        f,axarr=plt.subplots(9, sharex=True, sharey=True)
        axarr[0].plot(dist_2d,exx,'b-',label='exx' )
        axarr[1].plot(dist_2d,eyy, 'r-', label='eyy')
        axarr[2].plot(dist_2d,ezz, 'k-', label='ezz')
        axarr[3].plot(dist_2d,exy, 'g-', label='exy')
        axarr[4].plot(dist_2d,eyx, 'g-', label='exz')
        axarr[5].plot(dist_2d,exz, 'm-', label='eyz')
        axarr[6].plot(dist_2d,ezx, 'm-', label='exy')
        axarr[7].plot(dist_2d,eyz, 'c-', label='exz')
        axarr[8].plot(dist_2d,ezy, 'c-', label='eyz')

       # axarr[0].set_legend(loc='upper left')
    if type==2:
        f,axarr=plt.subplots(6, sharex=True,sharey=True)
        axarr[0].plot(dist_2d,avg_exy, 'g-', label='(exy+eyx)/2')
        axarr[1].plot(dist_2d,avg_exz, 'm-', label='(exz+ezx)/2')
        axarr[2].plot(dist_2d,avg_eyz, 'c-', label='(eyz+ezy)/2')
        axarr[3].plot(dist_2d, avg_exx_eyy, 'k-', label='(exx+eyy)/2')
        axarr[4].plot(dist_2d, diff_exx_eyy, 'b-', label='(exx-eyy)/2' )
        axarr[5].plot(dist_2d, shear_strain_error, 'y-', label='shear strain error' )
        axarr[0].legend(bbox_to_anchor=(0.38, 1.13))
        axarr[1].legend(loc='upper left')
        axarr[2].legend(loc='upper left')
        axarr[3].legend(loc='upper left')
        axarr[4].legend(loc='upper left')
        axarr[5].legend(bbox_to_anchor=(0.4, 1.13))
   #     axarr[5].set_ylim([0.000, 0.05])
    plt.xlim([0,Rcut])
#    plt.ylim([-0.0019,0.0021])
#    plt.yticks(np.arange(-0.002,0.0021, step=0.002))
#    plt.ylim([-0.1,0.1])
    plt.ylim([-0.0009, 0.0031])
    plt.yticks(np.arange(-0.001,0.0031, step=0.001))
   # plt.legend(loc='upper left')
    axarr[0].set_title('disl: Rcut= %d' %Rcut,fontsize=20)

    plt.xlabel('distance from core ($\AA$)', fontsize=20)
    f.text(0.02, 0.5, 'strain error', va='center', rotation='vertical', fontsize=20)

#    axarr[0].ylabel('strain error',fontsize=20)
    plt.show()

def get_strain(atom, strain, type,R_clamp,reverse):
    dist_2d, x,y,z,max_disp=zip(*atom)
    if reverse=='True':
        dist_2d=[R_clamp-dist for dist in dist_2d]
#    print len(atom)
    print 'max,disp:' , np.max(max_disp), np.min(max_disp)

    if type=='xx':
        return dist_2d, [abs(strain[3*i][0]) for i in range(0,len(atom))]
    if type=='yy':
        return dist_2d, [abs(strain[3*i+1][1]) for i in range(0,len(atom))]
    if type=='zz':
        return dist_2d,[abs(strain[3*i+2][2]) for i in range(0,len(atom))]
    if type=='xy':
        return dist_2d,[abs(strain[3*i][1]) for i in range(0,len(atom))]
    if type=='volumetric':
        return dist_2d, [abs(strain[3*i][0]+strain[3*i+1][1]+strain[3*i+2][2])/3.0 for i in range(0,len(atom)) ]
    if type=='shear':
        exx=[strain[3*i][0] for i in range(0,len(atom))]
        eyy=[strain[3*i+1][1] for i in range(0,len(atom))]
        ezz=[strain[3*i+2][2] for i in range(0,len(atom))]
        exy=[strain[3*i][1] for i in range(0,len(atom))]
        eyx=[strain[3*i+1][0] for i in range(0,len(atom))]
        exz=[strain[3*i][2] for i in range(0,len(atom))]
        ezx=[strain[3*i+2][0] for i in range(0,len(atom))]
        eyz=[strain[3*i+1][2] for i in range(0,len(atom))]
        ezy=[strain[3*i+2][1] for i in range(0,len(atom))]
        return dist_2d, [np.sqrt((exx[i]-eyy[i])**2+(eyy[i]-ezz[i])**2+(ezz[i]-exx[i])**2+
                                 3*(exy[i]**2+eyz[i]**2+ezx[i]**2+eyx[i]**2+ezy[i]**2+exz[i]**2)) for i in range(0,len(atom))]
 #       print dist,m
        return dist,m

def func(x,y,b):
    x_log=[np.log10(xx) for xx in x if xx>1 and xx< 500]
    y_log=[np.log10(yy) for xx,yy in zip(x,y) if xx>1  and xx< 500]

    a,c= np.polyfit(x_log,y_log,1)
    print a,c
    x_fit=[xx for xx in x if xx>1 and xx<500]
    y_fit=[10**c*x_val**a for x_val in x_fit]

    return x_fit,y_fit

if __name__ =='__main__':

    """READ STRAIN ERRORS OF BULK/DISLOCATION GEOMETRY"""
    path="./TiMo_strain_error/"
    geom='bulk'  ## either bulk or disl
    type='shear'
    boundary=1
    atom_R100, strain_R100 =  read_strain(open(path+'strain_error_TiMo_R100.txt','r').read())
 #   atom_R200, strain_R200 =  read_strain(open(path+'strain_error_list_%s_R200.txt'%geom,'r').read())
    atom_R300, strain_R300 =  read_strain(open(path+'strain_error_TiMo_R300.txt','r').read())
#    atom_R400, strain_R400 =  read_strain(open(path+'strain_error_list_%s_R400.txt'%geom,'r').read())
    atom_R500, strain_R500 =  read_strain(open(path+'strain_error_TiMo_R500.txt','r').read())

    dist_R100,e_R100= get_strain(atom_R100, strain_R100,type,100,reverse='True')
#    dist_R200,e_R200= get_strain(atom_R200, strain_R200,type,200,reverse='True')
    dist_R300,e_R300= get_strain(atom_R300, strain_R300,type,300,reverse='True')
#    dist_R400,e_R400= get_strain(atom_R400, strain_R400,type,400,reverse='True')
    dist_R500,e_R500= get_strain(atom_R500, strain_R500,type,500,reverse='True')

    axlog=plt.subplot(111)

    axlog.plot(dist_R500,e_R500,color='black', label='R_clamp=500$\AA$')
#    axlog.plot(dist_R400,e_R400, color='magenta',label='R_clamp=400$\AA$')
    axlog.plot(dist_R300,e_R300, color='blue',label='R_clamp=300$\AA$')
#    axlog.plot(dist_R200,e_R200,color='orange',label='R_clamp=200$\AA$')
#    axlog.plot(dist_R100,e_R100,color='green',label='R_clamp=100$\AA$')

    axlog.legend(loc='upper right',fontsize=15)

    for R in [500,300,100]:  #,300,100
        if R==500:
            x_fit,y_fit=func(dist_R500,e_R500,boundary)
        if R==300:
            x_fit,y_fit=func(dist_R300,e_R300,boundary)
        if R==100:
            x_fit,y_fit=func(dist_R100,e_R100,boundary)

        axlog.plot(x_fit,y_fit,'r--')

    axlog.set_xscale('log')
    axlog.set_yscale('log')
    axlog.set_xlim((boundary*10**0,10**3))
    axlog.spines['left'].set_visible(False)
    axlog.yaxis.set_ticks_position('right')
    axlog.yaxis.set_visible(False)
    ax2 = axlog.twinx()
    ax2.set_yscale('log')
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(axis='y',which='both',labelright='off')

    divider = make_axes_locatable(axlog)
    axLin = divider.append_axes("left", size=1.2, pad=0, sharey=axlog)
    axLin.plot(dist_R500,e_R500,color='black')
#    axLin.plot(dist_R400,e_R400, color='magenta')
    axLin.plot(dist_R300,e_R300, color='blue')
#    axLin.plot(dist_R200,e_R200,color='orange')
    axLin.plot(dist_R100,e_R100,color='green')

    axLin.set_xscale('linear')
    axLin.set_yscale('log')
    axLin.set_xlim((-8, boundary*10**0))
    axLin.spines['right'].set_visible(False)
    axLin.yaxis.set_ticks_position('left')
    axLin.set_xticks([-8,-6,-4,-2])

#    plt.setp(axLin.get_xticklabels(), visible=True)
#    plt.setp(axLin.get_xticks([-7,-5,-3,-1]),visible=True)
    plt.xlabel('R_clamp-distance ($\AA$)', fontsize=18,horizontalalignment='right', position=(3.4,25))
    plt.ylabel('|strain error|', fontsize=18)
#    plt.legend(loc='upper right',fontsize=18)
    if geom=='bulk':
        plt.title('Bulk, e_%s'%type,fontsize=18,horizontalalignment='right', position=(3.3,1))
    if geom=='disl':
        plt.title('Dislocation, e_%s'%type,fontsize=18,horizontalalignment='right', position=(3.4,1))

 #   plt.xlim([10**-1,10**3])
    plt.ylim([10**-6,10**0])

    plt.show()
