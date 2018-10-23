__author__ = 'hyojungkim'

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib import font_manager,rc
from scipy.signal import hilbert, chirp

def read_strain(s):

    grid = []
    atom=[]
    for line in s.splitlines()[:]: ## skip first 2 lines of xyz file which are comments
        if line != '':
            entries = line.split()
            if len(entries)==6:
                ## atom=[distance from core(2D),x,y,z,max_disp between neighbors, max_angle between neighbors]
                atom.append([float(entries[0]),float(entries[1]),float(entries[2]), float(entries[3]),float(entries[4]),float(entries[5])])
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

def plot_strain_loglog(atom, strain, type, Rcut):
    ### strain matrix
    ### [exx,exy,exz],
    ### [eyx,eyy,eyz],
    ### [ezx,ezy,ezz]
    dist_2d, x,y,z,max_dist,max_angle=zip(*atom)
#    dist_2d=[np.linalg.norm([x[i],y[i]]) for i in range(0,len(atom))]

    exx=[abs(strain[3*i][0]) for i in range(0,len(atom))]
    eyy=[abs(strain[3*i+1][1]) for i in range(0,len(atom))]
    ezz=[abs(strain[3*i+2][2]) for i in range(0,len(atom))]
    exy=[strain[3*i][1] for i in range(0,len(atom))]
    eyx=[strain[3*i+1][0] for i in range(0,len(atom))]
    exz=[strain[3*i][2] for i in range(0,len(atom))]
    ezx=[strain[3*i+2][0] for i in range(0,len(atom))]
    eyz=[strain[3*i+1][2] for i in range(0,len(atom))]
    ezy=[strain[3*i+2][1] for i in range(0,len(atom))]

    avg_exy=[(strain[3*i][1]+strain[3*i+1][0])/2.0 for i in range(0,len(atom))]
    avg_exz=[(strain[3*i][2]+strain[3*i+2][0])/2.0 for i in range(0,len(atom))]
    avg_eyz=[(strain[3*i+1][2]+strain[3*i+2][1])/2.0 for i in range(0,len(atom))]

    avg_exx_eyy=[(exx[i]+eyy[i])/2.0 for i in range(0,len(atom))]
    diff_exx_eyy=[(exx[i]-eyy[i])/2.0 for i in range(0,len(atom))]
    shear_strain_error=[np.sqrt((exx[i]-eyy[i])**2+(exy[i]+eyx[i])**2)/2.0 for i in range(0,len(atom))]
    if type==1: ## plot exx, eyy, ezz
        f,axarr=plt.subplots(1, sharex=True, sharey=True)
      #  axarr.plot(dist_2d,exx,'b-',label='exx')
        axarr.loglog(dist_2d,exx, basex = np.e, basey = np.e)

      #  axarr.loglog(dist_2d,ezz)
 #   plt.xscale("log")
    plt.plot([np.e**0.1,np.e**6],[np.e**-9.2,np.e**-7.5],'r--')
    x_range=[np.e,np.e**2,np.e**3,np.e**4,np.e**5,np.e**6]
    x_ticks=[1,2,3,4,5,6]
    plt.xticks(x_range,x_ticks)
#    plt.xticks([10,50,100,150,200,250,300],[10,50,100,150,200,250,300])
    plt.yticks([np.e**-2,np.e**-4,np.e**-6,np.e**-8,np.e**-10,np.e**-12,np.e**-14],[-2,-4,-6,-8,-10,-12,-14])

    plt.xlabel('ln(distance)', fontsize=15)
    plt.ylabel('ln(|strain error|)', fontsize=15)
    plt.show()

def get_strain(atom, strain, type,R_clamp,reverse):
    dist_2d, x,y,z,max_disp,max_angle=zip(*atom)
    if reverse=='True':
        dist_2d=[R_clamp-dist for dist in dist_2d]
    print len(atom)
    print 'max,disp:' , np.max(max_disp), np.min(max_disp)
    print 'max,angle:' , np.max(max_angle), np.min(max_angle)

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
    if type=='analysis':
        m=[[i, strain[3*i][0], max_dist[i], x[i],y[i],z[i]] for i in range(0,len(atom)) if (dist_2d[i] <np.e**3.5) and (abs(strain[3*i+1][1]) > np.e**-7)]
        dist=[dist_2d[i] for i in range(0,len(atom)) if (dist_2d[i] <np.e**3.5) and (abs(strain[3*i+1][1]) > np.e**-7)]

        print dist,m
        return dist,m
if __name__ =='__main__':

    """SET UP"""
#    Rcut=300
#    atom, strain = read_strain(open('e_grid_%d.txt'%Rcut,'r').read())
 #   atom, strain =read_strain(open('disl_strain_error_300_sorted.txt','r').read())
#    plot_strain_exx(atom, strain, 1,Rcut)
 #   plot_strain_loglog(atom, strain, 1,Rcut)
    """BULK"""
#    atom_R100, strain_R100 =  read_strain(open('bulk_strain_error_100_sorted.txt','r').read())
#    atom_R300, strain_R300 =  read_strain(open('bulk_strain_error_300_sorted.txt','r').read())
#    atom_R500, strain_R500 =  read_strain(open('bulk_strain_error_500_sorted.txt','r').read())
    """DISLOCATION"""
    atom_R100, strain_R100 =  read_strain(open('disl_strain_error_100_sorted.txt','r').read())
    atom_R300, strain_R300 =  read_strain(open('disl_strain_error_300_sorted.txt','r').read())
    atom_R400, strain_R400 =  read_strain(open('disl_strain_error_400_sorted.txt','r').read())
    atom_R500, strain_R500 =  read_strain(open('disl_strain_error_500_sorted.txt','r').read())
    print strain_R500[-1]
    type='volumetric'
 #   threshold=2.54
    dist_R100,e_R100= get_strain(atom_R100, strain_R100,type,100,reverse='True')
    dist_R300,e_R300= get_strain(atom_R300, strain_R300,type,300,reverse='True')
    dist_R400,e_R400= get_strain(atom_R400, strain_R400,type,400,reverse='True')
    dist_R500,e_R500= get_strain(atom_R500, strain_R500,type,500,reverse='True')

    f,axarr=plt.subplots(1, sharex=True, sharey=True)
#    axarr.loglog(dist_R100,e_R100, basex = np.e, basey = np.e, color='green',label='Rcut=100$\AA$')
#    axarr.loglog(dist_R300,e_R300, basex = np.e, basey = np.e, color='blue',label='Rcut=300$\AA$')
#    axarr.loglog(dist_R400,e_R400, basex = np.e, basey = np.e, color='magenta',label='Rcut=400$\AA$')
#    axarr.loglog(dist_R500,e_R500, basex = np.e, basey = np.e, color='black',label='Rcut=500$\AA$')


    axarr.loglog(dist_R100,e_R100,basex =10, basey = 10, color='green',label='Rcut=100$\AA$')
    axarr.loglog(dist_R300,e_R300,basex =10, basey = 10, color='blue', label='Rcut=300$\AA$')
#    axarr.loglog(dist_R400,e_R400,basex =10, basey = 10, color='yellow', label='Rcut=400$\AA$')
    axarr.loglog(dist_R500,e_R500,basex =10, basey = 10,color='black',label='Rcut=500$\AA$')

    """(Rcut-Distance)"""
    if type=='shear':
        plt.plot([np.e**2, np.e**4.5],[np.e**-2.2, np.e**-2.2],'r--')
        plt.plot([np.e**2, np.e**5.8],[np.e**-4.1, np.e**-7.83],'r--')
        plt.plot([np.e**2, np.e**6.3],[np.e**-4.1, np.e**-8.83],'r--')
        print 'slope=0',(-7.83+4.1)/(5.8-2), (-8.83+4.1)/(6.3-2)
    if type=='xx':
        plt.plot([np.e, np.e**4.5],[np.e**-3.4, np.e**-3.4],'r--')
        plt.plot([np.e, np.e**5.8],[np.e**-3.8, np.e**-8.6],'r--')
        plt.plot([np.e, np.e**6.3],[np.e**-4.1, np.e**-9.5],'r--')
    if type=='xy':
        plt.plot([np.e, np.e**4.5],[np.e**-3.4, np.e**-3.4],'r--')
        plt.plot([np.e, np.e**5.8],[np.e**-3.8, np.e**-8.8],'r--')
        plt.plot([np.e, np.e**6.3],[np.e**-4.1, np.e**-9.7],'r--')


    """DISLOCATION
    if type=='xx':
        plt.plot([np.e**0.1, np.e**4],[np.e**-3.80, np.e**-3.30],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-9.0,np.e**-8.3],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-9.7,np.e**-9.05],'r--')
    if type=='yy':
        plt.plot([np.e**0.1, np.e**4],[np.e**-3.80, np.e**-3.30],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-9.5,np.e**-8.5],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-12,np.e**-10.0],'r--')
    if type=='zz':
        plt.plot([np.e**0.1, np.e**4],[np.e**-4.2, np.e**-3.30],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-11.8,np.e**-9.9],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-12.5,np.e**-10.6],'r--')
    if type=='xy':
        plt.plot([np.e**0.1, np.e**4],[np.e**-4.1, np.e**-3.25],'r--')
        plt.plot([np.e**0.1,np.e**5.5],[np.e**-9.2,np.e**-8.4],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-11.8,np.e**-9.5],'r--')
    if type=='volumetric':
        plt.plot([np.e**0.1, np.e**4],[np.e**-4.1, np.e**-3.60],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-11.3,np.e**-9.7],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-11.3,np.e**-10],'r--')
    if type=='shear':
        plt.plot([np.e**0.1, np.e**4],[np.e**-2.5, np.e**-2.10],'r--')
        plt.plot([np.e**0.1,np.e**5.5],[np.e**-8.4,np.e**-7.2],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-9.4,np.e**-8.3],'r--')

    """
    """BULK
    if type=='xx':
        plt.plot([np.e**0.1, np.e**4],[np.e**-9.05, np.e**-7.8],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-9.7,np.e**-9.05],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-9.7,np.e**-9.15],'r--')
    if type=='yy':
        plt.plot([np.e**0.1, np.e**4],[np.e**-9.05, np.e**-7.8],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-10.3,np.e**-9.2],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-10.4,np.e**-9.63],'r--')
    if type=='zz':
        plt.plot([np.e**0.1, np.e**4],[np.e**-11.5, np.e**-9.3],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-12,np.e**-10.3],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-12.3,np.e**-10.9],'r--')
    if type=='xy':
        plt.plot([np.e**0.1, np.e**4],[np.e**-8.1, np.e**-7.4],'r--')
        plt.plot([np.e**0.1,np.e**5.5],[np.e**-10.6,np.e**-9.0],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-10.6,np.e**-9.25],'r--')
    if type=='volumetric':
        plt.plot([np.e**0.1, np.e**4],[np.e**-11.3, np.e**-9.6],'r--')
        plt.plot([np.e**0.1,np.e**5],[np.e**-11.5,np.e**-10.3],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-11.5,np.e**-10.7],'r--')
    if type=='shear':
        plt.plot([np.e**0.1, np.e**4],[np.e**-6.4, np.e**-6.0],'r--')
        plt.plot([np.e**0.1,np.e**5.5],[np.e**-8.1,np.e**-7.6],'r--')
        plt.plot([np.e**0.1,np.e**6],[np.e**-9,np.e**-8.4],'r--')
    """
    x_range=[np.e,np.e**2,np.e**3,np.e**4,np.e**5,np.e**6]
    x_ticks=[1,2,3,4,5,6]
#    plt.xticks(x_range,x_ticks)
#    plt.xticks([100,300,500],[100,300,500])
#    plt.yticks([np.e**0, np.e**-2,np.e**-4,np.e**-6,np.e**-8,np.e**-10,np.e**-12,np.e**-14],[0,-2,-4,-6,-8,-10,-12,-14])

    plt.xlabel('Rcut-distance', fontsize=18)
#    plt.xlabel('Rcut-distance', fontsize=18)
    plt.ylabel('|strain error|', fontsize=18)
    plt.legend(loc='lower left',fontsize=18)
    plt.title('Dislocation, e_%s'%type,fontsize=18)  #, threshold=2.54$\AA$

    plt.xlim([10**-3, np.e**6.5])
    plt.ylim([10**-5,10**0])
#    plt.xlim([np.e**2, np.e**6.5])
#    plt.ylim([np.e**-14,np.e**1])
    plt.show()

    print np.log((np.e**-1.9)*np.e**-2.2)