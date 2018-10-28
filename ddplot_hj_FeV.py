import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager,rc


def grid_from_xyz(s,atomtypes):

    """
    Read from a string containing the data from the anisotropic dislocation geometry file.

    Parameters
    ----------
    s    : string containing the data from the anisotropic dislocation geometry setup file

    Returns
    -------
    grid : list of [atom index,region,m-coord,n-coord,t-coord,basis] for each atom
           'basis' is 0,1,... indicating the different atom types

    """

    grid = []
    for line in s.splitlines()[2:]: ## skip first 2 lines of xyz file which are comments
        if line != '':
            entries = line.split()
            grid.append([int(len(grid)+1),0,float(entries[1]),float(entries[2]),float(entries[3]),atomtypes.index(entries[0])])

    return grid


def min_image_t(tcoord_i,tcoord_j,t_mag,a0):

    """
    Finds the minimum image of point j, relative to point i, in the periodic t direction

    Parameters
    ----------
    tcoord_i : t-coord of point i
    tcoord_j : t-coord of point j
    t_mag: periodic length (slab thickness) (in units of a0)
    a0 : lattice constant (in Angstroms)

    Returns
    -------
    tcoord_j : t-coord of the min image of point j

    """

    if abs(tcoord_j-tcoord_i+t_mag*a0) < abs(tcoord_j-tcoord_i):
        tcoord_j += t_mag*a0
    elif abs(tcoord_j-tcoord_i-t_mag*a0) < abs(tcoord_j-tcoord_i):
        tcoord_j -= t_mag*a0

    return tcoord_j


def compute_diffdisp(grid_perf,grid_disl,t_mag,a0,cutoff,modulo,plot):

    """
    computes the differential displacement between pairs of 1st nn atoms

    Parameters
    ----------
    grid_perf : list of [atom index,region,m-coord,n-coord,t-coord,basis] for perfect slab
    grid_disl : list of [atom index,region,m-coord,n-coord,t-coord,basis] for dislocated slab
    t_mag: periodic length (slab thickness) (in units of a0)
    a0 : lattice constant (in Angstroms)
    cutoff : neighbor distance cutoff (in Angstroms); include only 1st nn in the slab
    modulo : differential displacement goes from -modulo/2 to modulo/2

    Returns
    -------
    ddlist : list of [atom i coords, atom j coords, diff. disp. vector between atoms i and j]
             for each pair of 1st nn atoms i,j

    """

    disps = [[atom2[2]-atom1[2],atom2[3]-atom1[3],min_image_t(atom1[4],atom2[4],t_mag,a0)-atom1[4]] for atom1,atom2 in zip(grid_perf,grid_disl)]

    ddlist = []
    for atom1_i,atom2_i,disp_i in zip(grid_perf,grid_disl,disps):
        for atom1_j,atom2_j,disp_j in zip(grid_perf,grid_disl,disps):
            if atom1_i[0] > atom1_j[0]:
                ## only consider diff. disp. between 1st nn within the slab
                if np.linalg.norm([atom2_j[2]-atom2_i[2],atom2_j[3]-atom2_i[3],min_image_t(atom2_i[4],atom2_j[4],t_mag,a0)-atom2_i[4]]) < cutoff:
                    diffdisp = [disp_j[0]-disp_i[0],disp_j[1]-disp_i[1],(disp_j[2]-disp_i[2])%modulo]
                    if plot == 'screw':
                    ## take care of apparent discontinuity in the screw component across the (negative) cut plane
                        if diffdisp[2] > modulo/2: diffdisp[2]  -= modulo
                        elif diffdisp[2]  < -modulo/2: diffdisp[2]  += modulo
                    ddlist.append([atom2_i[2:5],atom2_j[2:5],diffdisp])

    return ddlist


if __name__ == '__main__':


    """ SETUP """

    ## you'll change these inputs for your system...
    atomtypes = ['Fe','V']
    a0 = 2.940
    t_mag = np.sqrt(3)/2
    print 'b:', a0*t_mag
    layer_num=1
    random_num=16
#    fig_title='randomness_v%d: layer%d'%(random_num,layer_num)
    fig_title='randomness_v%d: flattened'%random_num

    rnge=10
    ## read in perfect slab
    #    grid_perf = grid_from_xyz(open('Fe-V_Perf_matched_%d_R300.txt'%layer_num,'r').read(),atomtypes)
#    grid_perf = grid_from_xyz(open('Fe-V_Perf_matched_%d_R300_random%d.xyz'%(layer_num,random_num),'r').read(),atomtypes)

    grid_perf= grid_from_xyz(open('Fe-V_Perf_matched_flattened_random%d.xyz'%random_num,'r').read(),atomtypes)
#    grid_perf = grid_from_xyz(open('/Users/hyojungkim/Desktop/anisotropic_2018/prac3_Fe-V_hardcore/perf.xyz','r').read(),atomtypes)

    ## read in dislocated slab
 #   grid_disl = grid_from_xyz(open('Fe-V_Disl_matched_%d_R300.txt'%layer_num,'r').read(),atomtypes)
 #   grid_disl = grid_from_xyz(open('Fe-V_Disl_matched_%d_R300_random%d.xyz'%(layer_num,random_num),'r').read(),atomtypes)

    grid_disl = grid_from_xyz(open('Fe-V_Disl_matched_flattened_random%d.xyz'%random_num,'r').read(),atomtypes)
#    grid_disl = grid_from_xyz(open('/Users/hyojungkim/Desktop/anisotropic_2018/prac3_Fe-V_hardcore/screw_10.xyz','r').read(),atomtypes)
    ## your dislocated slab may be in a different format out of lammps
    ## you'll have to convert it or write a function to read from the dump file

    ## you'll probably want to carve out a smaller region of interest, i.e. nearer the core
    ## otherwise to compute and plot the dd map for 10s of thousands of atoms would be slow


    """ COMPUTE DIFFERENTIAL DISPLACEMENTS """

    plot = 'screw'
    ## you'll change these inputs for your system... modulo _screw = length of burgers vector
    modulo_screw = a0*t_mag

    if plot == 'screw': ## the screw component of each partial
        ddlist = compute_diffdisp(grid_perf,grid_disl,t_mag,a0,cutoff=3.0,modulo=modulo_screw,plot=plot)

    """ PLOTTING """

    rc('font',family='Times New Roman',size=15)
    ticks_font = font_manager.FontProperties(family='Times New Roman', style='normal',
    size=15, weight='normal', stretch='normal')
    axis_font = {'fontname':'Times New Roman','fontsize':15}

    fig = plt.figure(figsize=(5.0,5.0))
    ax1 = fig.add_subplot(111)

    index,reg,m,n,t,basis = zip(*grid_disl)
    coordi,coordj,ddij = zip(*ddlist)
    coordi_m,coordi_n,coordi_t = zip(*coordi)
    coordj_m,coordj_n,coordj_t = zip(*coordj)
    ddij_m,ddij_n,ddij_t = zip(*ddij)


    """ PLOT DD MAP """

    if plot == 'screw':
#        ax1.set_title('layer%d'%layer_num,**axis_font)
        ax1.set_title('%s'%fig_title,**axis_font)

        m_chem1=[m[i] for i in range(0,len(basis)) if basis[i]==0]
        n_chem1=[n[i] for i in range(0,len(basis)) if basis[i]==0]
        m_chem2=[m[i] for i in range(0,len(basis)) if basis[i]==1]
        n_chem2=[n[i] for i in range(0,len(basis)) if basis[i]==1]

#        ax1.scatter(m_chem1,n_chem1,s=30,linewidth=0.5,facecolor='white',edgecolors='k')
        ax1.scatter(m_chem1,n_chem1,s=30,linewidth=0.5,facecolor='orange',edgecolors='k')
        ax1.scatter(m_chem2,n_chem2,s=30,linewidth=0.5,facecolor='gray',edgecolors='k')

#        ax1.scatter(m,n,s=30,linewidth=1.0,facecolor='None',edgecolors='k')
        for mi,ni,mj,nj,ddt in zip(coordi_m,coordi_n,coordj_m,coordj_n,ddij_t):
            midpt = np.array([(mi+mj)/2,(ni+nj)/2])
          #  scaling = ddt/modulo_screw
            scaling= ddt/(a0*t_mag/3)*0.9
            arrowvec = scaling*np.array([mj-mi,nj-ni])
            if abs(arrowvec[0]) == 0: arrowvec[0] = 1e-16  ## if I leave it to be zero, I think it gives a plotting error
            if abs(arrowvec[1]) == 0: arrowvec[1] = 1e-16
            arrowtail = midpt - 0.5*arrowvec
            ax1.arrow(arrowtail[0],arrowtail[1],arrowvec[0],arrowvec[1],
                      head_width=0.2,head_length=0.2,length_includes_head=True,fc='k',ec='k')

    ## labeling the axes, etc.
    ax1.set_xlim(-rnge,rnge)
    ax1.set_ylim(-rnge,rnge)
 #   ax1.set_xticks([-20,-10,0,10,20])
 #   ax1.set_yticks([-5,0,5])
#    ax1.set_xticklabels([r'$-20\rm{\AA}$',r'$-10\rm{\AA}$',r'$0$',r'$10\rm{\AA}$',r'$20\rm{\AA}$'],fontsize=10)
#    ax1.set_yticklabels([r'$-5\rm{\AA}$',r'$0$',r'$5\rm{\AA}$'],fontsize=10)
    ax1.set_xlabel(r'$[11\bar2]$',labelpad=7,**axis_font)
    ax1.set_ylabel(r'$[\bar110]$',labelpad=-13,**axis_font)
#    ax1.set_aspect('equal',adjustable='box')
#    fig.subplots_adjust(bottom=0.24,top=0.85,left=0.24,right=0.73,wspace=0.00)

    plt.show()
#    plt.savefig('fullrelax/Ni_Nye+dd_screw.pdf')

