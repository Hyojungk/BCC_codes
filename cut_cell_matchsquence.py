__author__ = 'Hyojung'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager,rc

def read_geom(filetype, s, atomtypes):
    """
    :param filetype: 'sourcecode' or 'xyz' works
    :param s: file name to open
    :param atomtypes: ['V', 'Fe'] or ['Ti','Mo']
    :return: [ index, region, m, n, t, basis] where index=1,2...,N , region is 0 as default, basis are atom type: V=0, Fe=1.
    """
    grid=[]
    if filetype=='sourcecode':
        for line in s.splitlines()[2:]:
            if line !='':
                entries = line.split()
                grid.append([int(len(grid)+1),0,float(entries[1]),float(entries[2]),float(entries[3]),atomtypes.index(entries[0])])
        return grid
    if filetype=="dump":
        for line in s.splitlines()[9:]:
            if line != '':
                entries = line.split()
                grid.append([int(len(grid)+1),int(entries[1])-1, float(entries[2]),float(entries[3]),float(entries[4]),int(entries[1])-1])
                          ## ,float(entries[5]),float(entries[6]),float(entries[7])])
        return grid

def min_image_t(tcoord_i,tcoord_j,l):

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

    if abs(tcoord_j-tcoord_i+l) < abs(tcoord_j-tcoord_i):
        tcoord_j += l
    elif abs(tcoord_j-tcoord_i-l) < abs(tcoord_j-tcoord_i):
        tcoord_j -= l

    return tcoord_j

def cut_cell(ref,cutoff,l):
    cut_ref=[]
    for atom1_i in ref:
        if np.linalg.norm([atom1_i[2],atom1_i[3]]) < cutoff:
      #          if atom1_i[4] >=0.0 and atom1_i[4]<l/4.0-0.01:
                    cut_ref.append(atom1_i)

    print "number of atoms in cut geomtery: ref VS clamped", len(cut_ref)
 #   plot_2D(cut_ref)
    return cut_ref

def flatten_position(list,unique_xy,l):
    ## [ index, region, m, n, t, basis]
    print 'unique_num', len(unique_xy)

    sorted_list=[]
    delta=1.1
    for unique in unique_xy:
        same_xy=[]
        for atoms in list:
            if np.abs(atoms[2]-unique[2]) < delta:
                if np.abs(atoms[3]-unique[3]) < delta:
                        same_xy.append(atoms)
        if len(same_xy)!=4 or len(same_xy)>4: print "not 4", len(same_xy)
        if len(same_xy)==4: sorted_list.append(same_xy)
  #  plot_2D(unique_xy)
    print sorted_list[0][1]
    final=[]
    for i in range(0,len(sorted_list)):
        ## compute the average of positions along z
        sum=0.0
        for k in range(0,4):
            sum+=sorted_list[i][k][4]
        final.append([sorted_list[i][0][0],sorted_list[i][0][1], sorted_list[i][0][2],sorted_list[i][0][3],sum/4.0,0])
 #   plot_2D(final)
    return final

def plot_2D(list):
    a,b, x,y,z,c=zip(*list)
    plt.plot(x,y,'bo')
    plt.plot([0],[0],'ro')
    plt.title('2D')
    plt.show()

def write_geom(list, title, atomtype):
    f=open('%s'%title,'w')
    f.write('%s'%len(list)+'\n')
    f.write('##%s'%title+'\n')
    for atom in list:
        f.write('%s  %.8f  %.8f  %.8f'%(atomtype[atom[5]],atom[2],atom[3],atom[4])+'\n')
    f.close()

def match_squence_xy(perf,disl,atomtype,random_num):
    perf_matched=[]
    disl_matched=[]
    delta=1.1
    for atom1 in perf:
        for atom2 in disl:
            if np.abs(atom1[2]-atom2[2]) < delta:
                if np.abs(atom1[3]-atom2[3]) < delta:
                        perf_matched.append(atom1)
                        disl_matched.append(atom2)

    print "number of final geom of bulk VS dislo", len(perf_matched), len(disl_matched)
    print perf_matched[0], disl_matched[0]

    write_geom(perf_matched, 'Fe-V_Perf_matched_flattened_random%d.xyz'%random_num,atomtype)
    write_geom(disl_matched, 'Fe-V_Disl_matched_flattened_random%d.xyz'%random_num,atomtype)
    return perf_matched, disl_matched

if __name__ =='__main__':

    """SET UP"""

    z_mag_Fe_V = 10.184458748504998   # Fe_V
    atomtypes=['V', 'Fe']
    random_num=16
    cutoff=30
    ### READ PERFECT BULK
    ### READ PERFECT BULK
#    grid_perf_unrelaxed=read_geom('sourcecode', open('/Users/Hyojung/anisotropic_2/prac4_Fe-V_t4/perf.xyz','r').read(),atomtypes)
#    grid_perf_relaxed=read_geom('dump', open('/Users/Hyojung/anisotropic_2/lmp_p4_R300/dump.last','r').read(), atomtypes)
#    grid_perf_unique=read_geom('sourcecode', open('Fe-V_Perf_matched_1_R300.txt','r').read(), atomtypes)
    grid_perf_relaxed=read_geom('dump', open('/Users/Hyojung/anisotropic_2/lmp_R300_random%d/dump.last'%random_num,'r').read(), atomtypes)
    grid_perf_unique=read_geom('sourcecode', open('Fe-V_Perf_matched_1_R300_random%d.xyz'%random_num,'r').read(), atomtypes)

    ### READ DISLOCATED GEOMETRY
#    grid_disl_unrelaxed= read_geom('sourcecode', open('/Users/Hyojung/anisotropic_2/prac4_Fe-V_t4/screw_10.xyz','r').read(),atomtypes)
#    grid_disl_relaxed = read_geom('dump',open('/Users/Hyojung/anisotropic_2/lmp_p4_dislo_R300/dump.last','r').read(),atomtypes)
#    grid_disl_unique=read_geom('sourcecode', open('Fe-V_Disl_matched_1_R300.txt','r').read(), atomtypes)
    grid_disl_relaxed = read_geom('dump',open('/Users/Hyojung/anisotropic_2/lmp_disl_R300_random%d/dump.last'%random_num,'r').read(),atomtypes)
    grid_disl_unique=read_geom('sourcecode', open('Fe-V_Disl_matched_1_R300_random%d.xyz'%random_num,'r').read(), atomtypes)

    print "read dislocation geometry!"

    cut_perf= cut_cell(grid_perf_relaxed,cutoff,z_mag_Fe_V)
    cut_disl= cut_cell(grid_disl_relaxed,cutoff,z_mag_Fe_V)
    perf=flatten_position(cut_perf,grid_perf_unique,z_mag_Fe_V)
    disl=flatten_position(cut_disl,grid_disl_unique,z_mag_Fe_V)
 #   write_geom(perf, 'flatted_perf_Fe_V_random1.xyz', atomtypes)
 #   write_geom(disl, 'flatted_disl_Fe_V_random1.xyz', atomtypes)

    match_squence_xy(perf,disl,atomtypes,random_num)
  #  plot_2D(cut_perf_ref)