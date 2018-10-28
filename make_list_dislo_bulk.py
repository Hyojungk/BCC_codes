__author__ = 'hyojungkim'

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
                grid.append([int(len(grid)+1),int(entries[0]), float(entries[2]),float(entries[3]),float(entries[4]),int(entries[1])-1])
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

def cut_cell(ref,clamp,cutoff):
    cut_ref=[]
    cut_clamp=[]
    for atom1_i in ref:
        if np.linalg.norm([atom1_i[2],atom1_i[3]]) < cutoff+1.5:
            cut_ref.append(atom1_i)
    for atom2_i in clamp:
        if np.linalg.norm([atom2_i[2],atom2_i[3]]) < cutoff+1.5:
            cut_clamp.append(atom2_i)
    print "number of atoms in cut geomtery: ref VS clamped", len(cut_ref), len(cut_clamp)

    return cut_ref, cut_clamp

def cut_cell_perf(ref,clamp,cutoff,l,layer_num):
    cut_ref=[]
    cut_clamp=[]

    for atom1_i in ref:
            if np.linalg.norm([atom1_i[2],atom1_i[3]]) < cutoff:
                if layer_num==1:
                    if atom1_i[4] >=0.0 and atom1_i[4]< (l/4.0-0.01): cut_ref.append(atom1_i)
                if layer_num==2:
                    if atom1_i[4] >=(l/4.0-0.01) and atom1_i[4]< (l*2/4.0-0.01): cut_ref.append(atom1_i)
                if layer_num==3:
                    if atom1_i[4] >=(l*2/4.0-0.01) and atom1_i[4]< (l*3/4.0-0.01): cut_ref.append(atom1_i)
                if layer_num==4:
                    if atom1_i[4] >=(l*3/4.0-0.01) : cut_ref.append(atom1_i)

    for atom2_i in clamp:
            if np.linalg.norm([atom2_i[2],atom2_i[3]]) < cutoff:
                    cut_clamp.append(atom2_i)

    print "number of atoms in cut geomtery: ref VS clamped", len(cut_ref), len(cut_clamp)
    plot_2D(cut_ref)
    return cut_ref, cut_clamp

def match(grid_unrelaxed, grid_relaxed, l):

    b=[]
    delta=1.3
    print len(grid_unrelaxed)
    for atom1 in grid_unrelaxed:
        count=0
        for atom2 in grid_relaxed:
             ## find the matching atom in the sequence of unrelaxed geometry
             if abs(atom1[2]-atom2[2]) < delta:
                 if abs(atom1[3]-atom2[3]) < delta:
                    if abs(min_image_t(atom2[4],atom1[4],l)-atom2[4]) < delta:
                 #   if abs(atom1[4]-atom2[4]) < 0.8:
                         b.append([atom1[0],atom2[1],atom2[2],atom2[3],min_image_t(atom1[4],atom2[4],l),atom2[5]])
                         count+=1
        if count==0: print 'no matching atom found'
        if count >1: print 'more than one', count
    print len(b)
    return b

def plot_2D(list):
    a,b, x,y,z,c=zip(*list)
    plt.plot(x,y,'bo')
    plt.plot([0],[0],'ro')
    plt.title('2D')
    plt.show()

def match_squence(perf,disl,atomtype,write,layer_num,random_num):
    perf_matched=[]
    disl_matched=[]
    for atom1 in perf:
        for atom2 in disl:
            if atom1[0]==atom2[0]:
                perf_matched.append(atom1)
                disl_matched.append(atom2)
    print "number of final geom of bulk VS dislo", len(perf_matched), len(disl_matched)
    print perf_matched[0], disl_matched[0]
    if write=='True':
        write_geom(perf_matched, 'Fe-V_Perf_matched_%d_R300_random%d.xyz'%(layer_num,random_num),atomtype)
        write_geom(disl_matched, 'Fe-V_Disl_matched_%d_R300_random%d.xyz'%(layer_num,random_num),atomtype)
    return perf_matched, disl_matched

def write_geom(list, title, atomtype):
    f=open('%s'%title,'w')
    f.write('%s'%len(list)+'\n')
    f.write('##%s'%title+'\n')
    for atom in list:
        f.write('%s  %.8f  %.8f  %.8f'%(atomtype[atom[5]],atom[2],atom[3],atom[4])+'\n')
    f.close()

if __name__ =='__main__':

    ### This code is to
    """SET UP"""

    atomic_mass=[50.9415,55.845]  ## mass of V, Fe
    z_mag_Fe_V = 10.184458748504998   # Fe_V
    atomtypes=['V', 'Fe']
    cutoff=30
    layer_num=1
    random_num=16

    ### READ PERFECT BULK
#    grid_perf_unrelaxed=read_geom('sourcecode', open('/Users/Hyojung/anisotropic_2/prac4_Fe-V_t4/perf.xyz','r').read(),atomtypes)
#    grid_perf_relaxed=read_geom('dump', open('/Users/Hyojung/anisotropic_2/lmp_p4_R300/dump.last','r').read(), atomtypes)
    grid_perf_unrelaxed=read_geom('sourcecode', open('/Users/Hyojung/anisotropic_2/prac8_Fe-V_4t_R300/perf','r').read(),atomtypes)
    grid_perf_relaxed=read_geom('dump', open('/Users/Hyojung/anisotropic_2/lmp_R300_random%d/dump.last'%random_num,'r').read(), atomtypes)

    ### READ DISLOCATED GEOMETRY
    grid_disl_unrelaxed= read_geom('sourcecode', open('/Users/Hyojung/anisotropic_2/prac8_Fe-V_4t_R300/screw_10','r').read(),atomtypes)
    grid_disl_relaxed = read_geom('dump',open('/Users/Hyojung/anisotropic_2/lmp_disl_R300_random%d/dump.last'%random_num,'r').read(),atomtypes)
    print "read dislocation geometry!"

    """MATCH WITHIN CUTOFF"""
    grid_perf_unrelaxed, grid_perf_relaxed=cut_cell_perf(grid_perf_unrelaxed, grid_perf_relaxed, cutoff, z_mag_Fe_V,layer_num)
    grid_disl_unrelaxed, grid_disl_relaxed =cut_cell(grid_disl_unrelaxed, grid_disl_relaxed, cutoff)
    print "cut geometry"

    grid_perf_relaxed_squence=match(grid_perf_unrelaxed,grid_perf_relaxed, z_mag_Fe_V)
    grid_dislo_relaxed_squence=match(grid_disl_unrelaxed,grid_disl_relaxed, z_mag_Fe_V)

    ### WRITE THE LIST OF ATOMS IN THE SAME SQUENCE FOR THE BULK AND DISLOCATION

    match_squence(grid_perf_relaxed_squence,grid_dislo_relaxed_squence, atomtypes, 'True',layer_num,random_num)
#    match_squence(grid_perf_relaxed,grid_disl_relaxed,atomtypes,'True',layer_num,random_num)

