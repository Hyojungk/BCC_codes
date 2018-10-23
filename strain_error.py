__author__ = 'hyojungkim'

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
#            grid.append([int(len(grid)+1),0,float(entries[1]),float(entries[2]),float(entries[3]),atomtypes.index(entries[0])])
            grid.append([int(len(grid)+1),0,float(entries[1]),float(entries[2]),float(entries[3]), 0])

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

def cell_list(grid,max_rad,cell_size):
    num_cells = 2*(int(np.ceil(max_rad/cell_size)))
    print num_cells
    cells = [[[] for cols in range(num_cells)] for rows in range(num_cells)]
    for i,atom in enumerate(grid):
        cells[int(np.floor(atom[2]/cell_size))][int(np.floor(atom[3]/cell_size))].append(i)
    return cells

def findatom(nb_mnt,threshold,mycell_row,mycell_col,cells,grid,t_mag):

    found = []
    for g in [-1,0,1]:
        for h in [-1,0,1]:
            ## search for neighbours only in the same cell and surrounding cells
            for atomj in cells[mycell_row+g][mycell_col+h]:
                ## only consider atoms of the appropriate neighbour atom basis type
                t_dist = min(abs(grid[atomj][4]-nb_mnt[2]),abs(grid[atomj][4]-nb_mnt[2]+t_mag),abs(grid[atomj][4]-nb_mnt[2]-t_mag))
                dist=t_dist**2 + (grid[atomj][2]-nb_mnt[0])**2 + (grid[atomj][3]-nb_mnt[1])**2
                if (dist< threshold**2) and (dist>0.00001):
                        found.append(atomj)

    return found

##not used##
def compute_disp_error(grid_perf,grid_disl,cutoff,l):
    count=0
    ddlist = []
    neighbor_num=[]
    for atom1_i,atom2_i in zip(grid_perf,grid_disl):
        if np.linalg.norm([atom2_i[2],atom2_i[3],atom2_i[4]]) < 295:
            count+=1
            n_count=0
            for atom1_j,atom2_j in zip(grid_perf,grid_disl):
                if n_count==8: break;
                if atom1_i[0] != atom1_j[0]:
                    ## only consider diff. disp. between 1st nn within the slab
                    if np.linalg.norm([atom2_j[2]-atom2_i[2],atom2_j[3]-atom2_i[3],min_image_t(atom2_i[4],atom2_j[4],l)-atom2_i[4]]) < cutoff:
                        clamp_d=([atom2_j[2]-atom2_i[2],atom2_j[3]-atom2_i[3],min_image_t(atom2_i[4],atom2_j[4],l)-atom2_i[4]])
                        ref_d=([atom1_j[2]-atom1_i[2],atom1_j[3]-atom1_i[3],min_image_t(atom1_i[4],atom1_j[4],l)-atom1_i[4]])
                        ddlist.append([atom2_i[2:5],atom2_j[2:5], clamp_d, ref_d])
                        n_count+=1
            neighbor_num.append(n_count)

       #     print neighbor_num, ddlist
          #  if count==2: continue

 #   print len(ddlist)
 #   print ddlist[0]
 #   print ddlist[1]
    return ddlist, neighbor_num
###---####

def compute_disp_error_neighbor(grid_perf,grid_disl,cutoff,l,cell_size,max_rad):
    ddlist = []
    neighbor_num=[]
    cells =cell_list(grid_disl,max_rad,cell_size)
    print 'dislo: R%d: completed cell list'%max_rad
    for atom1_i,atom2_i in zip(grid_perf,grid_disl):
  #    if np.linalg.norm([atom1_i[2],atom1_i[3],atom1_i[4]]) < 10:
        mycell_row = int(np.floor(atom2_i[2]/cell_size))
        mycell_col = int(np.floor(atom2_i[3]/cell_size))
        nb_mnt = np.array([atom2_i[2],atom2_i[3],atom2_i[4]])
        ### version1: find neighbors based on the pre-relaxed geometry
    #    found = findatom(nb_mnt,cutoff,mycell_row,mycell_col,cells,grid_perf,l)
        ### version2: find neighbors based on the post-relaxed geometry
        found = findatom(nb_mnt,cutoff,mycell_row,mycell_col,cells,grid_disl,l)

        for n_index in found:  ## neighbor_index
                    ## only consider diff. disp. between 1st nn within the slab
                    clamp_d=([grid_disl[n_index][2]-atom2_i[2],grid_disl[n_index][3]-atom2_i[3],min_image_t(atom2_i[4],grid_disl[n_index][4],l)-atom2_i[4]])
                    ref_d=([grid_perf[n_index][2]-atom1_i[2],grid_perf[n_index][3]-atom1_i[3],min_image_t(atom1_i[4],grid_perf[n_index][4],l)-atom1_i[4]])
                    ddlist.append([atom2_i[2:5],grid_disl[n_index][2:5], clamp_d, ref_d])

        neighbor_num.append(len(found))
    if np.sum(neighbor_num) !=len(ddlist) : print 'ERROR: # of neighbor does not match to ddlist'
    print 'num of atoms:', len(neighbor_num)
 #   print found
 #   print ddlist[0]
 #   print ddlist[1]
 #   print neighbor_num
    return ddlist, neighbor_num

def square_fit(list):
    ## center atom, neighbors, displacement in clamped geometry, displacement in reference (large R) geometry
    atom, atom_n, clamp_d, ref_d=zip(*list)

    ## compute maximum displacement between reference and clamped geometry
    disp=[np.linalg.norm([atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2]]) for atom1,atom2 in zip(clamp_d,ref_d)]
    max_disp=np.max(disp)
    ## compute maximum angle between the displacement vectors in reference and clamped geometry
    angle=[np.arccos(np.dot(a,atom)/(np.linalg.norm(a)*np.linalg.norm(atom)))*180/np.pi for a,atom in zip(clamp_d,ref_d)]
    min_angle=np.min(angle)

 #   print e, residuals
 #   print 'max_disp, min_angle:', max_disp, min_angle
    if min_angle==0: print 'min_angle==0, and min_disp:', np.min(disp)
    clamp_d=np.array(clamp_d)
    ref_d=np.array(ref_d)
    ## compute (linear fit) strain error
    e, residuals, rank, sv= np.linalg.lstsq(clamp_d, ref_d)
    return e, max_disp, np.max(angle)

def square_fit_threshold(list):
    ## center atom, neighbors, displacement in clamped geometry, displacement in reference (large R) geometry
    atom, atom_n, clamp_d, ref_d=zip(*list)

    ## compute maximum displacement between reference and clamped geometry
    disp=[np.linalg.norm([atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2]]) for atom1,atom2 in zip(clamp_d,ref_d)]
    max_disp=np.max(disp)
    ## compute maximum angle between the displacement vectors in reference and clamped geometry
    angle=[np.arccos(np.dot(a,atom)/(np.linalg.norm(a)*np.linalg.norm(atom)))*180/np.pi for a,atom in zip(clamp_d,ref_d)]
    min_angle=np.min(angle)
    clamp=[clamp_d[i] for i in range(0,len(disp)) if disp[i]<0.38]
    ref=[ref_d[i] for i in range(0,len(disp)) if disp[i]<0.38]

    if min_angle==0: print 'min_angle==0, and min_disp:', np.min(disp)
    clamp=np.array(clamp)
    ref=np.array(ref)
    ## compute (linear fit) strain error
    e, residuals, rank, sv= np.linalg.lstsq(clamp, ref)
    return e, max_disp, np.max(angle)

def compute_strain_error(ddlist, neighbor_num):
    index=0
    r_error=[]
#    f=open('bulk_sort_prac.txt', 'w')
    f=open('bulk_strain_error_%d_sorted_threshold.txt'%R_clamp, 'w')
#    f=open('disl_strain_error_%d_sorted.txt'%R_clamp, 'w')
    for i in range(0,len(neighbor_num)):
        N=neighbor_num[i]
        if N<3:
            print "# of neighbor less than 3:" ,  ddlist[index][0]
            break
        if N>9:
            print "# of neighbor larger than 8:" ,  ddlist[index][0]
       #     break

      #  e, max_disp,max_angle=square_fit(ddlist[index:index+N])
        e, max_disp,max_angle=square_fit_threshold(ddlist[index:index+N])
        if max_disp==0: break
        dist=np.linalg.norm([ddlist[index][0][0],ddlist[index][0][1]])
        e_avg=np.sum([e[0][0]-1.0,e[1][1]-1.0],e[2][2]-1.0)/3.0
        r_error.append([dist, e_avg])
        ## save position of center atom x,y,z, the maximum distance and angle between neighbors in clamped geometry
        f.write('%.6f %.6f %.6f %.6f %.6f %.6f'%(dist,ddlist[index][0][0],ddlist[index][0][1],ddlist[index][0][2], max_disp, max_angle) + '\n')
        ## save 3X3 epsilon components : e= I+epsilon where I is identity matrix
        f.write('%.6f %.6f %.6f' %(e[0][0]-1.0,e[0][1],e[0][2])+'\n')
        f.write('%.6f %.6f %.6f' %(e[1][0],e[1][1]-1.0,e[1][2])+'\n')
        f.write('%.6f %.6f %.6f' %(e[2][0],e[2][1],e[2][2]-1.0)+'\n')
        index+=N
    plot(r_error)
    f.close()
    return r_error

def plot(list):
    x,y=zip(*list)
    plt.plot(x,y,'b-')
    plt.xlabel('distance from dislocation core ($\AA$)', fontsize=20)
    plt.ylabel('avg. of strain error', fontsize=20)

    plt.show()

if __name__ =='__main__':

    """SET UP"""

    atomic_mass=[50.9415,55.845]  ## mass of V, Fe
    z_mag_Fe_V = 10.184458748504998   # Fe_V
    atomtypes=['V', 'Fe']
    cutoff=2.8
    R_clamp=300
    ## read in reference: relaxed dislocation with large ~700A cutoff
#    grid_ref = grid_from_xyz(open('n_ref_dislo_R%d.txt'%R_clamp,'r').read(),atomtypes)
    grid_ref = grid_from_xyz(open('ref_bulk_R%d_sorted.txt'%R_clamp,'r').read(),atomtypes)
#    grid_ref = grid_from_xyz(open('ref_bulk_R%d_sorted.txt'%R_clamp,'r').read(),atomtypes)

    ## read in dislocated slab: relaxed dislocation with small cutoff
#    grid_disl = grid_from_xyz(open('n_grid_dislo_R%d.txt'%R_clamp,'r').read(),atomtypes)
    grid_disl = grid_from_xyz(open('grid_bulk_R%d_sorted.txt'%R_clamp,'r').read(),atomtypes)
#    grid_disl = grid_from_xyz(open('grid_bulk_R%d_sorted.txt'%R_clamp,'r').read(),atomtypes)

    """ COMPUTE STRAIN ERRORS BETWEEN CLAMPED GEOMETRY AND REFERENCE (LARGE R) """

 ##   ddlist, neighbor_num = compute_disp_error(grid_ref,grid_disl,cutoff,z_mag_Fe_V)
    cell_size=15
    ddlist, neighbor_num =compute_disp_error_neighbor(grid_ref,grid_disl,cutoff,z_mag_Fe_V,cell_size,R_clamp)
    compute_strain_error(ddlist, neighbor_num)

