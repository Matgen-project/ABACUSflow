"""
This is a simple model database.

Created on Tue May  1 22:36:18 2018

@author: Shen Zhen-Xiong
"""

def graphene():
    import numpy as np
    from ase import Atoms
    p = np.array([[0, 0, 0],
                  [0.0, 0.0, 3.35],
                  [0.0, 2.84, 0.0],
                  [0.0, 1.42, 3.35]])
    c = np.array([[1.23, 2.13, 0.0],
                  [-1.23, 2.13, 0.0],
                  [0, 0, 6.7]])
    print(c)
    return  Atoms('4C',positions=p, cell=c, pbc=(1,1,1))

def magic_angle_graphene(index = 0):
    index = index
    lattice_a0 = 2.46
    lattice_c0 = 3.35
####################################################################################################
##########  Calculate length, angle,number of atoms 
####################################################################################################    
    import math
    cos_cta = (int(index)**2 * 3 + 3 * int(index) + 0.5) / (int(index)**2 * 3 + int(index)*3+1)
    angle_deg = math.degrees(math.acos(cos_cta))
    superlattice_length = math.sqrt(int(index)**2 * 3 + int(index)*3+1) * lattice_a0
    natom = 4 *(3 * index ** 2 + 3 * index + 1)
    print('#'*70)
    print('-'*5,'Index', '-'*5,'Length','-'*13,'Angle(degree)','-'*6,'Atoms')
    print('-'*5,index,'-'*6,superlattice_length,'-'*3,angle_deg,'-'*3,natom)
    print('#'*70)
    print('\n\n')
          
####################################################################################################
##########   Create  unit cells of two layers
####################################################################################################
    import numpy as np
    from ase import Atoms
    
    pos_layer1_unit = np.array([[0, 0, 0],
                            [0.0, 2.84, 0.0]])
    pos_layer2_unit = np.array([[0.0, 0.0, 3.35],
                            [0.0, 1.42, 3.35]])
    cell_unit = np.array([[-1/2 * lattice_a0, math.sqrt(3)/2 * lattice_a0, 0],
                          [1/2 * lattice_a0, math.sqrt(3)/2 * lattice_a0, 0], 
                          [0.0,0.0,lattice_c0 * 2]])

    layer1_unit = Atoms('2C',positions=pos_layer1_unit, cell=cell_unit, pbc=[1,1,0])
    layer2_unit = Atoms('2C',positions=pos_layer2_unit, cell=cell_unit, pbc=[1,1,0])

####################################################################################################
##########   Create  super cells of two layers
####################################################################################################
    from ase.build import make_supercell
    cell_super_layer1 = np.array([[index, index + 1, 0],
                                  [-(index + 1),(2 * index + 1),0],
                                  [0,0,1]]) 
    cell_super_layer2 = np.array([[index + 1, index, 0],
                                  [-index ,(2 * index + 1),0],
                                  [0,0,1]])

    super_layer1 = make_supercell(layer1_unit,cell_super_layer1)
    super_layer2 = make_supercell(layer2_unit,cell_super_layer2)
      
####################################################################################################
##########   Create  magic angle twisted bilayer graphene
####################################################################################################   
    
    mat_cell_unit = np.mat(cell_unit)

    mat_cell_super_layer1 = np.mat(cell_super_layer1) * mat_cell_unit
    mat_cell_super_layer2 = np.mat(cell_super_layer2) * mat_cell_unit

    cell_super = np.array([[-0.5 * superlattice_length,math.sqrt(3)/2 * superlattice_length,0],
                     [0.5 * superlattice_length, math.sqrt(3)/2 * superlattice_length,0],
                     [0,0,6.7]])

    mat_cell_super = np.mat(cell_super)
    mat_pos_layer1_super = np.mat(super_layer1.get_positions())
    mat_pos_layer2_super = np.mat(super_layer2.get_positions())

    mat_pos_layer1_super_newcell = mat_pos_layer1_super * mat_cell_super_layer1.I *mat_cell_super
    mat_pos_layer2_super_newcell = mat_pos_layer2_super * mat_cell_super_layer2.I *mat_cell_super

    super_layer1.set_cell(np.array(mat_cell_super))
    super_layer1.set_positions(np.array(mat_pos_layer1_super_newcell))
    super_layer2.set_cell(np.array(mat_cell_super))
    super_layer2.set_positions(np.array(mat_pos_layer2_super_newcell))

    interface = super_layer1.copy()
    interface.extend(super_layer2)
    interface.center(vacuum = 20,axis = 2)

    return interface
    

def magic_angle_graphene_cluster(index = 1,corner=0):
    n = index
    corner = corner
    c = 1.59
    
    from ase import Atoms,Atom
####################################################################################################
##########   Create  unit cell
####################################################################################################    
    graphene_unit = Atoms('4C',positions=[[0.0,0.0,0.0],
                                          [0.0,0.0,3.35],
                                          [0.0,1.42,0.0],
                                          [0.0,1.42,3.35]],
                              cell = [[1.23,2.13,0.0],
                                      [-1.23,2.13,0.0],
                                      [0.0,0.0,6.7]],
                         pbc=(1,1,0))
    

####################################################################################################
##########   Create  super cell
####################################################################################################              
    graphen_natom_C = 12 * n * n
    graphen_natom_H = 12 * n
    graphen_cell = [2.46*(2*n-1)+15,0.71*(6*n-4)+15+2*c,21]
    graphen_cell_center=((2.46*(2*n-1)+15)*0.5,(0.71*(6*n-4)+15+2*c)*0.5,12.175)
    
    graphene_super = graphene_unit.repeat((2*n+1,2*n+1,1))
    #view(graphene_super)
    
    print(' Number of carbon:   ',graphen_natom_C,'\n',
      'Number of hydrogen: ', graphen_natom_H,'\n',
      'Number of atoms:    ',graphen_natom_C+graphen_natom_H)
    
####################################################################################################
##########   Create  cluster without passivation
####################################################################################################     
    graphen_cluster=Atoms()
    for atom in graphene_super:
        if atom.x >(-1.23*(2*n-1)-0.01) and atom.x <(1.23*(2*n-1)+0.01):
            if atom.y >(0.71*(3*n+4)-0.01) and atom.y <(0.71*(9*n+4)+0.01):
                line1 = 4.23+2.13/1.23*atom.x-atom.y
                line2 = 4.23-2.13/1.23*atom.x-atom.y
                if line1 < 0 and line2 < 0:
                    graphen_cluster.append(atom)

####################################################################################################
##########   Create  passivation
####################################################################################################     
    graphen_passivation=Atoms()
    for atom in graphen_cluster:
        if atom.y< 2.13/1.23*atom.x+4.26+0.02:
            graphen_passivation.append(Atom('H',
                                        position=(atom.x+c*0.866,atom.y-c*0.5,atom.z)))
        elif atom.y<-2.13/1.23*atom.x+4.26+0.02:
            graphen_passivation.append(Atom('H',
                                        position=(atom.x-c*0.866,atom.y-c*0.5,atom.z)))
        elif atom.y>2.13/1.23*atom.x+4.26+0.71*(12*n-4)-0.02:
            graphen_passivation.append(Atom('H',
                                        position=(atom.x-c*0.866,atom.y+c*0.5,atom.z)))
        elif atom.y>-2.13/1.23*atom.x+4.26+0.71*(12*n-4)-0.02:
            graphen_passivation.append(Atom('H',
                                        position=(atom.x+c*0.866,atom.y+c*0.5,atom.z)))
        elif atom.y<0.71*(3*n+5)+0.02:
            graphen_passivation.append(Atom('H',
                                        position=(atom.x,atom.y-c,atom.z)))
        elif atom.y>0.71*(9*n+3)-0.02:
            graphen_passivation.append(Atom('H',
                                        position=(atom.x,atom.y+c,atom.z)))
            
    graphen_cluster.extend(graphen_passivation)
    
####################################################################################################
##########   Create  magic angle graphene cluster
#################################################################################################### 
    graphen_cluster.set_cell(graphen_cell)
    graphen_cluster.center()
    graphen_cluster_corner = graphen_cluster.copy()
    graphen_cluster_layer2 = graphen_cluster.copy()
    del graphen_cluster_corner[[atom.index for atom in graphen_cluster_corner if atom.z > 10.5]]
    del graphen_cluster_layer2[[atom.index for atom in graphen_cluster_layer2 if atom.z < 10.5]]
    graphen_cluster_layer2.rotate(a=corner,v='z',center=graphen_cell_center)
    graphen_cluster_corner.extend(graphen_cluster_layer2)
    
    return graphen_cluster_corner
    




    







          
         
    
    