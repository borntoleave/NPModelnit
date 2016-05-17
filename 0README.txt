1.	Run gen_np_trna with parameters
Format: ./gen_np_trna <layer> <molecules per layer>
e.g. ./gen_np_trna 4 3
The tRNAs on each layer are distributed on the vertex of polygon.
The total beads on surface is related to the total charge of the nano particle. While total beads of inbox surface provides an exclusive volume interaction. 
surf_coords.xyz records the coordinates of all the nano particles.
surf_inbox_coords.xyz records the nanoparticle surface within simulation box.
mol_coords.xyz records the coordinates of all the biomolecules.
./fmet_coords.xyz is used to get the biomolecule structure. You can move the filename to command line argument. And may need to change #define MOL_BEAD 300

2.	Run merge_traj with parameters
Format: ./merge_traj <trajectory common title> <division n on theta(0~Pi)>
e.g. ./merge box 6
For convenience I’m using the same file of a cubic box as the trajectory. 
After running the simulation, write a normal expression of the trajectory files, so that pattern pre*post will include all the files. The number will be shown on screen. And transformed coordinates will be saved to TranCoords. The program needs enough trajectories to cover the sphere with certain divisions. So check whether you have more trajectories than required if error occurs. 
division n defines how many sectors you want divide the theta angle. The program will try to divide the phi angle on each layer with the same interval. 

3.	Generate movie
When generating movie, import the surf_coords and the transformed trajectories to VMD. 
If you don’t have a program to convert coordinates to pdb file, you can use my gen_movie. 
Format: ./gen_movie <CGM pdb> <trajectory coordinates> <output pdb> <region[p,r,a]>
If you use this program, you need to generate the trajectories separately. The nano partile never moves, and no frame stepping, so just open it using XYZ format. Each of the trajectories have to be converted to pdb movies. Then import them all in VMD.
