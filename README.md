This is the base for a project on packing charged multibody systems with implicit
solvation. This scripts are used for the extraction of statistical dada from large samples of 
crystal structures, selected from the PDB. The generated statistics are used on  two type of libraries
for multibody dynamics: a library with statistical potentials for pairs of dihedral angles 
in the main chain local conformation,  and a rotamer libraries. This  libraries are defined as a 
set of discrete probability distribution, and used to generat
internal forces multibody systems.

Examples of statistical potentials generated from two dihedral angles on a protein main chain:
![Image](../images/IAGE2.gif)

Examples of statistical potentials generated between dihedral angles and rotamers in an amino acid:
![Image](../images/IAGE1.gif)
