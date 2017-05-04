---
layout: post
title:  "Coarse-grained models"
date:   2017-04-30 01:02:47 +0100
categories: update
---
#Protein structure

The method described in this work will be evaluated for the kinematic analysis of polypeptide chains. 
Biological proteins are polymeric chains build from amino acid monomers. These amino acids contain 
five chemical components: a central $\alpha$-carbon ($C_\alpha$), an $\alpha$-proton ($H$), an
 \textit{amino functional} group ($-NH$), a \textit{carboxylic acid functional} group ($-COOH$), a \textit{side chain} 
group ($R$) (see Figure \ref{fig1}) \cite{McCammon 1984}. The residual side chain group differentiates the 
common biological amino acids, and is the main factor of the peptide chains local stable conformations.
 These amino acids combine to become proteins through an energy-driven combination. This result in the
 creation of a peptide bond between the two amino acids, and repeating the process creates a polypeptide
 containing several peptide bonds. These peptide bonds behave like a partial double bond, which have restricted
 rotation about the bond. This restriction results in a stable peptide plane. These peptide planes are repeating units 
that exhibit constant structures in the protein and reduce the number of degrees of freedom. The polypeptide chain 
is intrinsically flexible because many of the covalent bonds that occur in its backbone and side chains are rotationally 
permissible. A protein can be defined by one or more polypeptide chains.

![Image](../../../../images/aminoAcid.eps)
a) Amino acids are molecules containing an amine group ($H_2N-$), a carboxylic acid group ($-COOH$), 
and a side chain ($R$) that is pecific to each amino acide. The first carbon that attaches to a functional group is named 
apha-carbon ($C_\alpha$).

![Image](../../../../images/peptide.eps)
(b) Every peptide has a $N$-terminus residue and a $C$-terminus residue on the ends of the 
peptide (Source Wikipedia).

Geometric relationship involving atoms in the polypeptide fully define a thee-dimensional proteins stable conformation. 
The relationships consist of bond lengths, bond angles, dihedral angles and improper dihedral angles. The primary 
contributions from these parameters, which determine overall polypeptide structure, are the dihedral angles. Typically, 
the peptide plane remains relatively rigid during protein dynamics such that the bond lengths and bond angles remain 
constant, due the large energy cost for its deformation. As a result, the dihedral angles are the essential degrees of 
freedom that dictate the position of the polypeptide backbone atoms, defining the protein secondary structure 
\cite{McCammon 1984}. 

![Image](../../../../images/Trialanine.eps)
Backbone dihedral angles in the molecular structure of trialanine (Altis 2008)


#Coarse-grained models

In a polypeptide chain the torsional motion is predominantly local in character. Therefore its model is here simplified as
 a constrained multibody system, and the overall dynamic is described by backbone dihedral angles, and possible linear
 elastic deformation allowing only covalent bonds lengths fluctuations. This type of model is not new, it have been used 
on molecular coarse-grained simulations, e.g. in systems like GROMACS \cite{Hess 2008},  widely used for long time 
simulations, where however the system dynamics is constrain by force fields like MARTINI force field \cite{Bekker 1990}. 
The most important feature of this type of approach is the possibility of modelling protein dynamics without explicitly
 treating every atom in the problem. Using this quasi-continuum approach, must degrees of freedom are eliminate, and
 force or energy calculations are largely expedited. Here however we are interested in analysing the possible 
advantages of imposing kinematic constrains on the system using kinematic joints. Since the definition of such 
constrains increase the simulation computational complexity, we are concern on its numerical stability for long 
time simulations \cite{Bekker 1990}. 

![Image](../../../../images/multibody.eps)
Refinement for a protein with $n$ amino acids split in $2n+1$ bodies, with its $2n$ dihedral angles identified with
 $\beta_1,\cdots,\beta_{2n}$


The methodology described bellow was evaluated on two types of  coarse-grained models, defined using 
different kinematic joints. The first is defined by an open chain of rigid bodies, linked together using revolute joints.
 For the second model we used the same open chain system but, to allow some fluctuations in the system structure,
 the components of the system are linked using cylindrical joints. In Fig. \ref{fig3}, presents the coarse-grained model 
for a generic polypeptide chain defined by $n$ amino acid, using a constrained multibody system with $N=2n+1$ bodies 
$b_1,b_2,\dots,b_{2n+1}$, with dihedral angles $\beta_1,\beta_2,\dots,\beta_{2n}$. Note that, each amino acid is split
 between three bodies, and each body is linked to the next by a covalent bonds. This gives great flexibility to them and
 its conformation can be characterized by dihedral angles on this bonds. The values of these angles are not uniformly 
distributed, they have strong local propensity  and are highly correlated \cite{Betancourt 2004}. 

The use of simplified models reduces the complexity of the interactions and hence reduces the amount of computation 
involved. This shouldn't be seen as a limitation. Software tools for reconstructing all-atoms structure from backbone
 structures (e.q. PHOENIX\footnote{http:$\backslash\backslash$cbsu.tc.cornell.edu/software/photarch}  and
 Maxsprout \footnote{http:$\backslash\backslash$www.ebi.ac.uk/maxsprout} ) are often employed to complete
 folding trajectory with side chain information.
