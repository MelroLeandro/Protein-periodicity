---
layout: post
title:  "Protein local propensity"
date:   2017-04-30 01:02:47 +0100
categories: update
---
# Protein local propensity

Polypeptide chains are known to have strong local propensities, well characterized by the existence of specific amino 
acid patterns in chain with predefined ranges of dihedral angles. A commune strategy to study the angular propensity 
is through \textit{Ramachandran maps} used to produce distributions of the dihedral angles and their probabilities, 
extracted from statistical libraries with local propensities generated using the PDB. Studies of this plots 
\cite{McCammon 1984} show that they reflect the local interactions of free energy. The type of conformation are 
determined from a balance between local interactions (those closed to the sequence) and non-local ones. 
Effective statistical potentials can be extract from these populations have been studied for over 40 years. 
These potentials are based on the correlation of the observed frequency of a structure with its associated 
free energy. Thus, those potentials have a global minimum corresponding to the must frequency observed 
native informations (or collections of substructures found must often in native conformations). Such local potentials 
have been combined with non-local potential to predict protein folding structure, used  by the state-of-the-art systems
 for conformation perdition like ROSETTA.  

# Van Mises propensity distributions

Ramachandran maps play a central role in developing empirical energy functions for structure prediction and simulation 
\cite{amakrishann 1965}\cite{Shapovalov 2011}. Those are used as a probability density function gives the probability of 
finding an amino acid conformation in a specific range of $\phi$, $\psi$ values. For this work these functions are given on 
a $2^\circ\times 2^\circ$ grid from $-180^\circ$ to $180^\circ$ in $\phi$, $\psi$ values. Such distributions was derived for
 each amino acid types, but since they are very irregular we used the van Mises density estimations to have soother 
surfaces. The quality and quantity of the data are crucial in determining distributions to approximate the system free 
energy and its gradients. Here all the distribution are computed using a datasets of $254$ globular proteins with a
 resolution cutoffs of 1.0 $\AA$, filters to reduce the number of amino acids angles evaluated over great structural 
stresses. Those distribution estimations were used to build a library of conformational propensities for each amino 
acid, applied for drive polypeptide chain kinematics. Examples of van Mises propensity distributions for two amino 
acids are presented in Figure \ref{fig10}. The library was developed by combining a prior estimate of the probabilities
 of each $(\phi,\psi)$ bin raw counts by amino acid. This standard probability distributions are quite bumpy in their 
variation, a result of using raw counts in the probability estimates and calculation of simple averages. In order to 
produce smooth and continuous estimates of the conformation probabilities, we use kernel density estimation. A 
kernel is a nonnegative symmetric function that integrates to $1.0$ and is centered on each data point. Density estimates
 at specific query points are determined by summing the values of the kernel functions centered on the data points.
 The smoothness of the density estimate is determined by the form of the kernel, in particular its bandwidth. For each 
amino acid, $aa$, we determine a probability density estimate, $P(\phi,\psi | aa)$. For that \textit{von Mises probability 
density function} (PDF) as the kernel are used since this density estimates are more appropriately for angles than the
 usual Gaussian or other nonperiodic kernels \cite{Mardia 2000}. Because Ramachandran probability density is defined
 for the backbone torsion angles $\phi$ and $\psi$ as two arguments, we use a nonadaptive kernel density estimators 
in two dimensions written as the sum over products of $\phi$ and $\psi$ van Mises kernels for $N_r$ data points of 
amino acid of type, $aa$:

$$P(\phi,\psi | aa)=\frac{1}{4\pi^2N_r}\sum_{i=1}^{N_r}\frac{1}{I_0(\kappa)^2}\exp(\kappa \cos(\phi-\phi_i)+\kappa \cos(\psi-\psi_i))$$

In this case, $\sqrt{\frac{1}{\kappa}}$ defines a radius of the two dimensional hump covering $60\%$ of the kernel 
density \cite{Shapovalov 2011}, $I_0$ is the Bessel function of the first kind of order $0$, normalizing the kernel to $1$. 


![Image](../../../../images/Rama1_test2ALA.jpg)
Amino acid Mises propensity distributions  ALA

![Image](../../../../images/Rama1_test2PHE.jpg)
 Amino acid Mises propensity distributions PHE 

To reduce the amount of computation involved on the protein kinematics analysis, the library of conformational
 propensities was completed with directional gradients, for each distribution. For that the central difference method
 was used.

Amino acid dihedral angles cannot take any arbitrary values due to atomic clashes and orientations 
\cite{Shapovalov 2011}. It is also consensual \cite{Bosco 2003} on the basis of the conformational enumeration of
 polypeptide chains and molecular dynamic simulations that the Ramachandran basin populations are affected by
 their nearest neighbours. The populations are affected, in particular, by the neighbour's conformation and their identity.
 There are a strong correlations between a residue's conformation  and that of its neighbours are responsible for
 cooperative effects. This works however uses only correlation defined by dihedral angles in the same amino acid and 
between two consecutive amino acid. On the presented extended mechanical model, the coupling band of springs,
 uses these potentials to find optimal conformations, by the generation of body torsions along each associated covalent
 bond.

The dependence of an amino acid conformation on the conformations of its adjacent residues involves too many 
variables to be captured in a single probability density function for the available data. Instead, we divided the 
probability densities in individual terms involving pairs of angles. In particular, for dihedral angles $\phi_i$, $\psi_i$, 
$\phi_{i+1}$, and $\psi_{i+1}$ in consecutive amino acids
we looked at the density plots involving $(\phi_i,\psi_{i+1})$, $(\psi_i,\phi_{i+1})$ and $(\phi_i,\phi_{i+1})$.  For those 
von Mises distribution evaluate for each pair of amino acids, and the correspondent directional gradient, where added
 to our propensity library.
