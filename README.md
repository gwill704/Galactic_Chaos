# Galactic_Chaos

The objective of this repository is to integrate the trajectories of stars in galaxies, as we can study by variating two parameters b and c the chaos generated in the solutions. 

The programming language used is **Python3.10**. 
We include a virtual environment and a *requirements.txt* to install all the external packages. We used:
- Matplotlib
- Numpy
- Random

## CHAOS.py file
This is the main file of the repository. In it we find the definition of the class **Chaos** in which we find all the necesary tools for the study.

**Input**
* *Initial velocities*: vx0, vy0, vz0 in all three cartesian axes.
* *Time-step of integration*: dt
* *Asymetry parameters*: b and c. We recommend for physical cases to maintain the criteria c < b and both between 0 and 1. The parameter b introduces an asymetry in the y coordinate and the c in the z component. Both act on the potential $\Phi$.
* *Reference length*: rc. It will act on the $\Phi$ aswell. For bound-states we want rc at least bigger than the initial module $|\textbf{r}| = \sqrt{x^2+y^2+z^2}$

**Methods**
* *Chaos.v(vx = float, vy = float, vz = float)*: Returns the module of the velocity.
* *Chaos.phi(x = float, y = float, z = float)*: 
