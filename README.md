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
* *Chaos.phi(x = float, y = float, z = float)*: Returns the value of the logaritmic potential.
$$\Phi = \frac{v_0^2}{2}\log|\frac{x^2+y^2+z^2}{r_c^2}|$$
* *Chaos.a(x = float, y = float, z = float)*: Returns the value of acceleration derived from the logaritmic potentioal.
$$\textbf{a} = -\nabla\Phi$$
* *Chaos.verlet(x0 = float, y0 = float, z0 = float, k = int)*: Input of the initial position and returns the integrated trajectories in the phase space. It also returns energy and time, and the Poincare section of $x$ and $\dot{x}$. The integration method is based on Verlet integration:

$$\textbf{x}_{i+1} = \textbf{x}_i + \textbf{v}_i*dt + \frac{1}{2}\textbf{a}(\textbf{x}_i)*dt^2$$

$$\textbf{v}_{i+1} = \textbf{v}_i + \textbf{a}(x_i)*dt$$

This method returns:
  * *pos*(list): This is a list that saves all data points of the three coordinates of space. [*Ex:* pos[0] = list of x points]
  * *vel*(list): A list of the velocity in each coordinate. [*Ex*: vel[2] = list of $v_z$ points]
  * *energy* (list): List of energy for each point in \{ $\textbf{x},\textbf{v}$ \}. It also returns time and $L_z$ for the values of the first index 1 and 2 respectively.
  * *poin* (list): It is a list of the dots for the Poincare's section in y = 0, z > 0. [*Ex*: poin[0]: positions in the Poincare's section. poin[1]: velocities in the Poincare's section]
