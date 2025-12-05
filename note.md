# Technical Note

In this note, we will go over the approach used in the code to calculate the radial distribution function and the mean squared displacement, as well as the optimization strategy that was chosen. Let us first recall the definition and purpose of the RDF and the MSD.

## Physical meaning and utility of MSD and RDF

### Mean squared displacement

The MSD tracks how far, on average, a particle drifts from its starting position over a given time lag $\tau$. Plotting the relation between the time lag and the mean squared displacement can for example help in finding the diffusivity constant $D$.

The mean squared displacement of a single particle $i$, averaged over time is given by:

$$\Delta_i^2(m)=\frac{1}{N_t-m}\sum_{k=0}^{N_t-m-1}\Vert\vec{r_i}(k+m)-\vec{r_i}(k)\Vert^2$$

with $\vec{r_i}(k)$ the position vector of the particle $i$ at time $k\Delta t$ (assuming a constant time step), $N_t$ the total number of time steps in the trajectory and $m$ the index such that the time lag $\tau=m\Delta t$ ([source](https://www.neutron-sciences.org/articles/sfn/pdf/2011/01/sfn201112010.pdf)).

### Radial distribution function

Given a certain atom, the RDF gives the probability of finding another atom at a distance $r$. It represents the local structure and packing of the fluid.

For distances close to zero, the probabilty of finding another atom is usually zero because atoms cannot overlap. The first peak of the function represents the equilibrium bond length, while the valleys are the regions where it is energetically unfavorable for atoms to sit. At long distances, the atom at the origin no longer influences positions so the RDF becomes constant.

The RDF can be used to identify the state of the matter. For instance, sharp and distinct peaks that persist over long distances typically come from solids, while a strong first peak with weaker following peaks that quickly come to oscillate around a constant value typically come from liquids. Gas on the other hand would show a flatter line. These observations obviously strongly depend on the temperature and other factors.

The formula of the radial distribution function is:

$$g(r)=\frac{n(r)}{\rho \times 4\pi r^2\times\Delta r}$$

with $n(r)$ the mean number of atoms in a shell of width $\Delta r$ at distance $r$ from some central atom, $\rho=N_{\mathrm{atoms}}\,/\,V$ the mean density ([source](https://people.bath.ac.uk/chsscp/teach/md.bho/Theory/rdf.html)).

## Algorithmic approach for the MSD and RDF

### Assumptions

Below is a list of the main assumptions made regarding the code:

- Origin of the box at $(0,\,0,\,0)$
- Constant simulation box dimensions throughout a simulation
- Only cubic boxes (because the box dimensions for the simulations are not given as vectors)
- Constant number of atoms throughout a simulation
- `.xyz` files don't include atom forces or magnetic moments
- Unwrapped coordinates
- **Mean Squared Displacement**:
  - assumed a constant time step between the frames
  - computes a time **and** ensemble averaged MSD
- **Radial Distribution Function**:
  - cutoff distance shouldn't exceed half of the smallest side of the box
  - computes a time **and** ensemble averaged RDF

### Handling of the periodic boundary conditions

The particle coordinates are unwrapped, which means they can drift outside of the simulation box boundaries. This is advantageous to compute the MSD, but for the radial distribution function it is important to wrap the coordinates and apply periodic boundary conditions. To achieve this, let's look at how to implement the minimum image convention.

Say we have an atom $i$ with a position vector $\vec{r_i}$. Let's consider another atom $j\neq i$ with a position vector $\vec{r_j}$. Let's define $\vec{r_{ij}}=\vec{r_j}-\vec{r_i}$. Even if it means shifting all vectors, we can assume atom $i$ is sitting right in the middle of the simulation box, thus $\vec{r_i}=(L_x/2,\,L_y/2,\,L_z/2)^T$ with $L_x,\,L_y$ and $L_z$ the length of the sides of the box (assumed cubic). Even if atom $j$ drifted outside the box, the PBC would "project" its image inside the box, thus:

$$|\Delta x_{ij}|\leq\frac{L_x}{2},\,|\Delta y_{ij}|\leq\frac{L_y}{2},\,|\Delta z_{ij}|\leq\frac{L_z}{2}$$

This means that we can always find an "image" of atom $j$ such that the statement above is true. To do this, let's define the matrix $\mathbf{L}$ whose columns are the lattice vectors. In our case:

$$\mathbf{L}=\begin{bmatrix}
L_x & 0 & 0 \\
0 & L_y & 0 \\
0 & 0 & L_z
\end{bmatrix}$$

By definition, the determinant of $\mathbf{L}$ should be the volume of the simulation box (non-zero), thus $\mathbf{L}$ is invertible. Hence, the columns of this matrix define a basis. We can thus express any vector $\vec{r}$ in this basis: $\vec{r}=\mathbf{L}\vec{f}$. We usually refer to $\vec{f}$ as the vector of the fractional coordinates. Let's consider the fractional coordinates of $\vec{r_{ij}}$:

$$\vec{f_{ij}}=\mathbf{L}^{-1}\vec{r_{ij}}$$

Applying what we've seen before, we can translate the previous vector to make it so:

$$|f_{ij,\,x}|\leq\frac{1}{2},\,|f_{ij,\,y}|\leq\frac{1}{2},\,|f_{ij,\,z}|\leq\frac{1}{2}$$

To achieve this, we can simply remove the rounded coordinates from the original vector:

$$\vec{g_{ij}}=\vec{f_{ij}}-\vec{\mathrm{round}}\,(\vec{f_{ij}})$$

The new coordinates thus obtained with this method are:

$$\vec{r_{ij,\,\mathrm{pbc}}}=\mathbf{L}\vec{g_{ij}}=\vec{r_{ij}}-\mathbf{L}\vec{\mathrm{round}}\,(\mathbf{L}^{-1}\vec{r_{ij}})$$

## Optimization strategy

> `numba` is a library that compiles python code *just-in-time*, and makes use of parallelization.
The first time you run a function that uses `numba`, it is compiled once and stored in cache for as long as you keep your python kernel running.
`numba` does not require any preconfiguration, meaning you can use it right away without having to specify any parameters (like the number of cpu...).
>
> While the use of *just-in-time* compilation allows the code to run much faster, it can also be less memory-intensive than `Numpy`.
When dealing with very heavy simulations, storing and dealing with all the coordinates in a single `Numpy` array can overload the RAM.
Although this task doesn't require a lot of memory considering the given MD trajectory, I think the use of `numba` is still very good.
>
> The functions inside the modules `fastmsd` and `fastrdf` use `numba`.
>
> [Introduction of the documentation]

Below is a table with the mean execution times of the computation of the MSD and the RDF obtained using the `perf_counter` function of the `time` module. Note that the first column gives the first execution time (which includes the compilation), while the second gives the true execution time once `numba` has compiled the function.

|         | Execution + compilation time | Execution time |
|:-------:|:----------------------------:|:--------------:|
| FastMSD |         1.44784 [s]          |   0.22044 [s]  |
| FastRDF |          5.4116 [s]          |   0.86802 [s]  |
