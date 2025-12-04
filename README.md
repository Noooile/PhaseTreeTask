# Molecular Dynamics Trajectory Analysis & Optimization Assignment

This codebase aims to fulfill the assignment I was given as part of my application to PhaseTree.

## Assumptions

- Constant simulation box dimensions throughout a simulation
- Constant number of atoms throughout a simulation
- `.xyz` files don't include atom forces or magnetic moments
- Cubic boxes, because the box dimensions for the simulations are not given as vectors
- **Mean Squared Displacement**:
  - assumed a constant time step between the configurations
  - computes a time **and** ensemble averaged MSD
- **Radial Distribution Function**:
  - cutoff distance shouldn't exceed half of the smallest side of the box

## Content

The Python code can be found under the `src` folder. The `res` folder contains the data file given for the task. The Python environment used to write the code can be cloned using the YAML file under the `config` folder.

## References

[MSD formula](https://www.neutron-sciences.org/articles/sfn/pdf/2011/01/sfn201112010.pdf)

[RDF formula](https://people.bath.ac.uk/chsscp/teach/md.bho/Theory/rdf.html)
