Introduction
============

Overview
--------

PhaseTreeTask is a codebase that aims to fulfill the assignment I was given as part of my application to PhaseTree.
For a detailed rundown on how to use the code, please look at the :doc:`examples` section.

Modules
-------

The codebase consists of three modules:

- :doc:`trajectory`: contains the ``Trajectory`` class which may be instantiated with a single `.xyz` file path.
- :doc:`fastmsd`: contains a function to compute the **Mean Squared Displacement** of a trajectory.
- :doc:`fastrdf`: contains a function to compute the **Radial Distribution Function** of a trajectory.

Numba compiler
--------------

`numba` is a library that compiles python code *just-in-time*, and makes use of parallelization.
The first time you run a function that uses `numba`, it is compiled once and stored in cache for as long as you keep your python kernel running.
`numba` does not require any preconfiguration, meaning you can use it right away without having to specify any parameters (like the number of cpu...).

While the use of *just-in-time* compilation allows the code to run much faster, it can also be less memory-intensive than `Numpy`.
When dealing with very heavy simulations, storing and dealing with all the coordinates in a single `Numpy` array can overload the RAM.
Although this task doesn't require a lot of memory considering the given MD trajectory, I think the use of `numba` is still very good.

The functions inside the modules :doc:`fastmsd` and :doc:`fastrdf` use `numba`.

Requirements and versions
-------------------------

The list of all packages used in the environment can be found in a file called ``environment.yml`` under the ``config`` folder. Here are the main requirements:

- Python: ``3.13.2``
- Numpy: ``2.2.2``
- Matplotlib: ``3.10.0``
- Numba: ``0.61.2``
