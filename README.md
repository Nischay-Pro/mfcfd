# Installation
The Meshfree Solver uses Python 3.x and cmake to setup the environment and compile the executable. 

## Prerequisites
Meshfree Solver requires

### Core Dependencies
* Python 3.x (for self-installer)
* HDF5 (For reading **.h5** files)

Apart from these dependencies the solver might require a number of transitive dependencies depending on the solver being compiled.

### Serial Solver Dependency
* PETSc, MPI (OpenMPI or MPICH)

## Package Outline 

* `src_mpi_serial` contains the code for the MPI version of the primal solver. 

In addition to these, we have
* `install` contains the Python script to generate the executable for main verseions of the solvers.
* `examples` contains the subdirectories for the versions of the code containing all the required files and inputs via which the executable will run.

## Building
If building on the local compute machine which has Spack, first do the following to load the dependencies:
```bash
spack load hdf5%gcc
spack load petsc%gcc
```

Then run the following command from the `install` directory:

`python3 install.py --mfcfd serial --cextra="-DPETSC_EXECUTABLE_RUNS=ON"`

After building, an executable `execname` will be generated.

## Environment Setup

Move the generated executable `execname`, created in the install directory to inside `mfcfd/examples/mpi_serial`
* Rename `case.in.example` to `case.in` 

Finally in the same directory that the excutable has been shifted to, create a subdirectory called `point` and store the required HDF5 grid file to be run in that directory.  
For example `mkdir point` such that the path is `mfcfd/examples/mpi_serial/point` and store `point.h5` in this newly created subdirectory such that its path will be `mfcfd/examples/mpi_serial/point/point.h5`


## Execution

```
mpirun -np x ./execname
# X is the number of CPU cores on which the  solver will run
```

## Additional Information

* Major parameters like the `number of iterations`, `cfl`, `venkat_limter constant` can be changes in the required `case.in` files 
* Most of the data structures and global variables are declared in the `data_structure.F90` or similarly named files
* Most of the file reading and associated memory alllocation subroutines can be found in the `point_preprocessor.F90` files.
* Most of the output related subroutines can be found in `post_processing.F90` files.
* Primary file for the solvers is the `meshfree_solver.F90` file.
* The actual computational subroutines can be found or are called in `q_lskum.F90` or similar named files and `fpi_solver.F90` files. 
* MPI/PETSc related communication calls and structures can be found in `petsc_data_structure.F90`
