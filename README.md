#overview

PHASTA/phasta-next/master is a superset of PHASTA/phasta/master.  As such:
 * there will be periodic merges from PHASTA/phasta/master to PHASTA/phasta-next/master, and 
 * tests that pass in PHASTA/phasta are expected to pass in PHASTA/phasta-next/master.

#development

Minor modifications to PHASTA (phSolver, phastaIO, etc...) should be made in PHASTA/phasta.  

New unpublished developments (i.e. new features) should be added in a branch off PHASTA/phasta-next/master.  To create a new branch `foo` and push it to github run the following commands:

```
git clone git@github.com:PHASTA/phasta-next.git  # clone the repo
git checkout -b foo  # create a branch of master
git push origin foo  # push the branch to github
```

Once the compilation of the branch is clean (no new warnings on the systems+compilers combinations we care about) and can pass the regression tests, then the feature branch should be merged into PHASTA/next/master.  Creating a tag marking the addition of the new feature will help later merges from PHASTA/next/master to PHASTA/phasta/master .

#build and test

    wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz .
    tar xzf phastaChefTests.tar.gz # use for CASES path below
    
Note, the following disables the SVLS and PETSC solvers and relies on LESLIB for the incompressible solver and the native compressible solver.

    cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPHASTA_INCOMPRESSIBLE=ON \
    -DPHASTA_COMPRESSIBLE=ON \
    -DPHASTA_USE_LESLIB=ON \
    -DLESLIB=/path/to/libles.a \
    -DPHASTA_USE_SVLS=OFF \
    -DPHASTA_USE_PETSC=OFF \    
    -DPHASTA_TESTING=ON \
    -DCASES=/path/to/phastaCases/ \
    ..

    make

    ctest
