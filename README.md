# 2016_CIG_SPECFEM3D

Getting started
Connect to cluster: 
  ssh -X id@mcmillan2 

Download the SPECFEM3D_Cartesian software package: 
  git clone --recursive --branch devel https://github.com/geodynamics/specfem3d.git

Load appropriate fortran and MPI compilers: 
  module load intel/16.0/64/16.0.1.150 
  module load openmpi/intel-16.0/1.10.2/64 
  
Load python distribution:  
  module load anaconda3/2.5.0 
Configure SPECFEM3D_Cartesian for your system from the SPECFEM3D root directory: 
  cd specfem3d 
  ./configure FC=ifort MPIFC=mpif90
  
Compile all the source code  
  make all 
Check generated executable files: 
  ls ./bin 
