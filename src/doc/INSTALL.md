
If you have trouble during the installation (which is not unlikely), contact the developers via the Google group natrium-lbm
or via email: kraemer.research@gmail.com or wilde.aerospace@gmail.com .

The following configuration works fine:
- Boost 1.76
- Trilinos 13.0.1
- P4est 2.2
- Deal 9.3.3


# INSTALL REQUIRED RESSOURCES

**If you already installed openLB, you already have g++, openmpi-bin, openmpi-doc, libopenmpi-dev**, so just do
```
sudo apt-get install libblas-dev liblapack-dev
sudo apt-get install gfortran
conda install -c conda-forge openmpi cmake blas zlib libevent libblas liblapack
```

**Conda-forge's cxx-compiler automatically installs the correct c++ compiler (gcc for linux)**
```
conda install -c conda-forge libblas-dev liblapack-dev gfortran openmpi cxx-compiler cmake gfortran blas zlib libevent libblas liblapack
```

**From scratch, conda only**
For boost b2: C++11 compiler
```
conda install -c conda-forge cxx-compiler
```

For p4est: fortran77 compiler
```
conda install -c conda-forge fortran-compiler
```

# Set enviromental variables
## Go to the desired install folder and set environment:
```
mkdir .natrium
cd .natrium
export NATRIUM_BASE_DIR=$(pwd)
```

Environment
```
export BOOST_ROOT=$NATRIUM_BASE_DIR/libs/boost
export TRILINOS_DIR=$NATRIUM_BASE_DIR/libs/trilinos
export P4EST_DIR=$NATRIUM_BASE_DIR/libs/p4est
export DEAL_II_DIR=$NATRIUM_BASE_DIR/libs/deal.II
export NATRIUM_DIR=$NATRIUM_BASE_DIR/NATriuM
export NATRIUM_HOME=$NATRIUM_BASE_DIR/output
```

Write your environmental variables into a file "natriumrc" to reload them later:
```
cat > $NATRIUM_BASE_DIR/natriumrc <<EOF
export BOOST_ROOT=$BOOST_ROOT
export TRILINOS_DIR=$TRILINOS_DIR
export P4EST_DIR=$P4EST_DIR
export DEAL_II_DIR=$DEAL_II_DIR
export NATRIUM_DIR=$NATRIUM_DIR
export NATRIUM_HOME=$NATRIUM_HOME
export LD_LIBRARY_PATH=$BOOST_ROOT/lib:$LD_LIBRARY_PATH
export INCLUDE_PATH=$BOOST_ROOT:$INCLUDE_PATH
EOF
```

# Install resources
## via apt-get

## OR synaptic

## OR Manually)
1. Install Boost from https://www.boost.org/ **Go with boost 1.76.0, not 1.82.0!**	
    1.1 Download boost tar-file from www.boost.org.
    1.2 Extract file
    1.3 Go to folder
    1.4 Execute:
	```
./bootstrap.sh --prefix=$BOOST_ROOT --with-libraries=filesystem,program_options,graph,graph_parallel,iostreams,serialization,system,test,timer,thread
./b2
./b2 install
	```
2. p4est
    2.1 download tarball from p4est homepage (https://www.p4est.org/ here: version 2.2; **no need to untar**)
    2.2 get setup script from deal.II homepage (cf. documentation on installing deal.II with p4est)
    2.3 Set C and C++ compilers
	```
export CC=mpicc && export CXX=mpicxx
	```
	(somehow the configuration script does not detect the right compilers, otherwise)
    2.4 Execute Setup
    	```
./p4est-setup.sh <p4est tarball> $P4EST_DIR
	```
3. Trilinos
```
git clone https://github.com/trilinos/Trilinos.git

mkdir build_trilinos
cd build_trilinos
cmake 	-D Trilinos_ENABLE_Sacado=ON \
-D Trilinos_ENABLE_Stratimikos=ON \
-D Trilinos_ENABLE_MueLu=ON \
-D CMAKE_BUILD_TYPE=RELEASE \
-D CMAKE_CXX_FLAGS="-g -O3" \
-D CMAKE_C_FLAGS="-g -O3" \
-D CMAKE_FORTRAN_FLAGS="-g -O5" \
-D Trilinos_EXTRA_LINK_FLAGS="-lgfortran" \
-D CMAKE_VERBOSE_MAKEFILE=FALSE \
-D Trilinos_VERBOSE_CONFIGURE=FALSE \
-D TPL_ENABLE_MPI=ON \
-D TPL_ENABLE_Pthread=OFF \
-D BUILD_SHARED_LIBS=ON \
-D CMAKE_INSTALL_PREFIX:PATH=$TRILINOS_DIR \
../Trilinos*/

make -j8
make install
```

4. deal.ii
 	4.1 download and untar tarball from deal.ii homepage
		(rename directory if it has the name of your target directory)
        	(to get newest dealii version: git clone git://git@github.org/dealii/dealii.git dealii-git)
	4.2 Execute
	```
mkdir build_deal; 
cd build_deal; 
cmake -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR -DDEAL_II_WITH_PETSC=OFF -DDEAL_II_WITH_TRILINOS=ON -DDEAL_II_WITH_MPI=ON -DDEAL_II_COMPONENT_PARAMETER_GUI=OFF -DDEAL_II_WITH_BOOST=ON -DDEAL_II_ALLOW_BUNDLED=OFF -DBOOST_DIR=$BOOST_ROOT -DDEAL_II_WITH_THREADS=OFF -DBOOST_ROOT=$BOOST_ROOT -DP4EST_DIR=$P4EST_DIR -DDEAL_II_WITH_P4EST=ON -DDEAL_II_FORCE_BUNDLED_UMFPACK=ON -DDEAL_II_FORCE_BUNDLED_MUPARSER=ON -DDEAL_II_WITH_ZLIB=OFF ../dealii-9.3.3
	```
        4.3 Execute
	```
make -j 8 install
	```
	 (-j 8 enables parallel compilation on  processors; otherwise installation will take hours)
 
5. Check: Your $NATRIUM_BASE_DIR/libs folder should now contain Boost, Dealii, p4est, and Trilinos libraries! 

================================================================================================
   GET AND COMPILE NATRIUM CODE
================================================================================================
	0) cd $NATRIUM_BASE_DIR
	1) clone https://github.com/lettucecfd/NATriuM.git
	2) cd NATriuM; mkdir bin_debug; cd bin_debug
	3) cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug ../src/ -B.
	3) make -j8
	5) Load project via IDE (e.g., Eclipse or CLion (recommended))
	6) Repeat steps 1-5) for "bin_release" instead of "bin_debug" and "-DCMAKE_BUILD_TYPE=Release" instead of "-DCMAKE_BUILD_TYPE=Debug" to get a fast version of the program
	7) Check that the CMakeCache.txt: CMAKE_BUILD_TYPE:STRING= must be set to RELEASE! Otherwise the program will be really slow
CMAKE_BUILD_TYPE:STRING=Release



        
================================================================================================
   TEST THE CODE 
================================================================================================

Run unit tests:
   - go to NATriuM's bin directory
   - type ./test/NATriuM_UnitTest_exe

Run integration tests (takes a few minutes)
   - go to NATriuM's bin directory
   - type ./test/NATriuM_Test
   - The results will be written to natrium.html (in the bin directory)


================================================================================================
  GETTING STARTED
================================================================================================
To get started, take a look at the Mainpage of the technical Documentation (doc/html/index.html)
and navigate to the Examples section.
        
