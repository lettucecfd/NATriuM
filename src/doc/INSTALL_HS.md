
If you have trouble during the installation (which is not unlikely), contact the developers via the Google group natrium-lbm
or via email: kraemer.research@gmail.com or wilde.aerospace@gmail.com .

# Set enviromental variables
Go to the desired install folder and set environment:
```
export NATRIUM_BASE_DIR=$(pwd)
```

Environment
```
export BOOST_ROOT=$NATRIUM_BASE_DIR/libs/boost
export TRILINOS_DIR=$NATRIUM_BASE_DIR/libs/trilinos
export P4EST_DIR=$NATRIUM_BASE_DIR/libs/p4est
export DEAL_II_DIR=$NATRIUM_BASE_DIR/libs/deal.II
export NATRIUM_DIR=$NATRIUM_BASE_DIR/NATriuM
mkdir $NATRIUM_BASE_DIR/output
export NATRIUM_HOME=$NATRIUM_BASE_DIR/output
```

Write your environmental variables into a file "natriumrc" to reload them later:
```
cat > $NATRIUM_BASE_DIR/natriumrc <<EOF
export NATRIUM_BASE_DIR=$NATRIUM_BASE_DIR
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

If you have to interrupt your installation, make sure to reload the environment variables:
```
source <your natrium base dir>/natriumrc
```

# Install Required Resources

1. For boost b2: C++11 compiler `cxx-compiler`
2. For p4est: fortran77 compiler with compatible glibc `fortran-compiler`, `libgfortran5`
3. For trilinos: latest `cmake` (>=3.23), `openmpi`, `libhwloc`, `libevent`, `blas`, `liblapack`
4. For dealII: `zlib` (and `gsl` and `lapack` for Cluster) (and `cxx-compiler=1.5.2`)

Install Anaconda, if not already installed
```
cd $NATRIUM_BASE_DIR
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
bash Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
```

Update Conda
```
conda update conda
```
Create a dedicated environment
```
conda create -n "natrium"
```
Install required packages
```
conda activate natrium
conda install -c conda-forge cxx-compiler cmake libgfortran5 fortran-compiler openmpi libhwloc libevent blas liblapack zlib gsl lapack
```
Update all packages
```
conda update --all
```

**Note: On a server, you may need to specify the version of gfortran to 11.3**

Check version of cmake and, if below 3.23, install directly:
```
cd $NATRIUM_BASE_DIR
wget https://github.com/Kitware/CMake/releases/download/v3.25.3/cmake-3.25.3.tar.gz
tar -xf cmake-3.25.3.tar.gz
cd cmake-3.25.3
./bootstrap --prefix=$NATRIUM_BASE_DIR
make
make install
```

# Install resources

Alternatively, via apt-get or synaptic

### boost
from https://www.boost.org/ **Go with boost 1.76.0, not 1.82.0!**  

1. Download boost tar-file from www.boost.org.
```
cd $NATRIUM_BASE_DIR
mkdir .boost
cd .boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
```
2. Extract file and go to folder
```
tar -xf boost_1_76_0.tar.gz
cd boost_1_76_0
```
3. Install using the bootstrap shell:
```
./bootstrap.sh --prefix=$BOOST_ROOT --with-libraries=filesystem,program_options,graph,graph_parallel,iostreams,serialization,system,test,timer,thread
./b2
./b2 install
```

### p4est  
1. download **version 2.2** tarball from p4est homepage (https://www.p4est.org/; **no need to untar**) Later version had conflict with `cpp too many files`.
```
cd $NATRIUM_BASE_DIR
mkdir .p4est
cd .p4est
wget https://github.com/p4est/p4est.github.io/raw/master/release/p4est-2.2.tar.gz
```
2. get setup script from deal.II homepage
	   cf. documentation on installing deal.II with p4est
```
wget https://www.dealii.org/current/external-libs/p4est-setup.sh
```
3. Set C and C++ compilers (somehow the configuration script does not detect the right compilers, otherwise)
```
export CC=mpicc && export CXX=mpicxx
```
4. Make shell file executable and execute Setup 
```
chmod u+x p4est-setup.sh
./p4est-setup.sh p4est-2.2.tar.gz $P4EST_DIR
```

### Trilinos

Download and extraxt from https://github.com/trilinos/Trilinos/releases/tag/trilinos-release-13-0-1
```
cd $NATRIUM_BASE_DIR
wget https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-13-0-1.tar.gz
tar -xf trilinos-release-13-0-1.tar.gz
```
Then install trilinos.
```
mkdir build_trilinos
cd build_trilinos
cmake -D Trilinos_ENABLE_Sacado=ON \
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

### deal.ii  
1. download and untar tarball from deal.ii homepage https://github.com/dealii/dealii/releases
	(rename directory if it has the name of your target directory)
	*Did not work with latest Trilinos, so I downgraded to Trilinos 13.0.1*
	**deal.ii version 9.3.3 works**

**deal II is compiled without zlib, but runs a test compilation on mpicxx and mpicc, which fails in Siegen. You may need ot manually install/link it.**
In this case, search for the conda location and add this to options, e.g., `-D ZLIB_LIBRARY=~/miniconda3/pkgs/zlib-1.2.13-hd590300_5/lib/libz.so -D ZLIB_INCLUDE_DIR=~/miniconda3/pkgs/zlib-1.2.13-hd590300_5/include`.
```
cd $NATRIUM_BASE_DIR
wget https://github.com/dealii/dealii/releases/download/v9.3.3/dealii-9.3.3.tar.gz
tar -xf dealii-9.3.3.tar.gz
```

2. Setup installation (replace version of dealii!)
```
mkdir build_deal
cd build_deal
cmake -D CMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
-D DEAL_II_WITH_PETSC=OFF \
-D DEAL_II_WITH_TRILINOS=ON \
-D DEAL_II_WITH_MPI=ON \
-D DEAL_II_COMPONENT_PARAMETER_GUI=OFF \
-D DEAL_II_WITH_BOOST=ON \
-D DEAL_II_ALLOW_BUNDLED=OFF \
-D BOOST_DIR=$BOOST_ROOT \
-D DEAL_II_WITH_THREADS=OFF \
-D BOOST_ROOT=$BOOST_ROOT \
-D P4EST_DIR=$P4EST_DIR \
-D DEAL_II_WITH_P4EST=ON \
-D DEAL_II_FORCE_BUNDLED_UMFPACK=ON \
-D DEAL_II_FORCE_BUNDLED_MUPARSER=ON \
-D DEAL_II_WITH_ZLIB=OFF \
../dealii-*/
```

3. Install (-j 8 enables parallel compilation on  processors; otherwise installation will take hours)
```
make -j 8 install
```
 
### Check libraries

Your $NATRIUM_BASE_DIR/libs folder should now contain Boost, Dealii, p4est, and Trilinos libraries! 

# Get and compile NATriuM Code

1. Get repo and compile (compressible branch is more recent)

```
cd $NATRIUM_BASE_DIR
git clone https://github.com/lettucecfd/NATriuM.git
cd $NATRIUM_DIR
git checkout compressible
mkdir bin_debug
cd bin_debug
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug ../src/ -B.
make -j8
```

2. Load project via IDE (e.g., CLion (recommended, view as makefile project) or Eclipse for Embedded C/C++ Developers) 
3. Rebuild for for `bin_release` instead of `bin_debug` and `-DCMAKE_BUILD_TYPE=Release` instead of `-DCMAKE_BUILD_TYPE=Debug` to get a fast version of the program

```
cd $NATRIUM_DIR
mkdir bin_release
cd bin_release
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ../src/ -B.
make -j8
```

4. Check that the CMakeCache.txt contains `CMAKE_BUILD_TYPE:STRING=RELEASE`! Otherwise the program will be really slow

# Test the Code

## Run unit tests:

```
cd $NATRIUM_DIR/bin_release
./test/NATriuM_UnitTest_exe
```

## Run integration tests (takes a few minutes)
```
cd $NATRIUM_DIR/bin_release
./test/NATriuM_Test
```

The results will be written to natrium.html (in the bin directory)

# Getting started

To get started, take a look at the Mainpage of the technical Documentation (doc/html/index.html) and navigate to the Examples section.
        
