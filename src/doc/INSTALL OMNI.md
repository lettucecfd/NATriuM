
If you have trouble during the installation (which is not unlikely), contact the developers via the Google group natrium-lbm
or via email: kraemer.research@gmail.com or wilde.aerospace@gmail.com .

**Full description in [INSTALL HS](https://github.com/PhiSpel/NATriuM/blob/patch-1/src/doc/INSTALL_HS.md)**

# Set enviromental variables
Go to the desired install folder and set:

```
export NATRIUM_BASE_DIR=$(pwd)
export BOOST_ROOT=$NATRIUM_BASE_DIR/libs/boost
export TRILINOS_DIR=$NATRIUM_BASE_DIR/libs/trilinos
export P4EST_DIR=$NATRIUM_BASE_DIR/libs/p4est
export DEAL_II_DIR=$NATRIUM_BASE_DIR/libs/deal.II
export NATRIUM_DIR=$NATRIUM_BASE_DIR/NATriuM
mkdir $NATRIUM_BASE_DIR/output
export NATRIUM_HOME=$NATRIUM_BASE_DIR/output
mkdir $NATRIUM_BASE_DIR/installation
export NATRIUM_INSTALLATION_DIR=$NATRIUM_BASE_DIR/installation
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
export NATRIUM_INSTALLATION_DIR=$NATRIUM_INSTALLATION_DIR
EOF
```

If you have to interrupt your installation, make sure to reload the environment variables:
```
source <your natrium base dir>/natriumrc
```

# Install Conda Resources

Install Anaconda, if not already installed
```
cd $NATRIUM_INSTALLATION_DIR
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
bash Miniconda3-py310_23.3.1-0-Linux-x86_64.sh
```

Update Conda
```
conda update conda
conda create -n "natrium"

```
```
conda activate natrium
conda install -c conda-forge cmake blas libgfortran5 liblapack mkl

```
*In my case, zlib is also included in `conda list`. DealII is installed without, but the test program required it for testing the cxx compiler (mpicxx).*
```
conda update --all
```

# Install libraries

### Boost
Boost is available and may be loaded with `module load boost`, but this has not been tested.

Instead:
```
cd $NATRIUM_INSTALLATION_DIR
mkdir .boost
cd .boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
tar -xf boost_1_76_0.tar.gz
cd boost_1_76_0
./bootstrap.sh --prefix=$BOOST_ROOT --with-libraries=filesystem,program_options,graph,graph_parallel,iostreams,serialization,system,test,timer,thread
./b2
./b2 install
```

### p4est  

```
cd $NATRIUM_INSTALLATION_DIR
mkdir .p4est
cd .p4est
wget https://github.com/p4est/p4est.github.io/raw/master/release/p4est-2.2.tar.gz
wget https://www.dealii.org/current/external-libs/p4est-setup.sh
export CC=mpicc && export CXX=mpicxx
chmod u+x p4est-setup.sh
./p4est-setup.sh p4est-2.2.tar.gz $P4EST_DIR
```

### Trilinos

It's not as simple as loading the module `module load trilinos`, as apparently the installation cannot be found, then.

So, do it the long way:
```
cd $NATRIUM_INSTALLATION_DIR
wget https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-cnars-20191023.tar.gz
tar -xf trilinos-cnars-20191023.tar.gz
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

[//]: # (**deal II is compiled without zlib, but runs a test compilation on mpicxx and mpicc, which fails in Siegen. You may need ot manually install/link it.**
In this case, search for the conda location and add this to options, e.g., `-D ZLIB_LIBRARY=~/miniconda3/pkgs/zlib-1.2.13-hd590300_5/lib/libz.so -D ZLIB_INCLUDE_DIR=~/miniconda3/pkgs/zlib-1.2.13-hd590300_5/include`.)
**make sure `$BOOST_ROOT` is still set correctly**

```
cd $NATRIUM_INSTALLATION_DIR
wget https://github.com/dealii/dealii/releases/download/v9.3.3/dealii-9.3.3.tar.gz
tar -xf dealii-9.3.3.tar.gz
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
make -j 8 install
```
 
### Check libraries

Your $NATRIUM_BASE_DIR/libs folder should now contain Boost, Dealii, p4est, and Trilinos libraries! 

### Tidy up

You should be able to delete $NATRIUM_INSTALLATION_DIR

# Get and compile NATriuM Code

```
cd $NATRIUM_BASE_DIR
git clone https://github.com/lettucecfd/NATriuM.git
cd $NATRIUM_DIR
git checkout compressible
mkdir bin_debug
cd bin_debug
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug ../src/ -B.
make -j8
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
        

