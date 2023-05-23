
If you have trouble during the installation (which is not unlikely), contact the developers via the Google group natrium-lbm
or via email: kraemer.research@gmail.com or wilde.aerospace@gmail.com .

# Install Required Resources

1. For boost b2: C++11 compiler `cxx-compiler`
2. For p4est: fortran77 compiler `fortran-compiler`
3. For trilinos: latest `cmake` (>=3.23), `openmpi`, `libhwloc`, `libevent`, `blas`
4. For dealII: `zlib`

```
conda install -c conda-forge cxx-compiler fortran-compiler cmake openmpi libhwloc libevent blas zlib
```

# Set enviromental variables
Go to the desired install folder and set environment:
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

Alternatively, via apt-get or synaptic

### boost
from https://www.boost.org/ **Go with boost 1.76.0, not 1.82.0!**  

1. Download boost tar-file from www.boost.org.
2. Extract file
3. Go to folder
4. Execute:

```
./bootstrap.sh --prefix=$BOOST_ROOT --with-libraries=filesystem,program_options,graph,graph_parallel,iostreams,serialization,system,test,timer,thread
./b2
./b2 install
```

### p4est  
1. download **version 2.2** tarball from p4est homepage (https://www.p4est.org/; **no need to untar**) Later version had conflict with `cpp too many files`.
2. get setup script from deal.II homepage
	   cf. documentation on installing deal.II with p4est
3. Set C and C++ compilers (somehow the configuration script does not detect the right compilers, otherwise)

```
export CC=mpicc && export CXX=mpicxx
```

4. Execute Setup 

```
./p4est-setup.sh <p4est tarball> $P4EST_DIR`
```

### Trilinos

Download from https://github.com/trilinos/Trilinos/releases/tag/trilinos-release-13-0-1

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
2. Setup installation (replace version of dealii!)

```
mkdir build_deal
cd build_deal
cmake -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR -DDEAL_II_WITH_PETSC=OFF -DDEAL_II_WITH_TRILINOS=ON -DDEAL_II_WITH_MPI=ON -DDEAL_II_COMPONENT_PARAMETER_GUI=OFF -DDEAL_II_WITH_BOOST=ON -DDEAL_II_ALLOW_BUNDLED=OFF -DBOOST_DIR=$BOOST_ROOT -DDEAL_II_WITH_THREADS=OFF -DBOOST_ROOT=$BOOST_ROOT -DP4EST_DIR=$P4EST_DIR -DDEAL_II_WITH_P4EST=ON -DDEAL_II_FORCE_BUNDLED_UMFPACK=ON -DDEAL_II_FORCE_BUNDLED_MUPARSER=ON -DDEAL_II_WITH_ZLIB=OFF ../dealii-9.3.3
```

3. Install (-j 8 enables parallel compilation on  processors; otherwise installation will take hours)

```
make -j 8 install
```
 
### Check libraries

Your $NATRIUM_BASE_DIR/libs folder should now contain Boost, Dealii, p4est, and Trilinos libraries! 

# Get and compile NATriuM Code

1. Get repo and compile

```
cd $NATRIUM_BASE_DIR
git clone https://github.com/lettucecfd/NATriuM.git
cd NATriuM
mkdir bin_debug
cd bin_debug
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Debug ../src/ -B.
make -j8
```

2. Load project via IDE (e.g., CLion (recommended, view as makefile project) or Eclipse for Embedded C/C++ Developers) 
3. Rebuild for for `bin_release` instead of `bin_debug` and `-DCMAKE_BUILD_TYPE=Release` instead of `-DCMAKE_BUILD_TYPE=Debug` to get a fast version of the program

```
cd $NATRIUM_BASE_DIR/NATriuM
mkdir bin_release
cd bin_release
cmake -G "Eclipse CDT4 - Unix Makefiles" -DCMAKE_BUILD_TYPE=Release ../src/ -B.
make -j8
```

4. Check that the CMakeCache.txt contains `CMAKE_BUILD_TYPE:STRING=RELEASE`! Otherwise the program will be really slow

# Test the Code

## Run unit tests:

```
cd $NATRIUM_BASE_DIR/bin_release
./test/NATriuM_UnitTest_exe
```

## Run integration tests (takes a few minutes)
```
cd $NATRIUM_BASE_DIR/bin_release
./test/NATriuM_Test
```

The results will be written to natrium.html (in the bin directory)

# Getting started

To get started, take a look at the Mainpage of the technical Documentation (doc/html/index.html) and navigate to the Examples section.
        

