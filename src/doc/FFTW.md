## Install

Load environment variables
```
source ~/natrium/natriumrc

```

Download and extract tarball
```
cd $NATRIUM_INSTALLATION_DIR/
wget http://fftw.org/fftw-3.3.10.tar.gz
tar -xf fftw-*.tar.gz

```

Install
```
cd $NATRIUM_INSTALLATION_DIR/fftw-*/
./configure --prefix=$NATRIUM_BASE_DIR/libs/fftw --enable-threads --enable-openmp --enable-mpi
make
make install

```
