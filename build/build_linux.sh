#!/bin/sh

C_COMPILER=gcc
CXX_COMPILER=g++
LINKER=ld

# Default options are
# -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_GPU=OFF -DUSE_GTEST=ON

CPU_OPTIONS=
MIC_OPTIONS=

CPU_VENDOR=$(grep -m 1 'vendor_id' /proc/cpuinfo | sed 's/^.*: //')
NUM_CORES=$(grep processor /proc/cpuinfo | wc -l)
CPU_FLAGS=$(grep flags /proc/cpuinfo | uniq)

command_exists ()
{
    type $1 > /dev/null 2>&1;
}

if command_exists icc && command_exists icpc ; then
    C_COMPILER=icc
    CXX_COMPILER=icpc
    LINKER=icpc
fi

if [ "$CPU_VENDOR" = "GenuineIntel" ] && [ "$CPU_FLAGS" = *"avx"* ] && [ "$CPU_OPTIONS" != *"-DUSE_AVX="* ]; then
    CPU_OPTIONS="-DUSE_AVX=ON $CPU_OPTIONS"
fi

buildCPU=false
buildMIC=false
clean=false
if [ $# -eq 0 ]; then
    buildCPU=true
fi

script=$0
if [ $# -ge 1 ]; then
    case $1 in
        MIC)
            buildMIC=true
        ;;
        CPU)
            buildCPU=true
        ;;
        all)
            buildCPU=true
            buildMIC=true
        ;;
        clean)
            clean=true
        ;;
        *)
            echo "[ERROR]: Incorrect option $1"
            echo "Usage: $0 [CPU|MIC|all] [ON|OFF] [ON|OFF]"
            exit 1
        ;;
    esac
fi
if [ $# -ge 2 ] ; then
    CPU_OPTIONS="$CPU_OPTIONS -DUSE_MPI=$2"
fi
if [ $# -ge 3 ]; then
    CPU_OPTIONS="$CPU_OPTIONS -DUSE_OPENMP=$3"
fi

if $buildCPU ; then
    BUILD_DIR="unix_makefiles"
    if [ ! -d $BUILD_DIR ]; then
        mkdir -p $BUILD_DIR
    fi
    cd $BUILD_DIR

    CXX=$CXX_COMPILER CC=$C_COMPILER LD=$LINKER cmake -G "Unix Makefiles" $CPU_OPTIONS ../..
    make -j $NUM_CORES -k 2> /dev/null
    if [ $? -ne 0 ]; then
        make
    fi
    cd ..
fi

if $buildMIC ; then
    BUILD_DIR="unix_makefiles_mic"
    if [ ! -d $BUILD_DIR ]; then
        mkdir -p $BUILD_DIR
    fi
    cd $BUILD_DIR

    CXX=$CXX_COMPILER CC=$C_COMPILER LD=$LINKER cmake -G "Unix Makefiles" -DUSE_MIC=ON $MIC_OPTIONS ../..
    make -j $NUM_CORES -k 2> /dev/null
    if [ $? -ne 0 ]; then
        make
    fi
    cd ..
fi

if $clean ; then
    BUILD_DIR="unix_makefiles"
    if [ -d $BUILD_DIR ]; then
        rm -rf $BUILD_DIR
    fi

    BUILD_DIR="unix_makefiles_mic"
    if [ -d $BUILD_DIR ]; then
        rm -rf $BUILD_DIR
    fi
fi
