#!/bin/sh

C_COMPILER=gcc
CXX_COMPILER=g++
LINKER=ld

# Default options are
# -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_GPU=OFF -DUSE_GTEST=ON -DUSE_AVX=OFF

CPU_OPTIONS="-DCMAKE_CXX_FLAGS=-DGTEST_USE_OWN_TR1_TUPLE=1 -DCMAKE_MACOSX_RPATH=1"
MIC_OPTIONS=

CPU_VENDOR=$(sysctl -n machdep.cpu.vendor)
NUM_CORES=$(sysctl -n hw.physicalcpu)

function command_exists ()
{
    type $1 > /dev/null 2>&1;
}

if command_exists icc && command_exists icpc ; then
    C_COMPILER=icc
    CXX_COMPILER=icpc
    LINKER=icpc
fi

if [[ $CPU_VENDOR == "GenuineIntel" ]] && [[ $CPU_OPTIONS != "*-DUSE_AVX=*" ]]; then
    CPU_OPTIONS="-DUSE_AVX=ON $CPU_OPTIONS"
fi

buildCPU=false
buildMIC=false
clean=false
if [ $# -eq 0 ]; then
    buildCPU=true
fi

script=$0
while [ $# -ne 0 ]
do
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
            echo "Usage: $0 [CPU|MIC|all]"
            exit 1
        ;;
    esac
    shift
done

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
