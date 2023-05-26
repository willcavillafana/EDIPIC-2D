FC := mpif90
CC := mpicc
# Get the actual Fortran compiler invoked by mpif90. Need to use double $$ in awk.
BASE_FCOMP := $(shell mpif90 -show | awk '{print $$1}' | awk -F/ '{print $$NF}' | awk -F- '{print $$1}')

$(info ${BASE_FCOMP} )

# Make sure that the PETSC_DIR environment variable is set. Otherwise exit with error message.
ifeq ($(strip $(PETSC_DIR)),)
   $(error Please set the PETSC_DIR environment variable)
endif

# Git version 
GIT_HASH=$(shell git rev-parse HEAD)
# COMPILE_TIME=$(shell date +'%Y-%m-%d %H:%M:%S EST')
GIT_BRANCH=$(shell git branch | grep "^\*" | sed 's/^..//')
GIT_DATE=$(shell git log -1 --format=%cd)
# GIT_HASH=`git rev-parse HEAD`
# COMPILE_TIME=`date +'%Y-%m-%d %H:%M:%S EST'`
# GIT_BRANCH=`git branch | grep "^\*" | sed 's/^..//'`
# GIT_DATE=`git log -1 --format=%cd`
#export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\"" -DGIT_DATE="\"$(GIT_DATE)\""


ifeq ($(BASE_FCOMP),ifort)
   $(info "Using Intel compiler" )
   FFLAGS := -O2 -I${PETSC_DIR}/include #$(VERSION_FLAGS)
   CFLAGS := -O2 #$(VERSION_FLAGS)
else ifeq ($(BASE_FCOMP),gfortran)
   $(info "Using gcc/gfortran compilers" )
   FFLAGS := -O2 -ffree-line-length-none -I${PETSC_DIR}/include #$(VERSION_FLAGS)
   CFLAGS := -O2 #$(VERSION_FLAGS)
else
   FFLAGS := -O2 -I${PETSC_DIR}/include #$(VERSION_FLAGS)
   CFLAGS := -O2 #$(VERSION_FLAGS)
endif



# The only external functions explicitly called by EDIPIC2D are from MPI and PETSc.
# However, PETSc will internally call functions from HYPRE, LAPACK, and BLAS. In the
# case of the INTEL compiler, the MKL library can replace LAPACK and BLAS, although
# HYPRE includes its own copy of the LAPACK and BLAS source files.
# These library dependencies can be satisfied via shared libraries by listing them
# in the LD_LIBRARY_PATH environment variable.
# Before linking and running, make sure to set:
#    export LD_LIBRARY_PATH=${PETSC_DIR}/lib:${LD_LIBRARY_PATH}    # for bash
#    setenv LD_LIBRARY_PATH ${PETSC_DIR}/lib:${LD_LIBRARY_PATH}    # for csh

LIBS := -L${PETSC_DIR}/lib -lpetsc -lHYPRE

OBJ := WELL19937a_new.o \
       pic2d_Modules.o \
       mod_timers.o \
       mod_def_timers.o \
       pic2d_Snapshots.o \
       pic2d_SnapshotsAvg.o \
       pic2d_Prepare_FFTX_SYSY.o \
       pic2d_RandomNumberInterface.o \
       pic2d_TimeDependences.o \
       pic2d_VelocityDistributions.o \
       pic2d_Setup.o \
       pic2d_WallPotentials.o \
       pic2d_ParticleExchange.o \
       pic2d_BlockSetProc.o \
       pic2d_IonMoments.o \
       pic2d_Checkpoints.o \
       pic2d_IonWallCollisions.o \
       pic2d_LoadBalancing.o \
       pic2d_Diagnostics.o \
       pic2d_MainProgram.o \
       pic2d_ElectricFieldCalc_FFT_X.o \
       pic2d_Materials.o \
       pic2d_ElectronDynamics.o \
       pic2d_ElectronDynamicsWithVDFMoments.o \
       pic2d_ElectronExchange.o \
       pic2d_ElectronMoments.o \
       pic2d_ElectronWallCollisions.o \
       pic2d_enCollisionsGeneralProc.o \
       pic2d_enCollisionsSpecificProc.o \
       pic2d_ExternalFields.o \
       pic2d_HTSetup.o \
       pic2d_inCollisionsGeneralProc.o \
       pic2d_inCollisionsSpecificProc.o \
       pic2d_IonDynamics.o \
       pic2d_IonDynamicsWithVDFMoments.o \
       pic2d_PETSc_Solver.o \
       pic2d_ElectricFieldCalc_PETSc.o \
       pic2d_CurProblemValues.o \
       pic2d_SEEmission.o \
       pic2d_IonEmissions.o \
       pic2d_ExternalCircuits.o \
       pic2d_ExternalCircuit_PS_Voltages.o \
       pic2d_ElectronCoulombScattering_rescale_2.o

edipic2d : version_info $(OBJ) 
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIBS)
	@sed -e 's|\"$(GIT_BRANCH)\"|\"GIT_BRANCH\"|g' \
	-e 's|\"$(GIT_HASH)\"|\"GIT_HASH\"|g' \
	-e 's|\"$(GIT_DATE)\"|\"GIT_DATE\"|g' \
	pic2d_Modules.f90  > temp_version
	@mv temp_version pic2d_Modules.f90

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 

%.o : %.F90
	$(FC) $(FFLAGS) -c $< 

%.o : %.c
	$(CC) $(CFLAGS) -c $<

clean: clean_git_info
	rm -f *.mod *.o edipic2d

clean_git_info:
	@sed -e 's|\"$(GIT_BRANCH)\"|\"GIT_BRANCH\"|g' \
	-e 's|\"$(GIT_HASH)\"|\"GIT_HASH\"|g' \
	-e 's|\"$(GIT_DATE)\"|\"GIT_DATE\"|g' \
	pic2d_Modules.f90  > temp_version
	@mv temp_version pic2d_Modules.f90

version_info:
	@echo ""   
	@echo "     Branch: $(GIT_BRANCH)"   
	@echo "     Commit ID: $(GIT_HASH)"
	@echo "     Commit date: $(GIT_DATE)"
	@sed -e 's|\"GIT_BRANCH\"|\"$(GIT_BRANCH)\"|g' \
	-e 's|\"GIT_HASH\"|\"$(GIT_HASH)\"|g' \
	-e 's|\"GIT_DATE\"|\"$(GIT_DATE)\"|g' \
	pic2d_Modules.f90 > temp_version
	@mv temp_version pic2d_Modules.f90   
	@echo ""

   # @echo "     Branch: $(GIT_BRANCH)"; #echo "     Commit ID: $(GIT_HASH)": echo "     Commit date: $(GIT_DATE)"; echo "     Compilation time: $(COMPILE_TIME)";