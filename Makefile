# Makefile for BPLBM3D_NEW
# On the HPC, first run:
#   module load gcc/15.2.0
# Then compile with:
#   make

FC = gfortran
EXE ?= BPLBM3D

OBJDIR := build
MODDIR := $(OBJDIR)/mod

# Source files
FREE_SRCS  := MRT.f90 BPLBM3D.f90
FIXED_SRCS := nbsw3d.f

# Object files in correct dependency/link order
OBJS := $(OBJDIR)/MRT.o \
        $(OBJDIR)/nbsw3d.o \
        $(OBJDIR)/BPLBM3D.o

# Optimization level can be overridden, e.g. make OPT="-O3 -march=native"
OPT ?= -O2

# Common Fortran flags
COMMON_FFLAGS := $(OPT) -I$(MODDIR) -J$(MODDIR) -Wno-tabs

# Free-form and fixed-form source flags
FREE_FFLAGS  := -ffree-line-length-none
FIXED_FFLAGS := -ffixed-line-length-none

.PHONY: all clean veryclean debug run

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(COMMON_FFLAGS) -o $@ $^

# BPLBM3D uses module MRT, so it must compile after MRT.o/mrt.mod exists.
$(OBJDIR)/BPLBM3D.o: BPLBM3D.f90 $(OBJDIR)/MRT.o | $(OBJDIR) $(MODDIR)
	$(FC) $(COMMON_FFLAGS) $(FREE_FFLAGS) -c $< -o $@

$(OBJDIR)/MRT.o: MRT.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(COMMON_FFLAGS) $(FREE_FFLAGS) -c $< -o $@

$(OBJDIR)/nbsw3d.o: nbsw3d.f | $(OBJDIR) $(MODDIR)
	$(FC) $(COMMON_FFLAGS) $(FIXED_FFLAGS) -c $< -o $@

$(OBJDIR) $(MODDIR):
	mkdir -p $@

debug:
	$(MAKE) clean
	$(MAKE) OPT="-O0 -g -fcheck=all -fbacktrace -Wall -Wextra"

run: $(EXE)
	./$(EXE)

clean:
	rm -rf $(OBJDIR)

veryclean: clean
	rm -f $(EXE)
