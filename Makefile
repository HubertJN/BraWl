#####################################
# Makefile for Monte Carlo B-W code #
#####################################

# Specify a particular compiler with "make compiler=pgi", etc.
# Specify debugging flags with "make MODE=debug"
# Alternatively, these options can be specified as environment variables
# eg. "export compiler=gfortran" can be added to $HOME/.bashrc

SYSTEM=$(shell uname -s)

# Compiler specific flags

# Note: you MUST specify a compiler option. None are specified by default.

ifeq ($(strip $(compiler)),)
  MAKECMDGOALS = error
error:
	@echo '*** ERROR ***'
	@echo 'You MUST set a value for the compiler variable'
	@echo ' eg. "make compiler=intel"'
	@echo 'Alternatively, you can add "export compiler=intel" to $$HOME/.bashrc'
endif

# Intel
# =====
ifeq ($(strip $(compiler)),intel)
  FC = mpiifort
#  FFLAGS = -O3 
  FFLAGS = -O0 -warn all -check all
  FFLAGS += -module $(OBJDIR)
  CC = icc
  CFLAGS = -O3
endif

# gfortran
# ========
ifeq ($(strip $(compiler)),gfortran)
  FC = mpif90
  FFLAGS = -O0 -Wall -Wextra -fcheck=bounds
  FFLAGS += -I/usr/local/include -I$(OBJDIR) -J$(OBJDIR)
  LDFLAGS=-lgcc
  CC=gcc -I$(INCDIR)
  CFLAGS=-O3
endif

SRCDIR=src
BINDIR=bin
DATDIR=data
OBJDIR=obj
INCDIR=include

ifeq ($(SYSTEM),Darwin)
         FFLAGS +=-I$(shell nf-config --fflags)
         LDFLAGS += $(shell nf-config --flibs) \
                    -lnetcdf -lnetcdff
else
         FFLAGS +=$(shell nf-config --fflags)
         LDFLAGS += $(shell nf-config --flibs) 
endif

# Command to use for linking and executable
LD=$(FC)
EXE=alloy_mc.exe

SRCFILES=mt19937ar.c kinds.f90 mpi_shared_data.f90 io.f90 comms.f90 write_netcdf.f90 write_xyz.f90\
	write_diagnostics.f90 command_line.f90 c_functions.f90 display.f90\
        energetics.f90 analytics.f90 random_site.f90 metropolis.f90\
	initialise.f90 model_run.f90 mpi_main.f90

OBJFILES:=$(SRCFILES:.f90=.o)
OBJFILES:=$(OBJFILES:.c=.o)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(OBJDIR):$(INCDIR)

alloy: $(OBJFILES)
	$(FC) $(FFLAGS) -o $(EXE) $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)


# Purge build files and executable
clean :
	@rm -rf $(OBJDIR) $(BINDIR) $(EXE)

# Purge build files and executable
cleardata :
	@rm -r data/grids/*
	@rm -r data/correlations/*
	@rm -r data/order/*

# Rules for building object files
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $(OBJDIR)/$@ $<

# Rules for building object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $(OBJDIR)/$@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(OBJFILES): | $(OBJDIR)

################
# Dependencies #
################
mpi_shared_data.o: kinds.o
io.o: kinds.o mpi_shared_data.o
write_netcdf.o: kinds.o mpi_shared_data.o
write_xyz.o: mpi_shared_data.o kinds.o analytics.o
write_diagnostics.o: mpi_shared_data.o kinds.o analytics.o
command_line.o: kinds.o
c_functions.o: mt19937ar.o
display.o: kinds.o mpi_shared_data.o
energetics.o: kinds.o mpi_shared_data.o c_functions.o
analytics.o: mpi_shared_data.o kinds.o display.o
random_site.o: mpi_shared_data.o kinds.o c_functions.o analytics.o
metropolis.o: kinds.o mpi_shared_data.o c_functions.o energetics.o random_site.o analytics.o
initialise.o: kinds.o mpi_shared_data.o c_functions.o energetics.o random_site.o metropolis.o
mpi_main.o: initialise.o mpi_shared_data.o kinds.o c_functions.o write_netcdf.o\
	write_xyz.o write_diagnostics.o command_line.o display.o model_run.o
