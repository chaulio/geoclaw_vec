# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Version control:
VERSION?=1
# 0 = original code (without statistics)
# 1 = original code (with statistics)
# 2 = SOA 


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= 

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw
NEWCODE = ../src

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

ifeq ($(VERSION),1)
EXCLUDE_SOURCES = \
  $(GEOLIB)/amr2.f90 \

endif

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------


MODULES = \
  $(NEWCODE)/statistics.f90 \

# original code without statistics
ifeq ($(VERSION),0)
SOURCES = \
   $(CLAW)/riemann/src/rpn2_geoclaw.f \
   $(CLAW)/riemann/src/rpt2_geoclaw.f \
   $(CLAW)/riemann/src/geoclaw_riemann_utils.f \

endif

# original code WITH statistics
ifeq ($(VERSION),1)
SOURCES = \
  $(NEWCODE)/timer/rpn2_geoclaw.f \
  $(NEWCODE)/timer/amr2.f90 \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \

endif
  
  
  

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
.PHONY: topo all
topo:
	python maketopo.py

all: 
	$(MAKE) topo
	$(MAKE) .plots
	$(MAKE) .htmls

