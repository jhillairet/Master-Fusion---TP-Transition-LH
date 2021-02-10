.SUFFIXES: .o .f90 .F90 .f .swp
#
#  Makefile for gyro2D
#
SRCS_f90 = prec_const.f90 globals.f90 read_data.f90 initialisation.f90 \
           tools.f90 diffusionNL.f90 evolution.f90 \
           sauvegarde.f90  transitionLH.f90

SRCS_F90 = 

SRCS =	$(SRCS_f90) $(SRCS_F90)

OBJS =  prec_const.o globals.o read_data.o initialisation.o \
        tools.o diffusionNL.o evolution.o \
        sauvegarde.o  transitionLH.o

F77 = gfortran -c -fdefault-real-8
F95 = gfortran -c -fdefault-real-8
LD  = gfortran -fdefault-real-8
FFLAGS   = -O2  
#option a utiliser si pas de totalview 
F95FLAGS = -O2 
  
#LIB = -ldxml

.f90.o:	
	$(F95) $(F95FLAGS) $<
.f.o:	
	$(F77) $(FFLAGS) $<
vlasm:	$(OBJS)
	$(LD) $(F95FLAGS) -o transitionLH $(OBJS) $(LIB)
clean:	
	rm -f $(OBJS) *~ *.mod *.kmo *.l core
debug:	
	make F95FLAGS=-g FFLAGS=-g
