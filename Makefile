.DEFAULT_GOAL := all

DIR = .

TGT = main.exe

FC = mpif90

OBJDIR = obj
SRCDIR = src
INCLUDEDIR = inc

MOD = -J$(DIR)/$(INCLUDEDIR)

SRC =  param.f90 \
       mpi.f90 \
       basis.f90 \
       mesh.f90 \
       graph_partition.f90 \
       write_fields.f90 \
       user_defined.f90 \
       DG_wave.f90 \
       main_loop.f90 

#SOURCE = $(wildcard $(DIR)/$(SRCDIR)/$(SRC))
#SOURCE  = $(DIR)/$(SRCDIR)
SOURCE = $(patsubst %, $(DIR)/$(SRCDIR)/%, $(SRC))

OBJ = $(addprefix $(DIR)/$(OBJDIR)/, $(notdir $(SRC:.f90=.o)))


$(DIR)/$(OBJDIR)/$(TGT) : $(OBJ)
#	$(FC) $(MOD) -o $@ $^
	$(FC) $(MOD) -o $(TGT) $^
 
$(DIR)/$(OBJDIR)/%.o : $(DIR)/$(SRCDIR)/%.f90
	$(FC) $(MOD) -c $< -o $@



.PHONY : help run clean all 

all : $(DIR)/$(OBJDIR)/$(TGT)
#	$(DIR)/$(OBJDIR)/$(TGT)
	@echo "------------------------------"
	@echo "Makefile succeed"
	@echo "------------------------------"

run : $(TGT)
	mpirun -np 1 $(TGT)

help : 
	@echo "source : $(SOURCE)"
	@echo "src : $(SRC)"
	@echo "obj : $(OBJ)"


clean :
	rm -rf $(OBJ) 
#	rm -rf $(DIR)/$(OBJDIR)/$(TGT) 
	rm -rf $(DIR)/$(OBJDIR)/*.mod
	rm -rf *.dat *.txt $(TGT)
