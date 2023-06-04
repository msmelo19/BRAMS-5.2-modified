include ../../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(imogen_anlg_vals.o)  $(ARC)(imogen_constants.o) $(ARC)(imogen_run.o) \
      $(ARC)(imogen_time.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(imogen_anlg_vals.o): $(THISFILE)
$(ARC)(imogen_constants.o): $(THISFILE)
$(ARC)(imogen_run.o): $(THISFILE)
$(ARC)(imogen_time.o): $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

# End of file.
