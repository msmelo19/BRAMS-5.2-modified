include ../../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(imogen_clim.o) $(ARC)(imogen_drive_vars.o) $(ARC)(imogen_io_vars.o) \
      $(ARC)(imogen_progs.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(imogen_clim.o): $(THISFILE)
$(ARC)(imogen_drive_vars.o): $(THISFILE)
$(ARC)(imogen_io_vars.o): $(THISFILE)
$(ARC)(imogen_progs.o): $(ARC)(imogen_constants.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

# End of file.
