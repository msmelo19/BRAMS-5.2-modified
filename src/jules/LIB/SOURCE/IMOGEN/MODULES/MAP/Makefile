include ../../../../Makefile.common.mk
THISFILE = ./Makefile

all : $(ARC)(imogen_map.o) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(imogen_map.o): $(ARC)(imogen_constants.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

# End of file.
