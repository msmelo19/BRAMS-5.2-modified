.SUFFIXES: .a .o .f90 .F90

#.f.a:
#	$(FC) -c $(FFLAGS) $< $(MOD_PATH)
#	@$(AR) $@ $(%F:.f=.o)
#	@$(RM) $(%F:.f=.o)

.f90.a:
	$(FC) -c $(FFLAGS) $< $(MOD_PATH)
	@$(AR) $@ $(%F:.f90=.o)
	@$(RM) $(%F:.f90=.o)

.F90.a:
	$(FC) -c $(FFLAGS) $(FPP) $(FPP_DEFS) $< $(MOD_PATH) $(FPP_INC_PATH)
	@$(AR) $@ $(%F:.F90=.o)
	@$(RM) $(%F:.F90=.o)

######################################
## Architecture specific variables. ##
######################################
AR=ar r
EXTRACT=ar x
RM=rm -f
MOD_FSUF=.mod
