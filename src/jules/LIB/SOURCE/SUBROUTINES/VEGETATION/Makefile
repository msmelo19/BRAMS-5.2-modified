include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(compete.o)    $(ARC)(growth.o)   $(ARC)(phenol.o) \
      $(ARC)(soilcarb.o)   $(ARC)(tilepts.o)  $(ARC)(triffid.o)  \
      $(ARC)(veg-veg2a.o)  $(ARC)(decay.o)    $(ARC)(lotka.o)  \
      $(ARC)(pft_sparm.o)  $(ARC)(sparm.o)    $(ARC)(veg-veg1a.o) \
      $(ARC)(vegcarb.o)    $(ARC)(dpm_rpm.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(compete.o): $(ARC)(descent.o) $(ARC)(nstypes.o) $(ARC)(seed.o) $(THISFILE)
$(ARC)(decay.o): $(ARC)(descent.o) $(ARC)(csmin.o) $(THISFILE)
$(ARC)(dpm_rpm.o): $(ARC)(nstypes.o) $(THISFILE)
$(ARC)(growth.o): $(ARC)(descent.o) $(ARC)(pftparm.o) $(ARC)(trif.o) $(THISFILE)
$(ARC)(lotka.o): $(ARC)(nstypes.o) $(ARC)(pftparm.o) $(ARC)(seed.o) $(ARC)(sigm.o) \
                 $(ARC)(trif.o) $(THISFILE)
$(ARC)(pft_sparm.o): $(ARC)(pftparm.o) $(THISFILE)
$(ARC)(phenol.o): $(ARC)(pftparm.o) $(ARC)(trif.o) $(THISFILE)
$(ARC)(soilcarb.o): $(ARC)(nstypes.o) $(THISFILE)
$(ARC)(sparm.o): $(ARC)(nstypes.o) $(ARC)(blend_h.o) $(ARC)(nvegparm.o) \
                 $(ARC)(pftparm.o) $(ARC)(snow_param.o) $(ARC)(urban_param.o) \
                 $(THISFILE)
$(ARC)(tilepts.o): $(ARC)(nstypes.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(triffid.o): $(ARC)(nstypes.o)  $(ARC)(pftparm.o) $(ARC)(trif.o) $(THISFILE)
$(ARC)(veg-veg1a.o): $(ARC)(nstypes.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(veg-veg2a.o): $(ARC)(nstypes.o) $(ARC)(descent.o) $(ARC)(seed.o) \
                     $(ARC)(veg_param.o) $(ARC)(c_mdi.o) $(ARC)(switches.o) \
                     $(THISFILE)
$(ARC)(vegcarb.o): $(ARC)(pftparm.o) $(ARC)(trif.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.
