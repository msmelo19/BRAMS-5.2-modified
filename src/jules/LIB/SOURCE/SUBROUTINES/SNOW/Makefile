include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(canopysnow.o) $(ARC)(compactsnow.o) $(ARC)(layersnow.o)   \
      $(ARC)(relayersnow.o) $(ARC)(snow.o) $(ARC)(snowgrain.o)    \
      $(ARC)(snowpack.o) $(ARC)(snowtherm.o) $(ARC)(tridag.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(canopysnow.o):  $(THISFILE)
$(ARC)(compactsnow.o):  $(ARC)(ancil_info.o) $(ARC)(c_g.o) $(THISFILE)
$(ARC)(layersnow.o): $(ARC)(ancil_info.o) $(ARC)(snow_param.o) $(THISFILE)
$(ARC)(relayersnow.o): $(ARC)(ancil_info.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_perma.o) $(THISFILE)
$(ARC)(snow.o): $(ARC)(ancil_info.o) $(ARC)(c_lheat.o) $(ARC)(snow_param.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(snowgrain.o): $(ARC)(ancil_info.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_pi.o)   \
                     $(ARC)(rad_param.o)  $(THISFILE)
$(ARC)(snowpack.o): $(ARC)(ancil_info.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_densty.o) $(ARC)(c_perma.o) \
                    $(ARC)(c_lheat.o) $(ARC)(snow_param.o) $(ARC)(soil_param.o) $(THISFILE)
$(ARC)(snowtherm.o): $(ARC)(ancil_info.o) $(ARC)(c_perma.o) $(ARC)(c_rhowat.o) $(THISFILE)
$(ARC)(tridag.o): $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.
