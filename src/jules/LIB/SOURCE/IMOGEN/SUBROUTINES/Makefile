include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(clim_calc.o)          $(ARC)(day_calc.o)             $(ARC)(delta_temp.o)         \
      $(ARC)(diffcarb_land.o)      $(ARC)(drdat.o)                $(ARC)(gcm_anlg.o)           \
      $(ARC)(imogen_check.o)       $(ARC)(imogen_confirmed_run.o) $(ARC)(imogen_update_carb.o) \
      $(ARC)(imogen_update_clim.o) $(ARC)(invert.o)               $(ARC)(ocean_co2.o)          \
      $(ARC)(radf_co2.o)           $(ARC)(radf_non_co2.o)         $(ARC)(redis.o)              \
      $(ARC)(response.o)           $(ARC)(rndm.o)                 $(ARC)(solang.o)             \
      $(ARC)(solpos.o)             $(ARC)(sunny.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(clim_calc.o): $(THISFILE)
$(ARC)(day_calc.o): $(ARC)(c_pi.o) $(THISFILE)
$(ARC)(delta_temp.o): $(THISFILE)
$(ARC)(diffcarb_land.o): $(ARC)(c_pi.o) $(ARC)(earth_constants_mod.o) $(THISFILE)
$(ARC)(drdat.o): $(ARC)(file_utils.o) $(ARC)(inout.o) $(THISFILE)
$(ARC)(gcm_anlg.o): $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(imogen_map.o) \
                    $(ARC)(imogen_constants.o) $(THISFILE)
$(ARC)(imogen_check.o): $(THISFILE)
$(ARC)(imogen_confirmed_run.o): $(THISFILE)
$(ARC)(imogen_update_carb.o): $(ARC)(ancil_info.o) $(ARC)(trifctl.o) $(ARC)(prognostics.o)      \
                              $(ARC)(imogen_run.o) $(ARC)(imogen_clim.o) $(ARC)(imogen_progs.o) \
                              $(ARC)(imogen_constants.o) $(ARC)(imogen_time.o)                  \
                              $(ARC)(imogen_anlg_vals.o) $(ARC)(imogen_io_vars.o)               \
                              $(THISFILE)
$(ARC)(imogen_update_clim.o): $(ARC)(ancil_info.o) $(ARC)(trifctl.o) $(ARC)(prognostics.o) \
                              $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(aero.o) $(ARC)(timeconst.o) \
                              $(ARC)(imogen_run.o) $(ARC)(imogen_clim.o) $(ARC)(imogen_progs.o) \
                              $(ARC)(imogen_constants.o) $(ARC)(imogen_time.o)  \
                              $(ARC)(imogen_anlg_vals.o) $(ARC)(imogen_drive_vars.o) \
                              $(THISFILE)
$(ARC)(invert.o): $(THISFILE)
$(ARC)(ocean_co2.o): $(THISFILE)
$(ARC)(radf_co2.o): $(THISFILE)
$(ARC)(radf_non_co2.o): $(ARC)(file_utils.o) $(ARC)(inout.o) $(THISFILE)
$(ARC)(redis.o): $(THISFILE)
$(ARC)(response.o): $(THISFILE)
$(ARC)(rndm.o): $(THISFILE)
$(ARC)(solang.o): $(ARC)(c_pi.o) $(THISFILE)
$(ARC)(solpos.o): $(ARC)(c_pi.o) $(THISFILE)
$(ARC)(sunny.o): $(ARC)(c_pi.o) $(ARC)(timeconst.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

# End of file.
