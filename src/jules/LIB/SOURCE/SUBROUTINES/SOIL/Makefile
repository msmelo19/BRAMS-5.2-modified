include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(calc_baseflow_jules.o)   $(ARC)(darcy_ch.o)        $(ARC)(darcy_ic.o)      \
      $(ARC)(darcy_vg.o)        $(ARC)(heat_con.o)        $(ARC)(ice_htc.o)      \
      $(ARC)(sieve.o)           $(ARC)(calc_fsat.o)       $(ARC)(frunoff.o)      \
      $(ARC)(hyd_con_ch.o)      $(ARC)(hyd_con_ic.o)      $(ARC)(hyd_con_vg.o)     \
      $(ARC)(pdm.o)             $(ARC)(soil_htc.o)        $(ARC)(soilt.o)         \
      $(ARC)(ch4_wetl.o)        $(ARC)(gauss.o)           $(ARC)(hydrol.o)      \
      $(ARC)(soil_hyd.o)        $(ARC)(soilmc.o)      \
      $(ARC)(surf_hyd.o)        $(ARC)(routing.o)         $(ARC)(calc_zw.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(calc_baseflow_jules.o): $(ARC)(c_topog.o) $(THISFILE)
$(ARC)(calc_fsat.o): $(ARC)(c_topog.o) $(THISFILE)
$(ARC)(calc_zw.o): $(ARC)(c_densty.o) $(ARC)(c_topog.o) $(THISFILE)
$(ARC)(ch4_wetl.o): $(ARC)(c_ch4.o) $(THISFILE)
$(ARC)(darcy_ch.o): $(THISFILE)
$(ARC)(darcy_ic.o): $(ARC)(switches.o) $(THISFILE)
$(ARC)(darcy_vg.o): $(THISFILE)
$(ARC)(frunoff.o): $(THISFILE)
$(ARC)(gauss.o): $(THISFILE)
$(ARC)(heat_con.o): $(ARC)(snow_param.o) $(ARC)(soil_param.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(hyd_con_ch.o): $(THISFILE)
$(ARC)(hyd_con_ic.o): $(ARC)(switches.o) $(THISFILE)
$(ARC)(hyd_con_vg.o): $(THISFILE)
$(ARC)(hydrol.o): $(ARC)(c_topog.o) $(ARC)(soil_param.o) $(THISFILE)
$(ARC)(ice_htc.o): $(ARC)(snow_param.o) $(THISFILE)
$(ARC)(pdm.o): $(ARC)(c_pdm.o) $(ARC)(c_densty.o) $(ARC)(soil_param.o) $(THISFILE)
$(ARC)(routing.o): $(ARC)(ancil_info.o) $(ARC)(drive_io_vars.o) $(ARC)(earth_utils.o) \
                   $(ARC)(fluxes.o) $(ARC)(grid_utils.o)  $(ARC)(offline_diag.o) \
                   $(ARC)(route_mod.o) $(ARC)(switches.o) $(ARC)(top_pdm.o) $(ARC)(timeconst.o) $(THISFILE)
$(ARC)(sieve.o): $(THISFILE)
$(ARC)(soil_htc.o): $(ARC)(c_0_dg_c.o) $(ARC)(c_perma.o) $(ARC)(c_lheat.o) \
                    $(ARC)(c_densty.o) $(ARC)(snow_param.o) $(ARC)(soil_param.o) $(THISFILE)
$(ARC)(soil_hyd.o): $(ARC)(c_topog.o) $(ARC)(c_densty.o) $(ARC)(soil_param.o) $(THISFILE)
$(ARC)(soilmc.o): $(ARC)(soil_param.o) $(ARC)(c_densty.o) $(THISFILE)
$(ARC)(soilt.o): $(ARC)(soil_param.o) $(THISFILE)
$(ARC)(surf_hyd.o): $(ARC)(soil_param.o) $(ARC)(switches.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.

