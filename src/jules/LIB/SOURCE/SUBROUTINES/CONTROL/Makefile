include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(allocate_arrays.o) $(ARC)(control.o)    $(ARC)(deallocate_arrays.o)  \
      $(ARC)(do_areaver.o)      $(ARC)(jules.o)      $(ARC)(jules_final.o)        \
      $(ARC)(newtime.o)         $(ARC)(pre_areaver.o) $(ARC)(spin_check.o)         \
      $(ARC)(timer.o)           $(ARC)(zenith_mod.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(allocate_arrays.o): $(ARC)(aero.o) $(ARC)(ancil_info.o) $(ARC)(dust_param.o)  \
                           $(ARC)(c_z0h_z0m.o) $(ARC)(coastal.o) $(ARC)(drive_io_vars.o)   \
                           $(ARC)(fluxes.o) $(ARC)(forcing.o) $(ARC)(inout.o)          \
                           $(ARC)(misc_utils.o) $(ARC)(nstypes.o) $(ARC)(nvegparm.o)   \
                           $(ARC)(offline_diag.o) $(ARC)(orog.o) $(ARC)(p_s_parms.o)      \
                           $(ARC)(pftparm.o) $(ARC)(prognostics.o) $(ARC)(route_mod.o) $(ARC)(screen.o)   \
                           $(ARC)(spin_mod.o) $(ARC)(switches.o) $(ARC)(c_elevate.o)  \
                           $(ARC)(time_loc.o)  $(ARC)(trif.o) $(ARC)(trifctl.o)  \
                           $(ARC)(top_pdm.o) $(ARC)(u_v_grid.o) $(ARC)(veg_io_vars.o)  \
                           $(ARC)(soil_param.o) $(ARC)(surf_param.o) $(ARC)(snow_param.o)  \
                           $(ARC)(switches_urban.o) $(ARC)(urban_param.o) $(ARC)(ozone_vars.o) \
                           $(ARC)(imogen_time.o) $(ARC)(imogen_clim.o) $(ARC)(imogen_drive_vars.o) \
                           $(ARC)(jules_mod.o) $(THISFILE)
$(ARC)(control.o): $(ARC)(c_epslon.o) $(ARC)(csigma.o) $(ARC)(ancil_info.o)  \
                   $(ARC)(drive_io_vars.o) $(ARC)(inout.o)   \
                   $(ARC)(switches.o) $(ARC)(diag_swchs.o)  $(ARC)(pftparm.o) \
                   $(ARC)(screen.o) $(ARC)(p_s_parms.o) $(ARC)(forcing.o)  \
                   $(ARC)(prognostics.o) $(ARC)(fluxes.o) $(ARC)(aero.o) $(ARC)(nstypes.o)  \
                   $(ARC)(orog.o) $(ARC)(u_v_grid.o) $(ARC)(trifctl.o) $(ARC)(timeconst.o)  \
                   $(ARC)(top_pdm.o) $(ARC)(coastal.o) $(ARC)(sea_ice.o)  \
                   $(ARC)(snow_param.o) $(ARC)(surf_param.o) $(ARC)(zenith_mod.o) \
                   $(ARC)(urban_param.o) $(ARC)(ozone_vars.o) $(THISFILE)
$(ARC)(deallocate_arrays.o): $(ARC)(aero.o) $(ARC)(ancil_info.o) $(ARC)(coastal.o)  \
                             $(ARC)(drive_io_vars.o) $(ARC)(fluxes.o) $(ARC)(forcing.o)     \
                             $(ARC)(misc_utils.o) $(ARC)(offline_diag.o) $(ARC)(orog.o) $(ARC)(p_s_parms.o)  \
                             $(ARC)(prognostics.o) $(ARC)(route_mod.o) $(ARC)(screen.o) $(ARC)(spin_mod.o) \
                             $(ARC)(time_loc.o) $(ARC)(top_pdm.o) \
                             $(ARC)(trifctl.o) $(ARC)(u_v_grid.o) $(ARC)(veg_io_vars.o)  \
                             $(ARC)(soil_param.o) $(ARC)(surf_param.o) $(ARC)(snow_param.o)  \
                             $(ARC)(urban_param.o) $(ARC)(ozone_vars.o) \
                             $(ARC)(imogen_clim.o) $(ARC)(imogen_drive_vars.o) $(ARC)(jules_mod.o) \
                             $(THISFILE)
$(ARC)(do_areaver.o): $(ARC)(c_pi.o) $(THISFILE)
$(ARC)(jules.o): $(ARC)(inout.o) $(ARC)(output_mod.o) $(ARC)(spin_mod.o) $(ARC)(time_loc.o) \
                 $(ARC)(time_mod.o) $(ARC)(trifctl.o)  \
                 $(ARC)(update_mod.o) $(ARC)(veg_io_vars.o) $(THISFILE)
$(ARC)(jules_FINAL.o): $(ARC)(initial_mod.o) $(ARC)(inout.o) $(ARC)(output_mod.o) $(THISFILE)
$(ARC)(newtime.o): $(ARC)(drive_io_vars.o) $(ARC)(initial_mod.o) $(ARC)(inout.o) $(ARC)(misc_utils.o)  \
                 $(ARC)(route_mod.o) $(ARC)(spin_mod.o) $(ARC)(switches.o)  \
                 $(ARC)(time_loc.o) $(ARC)(time_mod.o)  \
                 $(ARC)(update_mod.o) $(ARC)(veg_io_vars.o) \
                 $(ARC)(imogen_time.o) $(THISFILE)
$(ARC)(pre_areaver.o): $(ARC)(c_pi.o) $(THISFILE)
$(ARC)(spin_check.o): $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(prognostics.o)  \
                      $(ARC)(route_mod.o)  \
                      $(ARC)(spin_mod.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(timer.o): $(THISFILE)
$(ARC)(zenith_mod.o): $(ARC)(switches.o) $(ARC)(c_pi.o) $(ARC)(ancil_info.o)  \
                  $(ARC)(time_loc.o) $(ARC)(timeconst.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.
