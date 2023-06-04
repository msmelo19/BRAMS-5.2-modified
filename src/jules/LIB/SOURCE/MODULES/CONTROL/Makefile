include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(aero.o)          $(ARC)(ancil_info.o)   $(ARC)(blopt8a.o)      \
      $(ARC)(fluxes.o)        $(ARC)(grid_utils.o)   $(ARC)(c_elevate.o)    \
      $(ARC)(orog.o)          $(ARC)(sea_ice.o)      $(ARC)(time_mod.o)     \
      $(ARC)(jules_netcdf.o)  $(ARC)(ancil_info.o)   $(ARC)(forcing.o)      \
      $(ARC)(output_mod.o)    $(ARC)(prognostics.o)                         \
      $(ARC)(top_pdm.o)       $(ARC)(coastal.o)      $(ARC)(p_s_parms.o)    \
      $(ARC)(trifctl.o)       $(ARC)(diag_swchs.o)   $(ARC)(initial_mod.o)  \
      $(ARC)(readwrite_mod.o) $(ARC)(spin_mod.o)     $(ARC)(update_mod.o)   \
      $(ARC)(drive_io_vars.o) $(ARC)(inout.o)        $(ARC)(route_mod.o)    \
      $(ARC)(switches.o)      $(ARC)(u_v_grid.o)     $(ARC)(file_utils.o)   \
      $(ARC)(misc_utils.o)    $(ARC)(screen.o)       $(ARC)(time_loc.o)     \
      $(ARC)(veg_io_vars.o)   $(ARC)(earth_utils.o)  $(ARC)(offline_diag.o) \
      $(ARC)(rwerr_mod.o)     $(ARC)(switches_urban.o) $(ARC)(ozone_vars.o) \
      $(ARC)(solinc_data.o)   $(ARC)(um_types.o)     $(ARC)(jules_mod.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(aero.o): $(THISFILE)
$(ARC)(ancil_info.o): $(THISFILE)
$(ARC)(blopt8a.o): $(THISFILE)
$(ARC)(coastal.o): $(THISFILE)
$(ARC)(c_elevate.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(diag_swchs.o):  $(THISFILE)
$(ARC)(drive_io_vars.o): $(ARC)(inout.o) $(THISFILE)
$(ARC)(earth_utils.o):  $(THISFILE)
$(ARC)(file_utils.o): $(ARC)(inout.o) $(ARC)(jules_netcdf.o) $(THISFILE)
$(ARC)(fluxes.o):  $(THISFILE)
$(ARC)(forcing.o): $(THISFILE)
$(ARC)(grid_utils.o):  $(ARC)(ancil_info.o) $(THISFILE)
$(ARC)(initial_mod.o): $(ARC)(ancil_info.o) $(ARC)(c_densty.o)   \
                   $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(soil_param.o) \
                   $(ARC)(misc_utils.o) $(ARC)(nvegparm.o) $(ARC)(p_s_parms.o)  \
                   $(ARC)(pftparm.o) $(ARC)(prognostics.o) $(ARC)(nstypes.o)   \
                   $(ARC)(route_mod.o) $(ARC)(switches.o) $(ARC)(readwrite_mod.o)  \
                   $(ARC)(spin_mod.o) $(ARC)(time_loc.o) $(ARC)(snow_param.o)   \
                   $(ARC)(top_pdm.o) $(ARC)(switches_urban.o) $(ARC)(urban_param.o) $(THISFILE)
$(ARC)(inout.o): $(THISFILE)
$(ARC)(jules_netcdf.o): $(ARC)(inout.o) $(ARC)(rwerr_mod.o) $(ARC)(switches.o) $(ARC)(time_loc.o) $(THISFILE)
$(ARC)(misc_utils.o): $(ARC)(inout.o) $(ARC)(file_utils.o) $(ARC)(switches.o) $(ARC)(time_loc.o)  \
                      $(ARC)(time_mod.o) $(ARC)(timeconst.o) $(THISFILE)
$(ARC)(offline_diag.o): $(THISFILE)
$(ARC)(orog.o):  $(THISFILE)
$(ARC)(output_mod.o): $(ARC)(ancil_info.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_densty.o)  \
                  $(ARC)(csigma.o) $(ARC)(file_utils.o)   \
                  $(ARC)(fluxes.o) $(ARC)(forcing.o) $(ARC)(grid_utils.o) $(ARC)(inout.o)  \
                  $(ARC)(misc_utils.o) $(ARC)(nvegparm.o) $(ARC)(nstypes.o)  \
                  $(ARC)(offline_diag.o)  \
                  $(ARC)(p_s_parms.o) $(ARC)(pftparm.o) $(ARC)(prognostics.o)   \
                  $(ARC)(readwrite_mod.o) $(ARC)(route_mod.o) $(ARC)(screen.o)  \
                  $(ARC)(snow_param.o) $(ARC)(spin_mod.o) $(ARC)(switches.o)  \
                  $(ARC)(timeconst.o) $(ARC)(surf_param.o) $(ARC)(soil_param.o)  \
                  $(ARC)(time_loc.o) $(ARC)(time_mod.o) $(ARC)(trifctl.o)     \
                  $(ARC)(top_pdm.o) $(ARC)(ozone_vars.o) $(THISFILE)
$(ARC)(ozone_vars.o): $(THISFILE)
$(ARC)(prognostics.o): $(THISFILE)
$(ARC)(p_s_parms.o):  $(THISFILE)
$(ARC)(readwrite_mod.o): $(ARC)(inout.o) $(ARC)(grid_utils.o) $(ARC)(jules_netcdf.o) $(THISFILE)
$(ARC)(rwerr_mod.o): $(ARC)(inout.o) $(THISFILE)
$(ARC)(route_mod.o): $(ARC)(inout.o) $(THISFILE)
$(ARC)(screen.o): $(THISFILE)
$(ARC)(sea_ice.o): $(THISFILE)
$(ARC)(spin_mod.o): $(THISFILE)
$(ARC)(switches.o):  $(THISFILE)
$(ARC)(time_loc.o): $(THISFILE)
$(ARC)(time_mod.o): $(ARC)(inout.o) $(ARC)(timeconst.o) $(THISFILE)
$(ARC)(top_pdm.o):  $(THISFILE)
$(ARC)(trifctl.o):  $(THISFILE)
$(ARC)(u_v_grid.o):  $(THISFILE)
$(ARC)(update_mod.o): $(ARC)(ancil_info.o) $(ARC)(drive_io_vars.o) $(ARC)(file_utils.o)  \
                      $(ARC)(fluxes.o)     $(ARC)(forcing.o)   $(ARC)(grid_utils.o)      \
                      $(ARC)(initial_mod.o) $(ARC)(inout.o) $(ARC)(misc_utils.o)         \
                      $(ARC)(pftparm.o)    $(ARC)(timeconst.o)  \
                      $(ARC)(p_s_parms.o)  $(ARC)(prognostics.o) $(ARC)(readwrite_mod.o)  \
                      $(ARC)(spin_mod.o) $(ARC)(switches.o) $(ARC)(time_loc.o)  \
                      $(ARC)(time_mod.o) $(ARC)(u_v_grid.o) $(ARC)(veg_io_vars.o)  \
                      $(ARC)(ozone_vars.o) $(ARC)(imogen_drive_vars.o) $(THISFILE)
$(ARC)(veg_io_vars.o): $(ARC)(inout.o) $(THISFILE)
$(ARC)(solinc_data.o): $(THISFILE)
$(ARC)(um_types.o): $(THISFILE)
$(ARC)(jules_mod.o): $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

# End of file.
