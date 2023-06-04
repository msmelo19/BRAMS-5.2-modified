include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(calc_fit_fsat.o)    $(ARC)(freeze_soil.o)       $(ARC)(init.o)               \
      $(ARC)(init_drive.o)       $(ARC)(init_hgt.o)          $(ARC)(init_misc.o)          \
      $(ARC)(init_out.o)         $(ARC)(init_out_check.o)    $(ARC)(init_out_map_mod.o)   \
      $(ARC)(init_out_time.o)    $(ARC)(init_out_var.o)      $(ARC)(init_out_varlist.o)   \
      $(ARC)(init_snow.o)        $(ARC)(init_trif.o)         $(ARC)(init_frac.o)          \
      $(ARC)(init_nonveg.o)      $(ARC)(init_parms.o)        $(ARC)(init_pdm.o)           \
      $(ARC)(init_soil.o)        $(ARC)(init_soil_lut.o)     $(ARC)(init_top.o)           \
      $(ARC)(init_vars_tmp.o)    $(ARC)(init_agric.o)        $(ARC)(init_grid_mod.o)      \
      $(ARC)(init_opts.o)        $(ARC)(init_route.o)        $(ARC)(init_time.o)          \
      $(ARC)(init_veg.o)         $(ARC)(init_urban.o)        $(ARC)(init_imogen.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(calc_fit_fsat.o): $(THISFILE)
$(ARC)(freeze_soil.o): $(ARC)(c_0_dg_c.o) $(ARC)(c_perma.o) $(ARC)(c_densty.o) $(THISFILE)
$(ARC)(init.o): $(ARC)(init_grid_mod.o) $(ARC)(initial_mod.o) $(ARC)(inout.o)   \
                $(ARC)(spin_mod.o) $(ARC)(veg_io_vars.o) $(THISFILE)
$(ARC)(init_agric.o): $(ARC)(ancil_info.o) $(ARC)(file_utils.o) \
                     $(ARC)(inout.o) $(ARC)(misc_utils.o) $(ARC)(readwrite_mod.o) \
                     $(ARC)(switches.o) $(ARC)(trifctl.o) $(THISFILE)
$(ARC)(init_drive.o): $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(update_mod.o)  \
                    $(ARC)(file_utils.o) $(ARC)(timeconst.o) \
                    $(ARC)(inout.o) $(ARC)(misc_utils.o)  \
                    $(ARC)(spin_mod.o) $(ARC)(switches.o) $(ARC)(time_loc.o) $(THISFILE)
$(ARC)(init_frac.o): $(ARC)(ancil_info.o) $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(misc_utils.o)  \
                     $(ARC)(nstypes.o) $(ARC)(readwrite_mod.o) $(ARC)(switches.o) \
                     $(ARC)(switches_urban.o)  $(THISFILE)
$(ARC)(init_grid_mod.o): $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(coastal.o)  \
                      $(ARC)(file_utils.o)  \
                     $(ARC)(grid_utils.o) $(ARC)(inout.o) $(ARC)(misc_utils.o) \
                     $(ARC)(nstypes.o) $(ARC)(readwrite_mod.o) $(ARC)(switches.o)   \
                     $(ARC)(time_loc.o) $(THISFILE)
$(ARC)(init_hgt.o):  $(THISFILE)
$(ARC)(init_misc.o): $(ARC)(aero.o) $(ARC)(csmin.o)   \
                     $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(seed.o)   \
                     $(ARC)(sigm.o) $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(init_nonveg.o): $(ARC)(c_z0h_z0m.o) $(ARC)(file_utils.o) $(ARC)(inout.o)  \
                       $(ARC)(misc_utils.o) $(ARC)(nstypes.o) $(ARC)(nvegparm.o)  \
                       $(ARC)(readwrite_mod.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(init_opts.o): $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(file_utils.o)  \
                     $(ARC)(inout.o) $(ARC)(misc_utils.o)  $(ARC)(nstypes.o )    \
                     $(ARC)(nvegparm.o) $(ARC)(pftparm.o)  \
                     $(ARC)(switches.o) $(THISFILE)
$(ARC)(init_out.o) : $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(file_utils.o)   \
                     $(ARC)(init_out_map_mod.o) $(ARC)(inout.o)     $(ARC)(misc_utils.o)   \
                     $(ARC)(offline_diag.o) \
                     $(ARC)(p_s_parms.o)       $(ARC)(spin_mod.o)   $(ARC)(switches.o)     \
                     $(ARC)(time_loc.o)        $(ARC)(time_mod.o)   $(ARC)(timeconst.o) $(THISFILE)
$(ARC)(init_parms.o): $(ARC)(ancil_info.o) $(ARC)(switches.o) $(ARC)(prognostics.o) $(ARC)(p_s_parms.o) \
                      $(ARC)(fluxes.o) $(ARC)(u_v_grid.o) $(ARC)(top_pdm.o) $(ARC)(coastal.o) $(THISFILE)
$(ARC)(init_pdm.o): $(THISFILE)
$(ARC)(init_route.o): $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(file_utils.o)  \
                      $(ARC)(grid_utils.o) $(ARC)(inout.o) $(ARC)(timeconst.o)  \
                      $(ARC)(misc_utils.o) $(ARC)(readwrite_mod.o) $(ARC)(route_mod.o)  \
                      $(ARC)(switches.o) $(ARC)(time_mod.o) $(THISFILE)
$(ARC)(init_snow.o): $(ARC)(c_0_dg_c.o) $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(p_s_parms.o)  \
                     $(ARC)(rad_param.o) $(ARC)(snow_param.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(init_soil.o): $(ARC)(ancil_info.o) $(ARC)(file_utils.o)  \
                     $(ARC)(inout.o) $(ARC)(misc_utils.o)    \
                     $(ARC)(p_s_parms.o) $(ARC)(readwrite_mod.o) $(ARC)(switches.o) $(THISFILE)
$(ARC)(init_soil_lut.o): $(THISFILE)
$(ARC)(init_time.o): $(ARC)(file_utils.o) $(ARC)(inout.o)   \
                     $(ARC)(prognostics.o) $(ARC)(switches.o) $(ARC)(time_loc.o)    \
                     $(ARC)(time_mod.o) $(ARC)(timeconst.o) $(ARC)(trifctl.o) $(THISFILE)
$(ARC)(init_top.o): $(THISFILE)
$(ARC)(init_trif.o): $(ARC)(file_utils.o) $(ARC)(inout.o)  \
                      $(ARC)(misc_utils.o) $(ARC)(nstypes.o) $(ARC)(readwrite_mod.o) $(ARC)(switches.o)  \
                     $(ARC)(trif.o) $(THISFILE)
$(ARC)(init_vars_tmp.o): $(ARC)(ancil_info.o) $(ARC)(top_pdm.o) \
                         $(ARC)(p_s_parms.o) $(ARC)(forcing.o) $(ARC)(prognostics.o)  \
                         $(ARC)(aero.o) $(ARC)(orog.o) $(ARC)(trifctl.o)              \
                         $(ARC)(coastal.o) $(ARC)(c_densty.o)          \
                         $(ARC)(soil_param.o) $(ARC)(c_rough.o) $(ARC)(switches.o)    \
                         $(ARC)(sea_ice.o) $(THISFILE)
$(ARC)(init_veg.o): $(ARC)(c_z0h_z0m.o) $(ARC)(update_mod.o)  \
                    $(ARC)(file_utils.o) $(ARC)(inout.o) $(ARC)(misc_utils.o)        \
                    $(ARC)(nstypes.o) $(ARC)(pftparm.o) $(ARC)(prognostics.o) $(ARC)(readwrite_mod.o)  \
                    $(ARC)(switches.o) $(ARC)(time_loc.o) $(ARC)(time_mod.o) $(ARC)(timeconst.o) \
                    $(ARC)(veg_io_vars.o) $(THISFILE)
$(ARC)(init_out_varlist.o) : $(ARC)(allocate_arrays.o) $(ARC)(inout.o) $(THISFILE)
$(ARC)(init_out_time.o) : $(ARC)(inout.o) $(ARC)(spin_mod.o) $(ARC)(time_loc.o)   \
                          $(ARC)(time_mod.o) $(ARC)(timeconst.o) $(THISFILE)
$(ARC)(init_out_map_mod.o): $(ARC)(allocate_arrays.o) $(ARC)(ancil_info.o) $(ARC)(file_utils.o) \
                            $(ARC)(inout.o)           $(ARC)(grid_utils.o) $(ARC)(misc_utils.o) \
                            $(ARC)(output_mod.o)      $(ARC)(readwrite_mod.o) $(ARC)(route_mod.o)  \
                            $(ARC)(time_loc.o) $(THISFILE)
$(ARC)(init_out_var.o) : $(ARC)(ancil_info.o) $(ARC)(inout.o) $(ARC)(nstypes.o)  \
                         $(ARC)(route_mod.o) $(THISFILE)
$(ARC)(init_out_check.o) : $(ARC)(inout.o) $(ARC)(route_mod.o) $(ARC)(switches.o) $(THISFILE)

$(ARC)(init_urban.o) : $(ARC)(allocate_arrays.o) $(ARC)(inout.o)          $(ARC)(file_utils.o) \
                       $(ARC)(urban_param.o)     $(ARC)(switches_urban.o) $(ARC)(nstypes.o)    \
                       $(ARC)(ancil_info.o)      $(ARC)(switches.o)       $(ARC)(misc_utils.o) \
                       $(ARC)(readwrite_mod.o)   $(THISFILE)
$(ARC)(init_imogen.o): $(ARC)(allocate_arrays.o) $(ARC)(switches.o) $(ARC)(file_utils.o) \
                       $(ARC)(inout.o) $(ARC)(time_loc.o) $(ARC)(time_mod.o) $(ARC)(spin_mod.o) \
                       $(ARC)(timeconst.o) $(ARC)(imogen_anlg_vals.o) $(ARC)(imogen_run.o) \
                       $(ARC)(imogen_time.o) $(ARC)(imogen_progs.o) $(ARC)(imogen_constants.o) \
                       $(ARC)(imogen_map.o) $(ARC)(imogen_clim.o) $(ARC)(imogen_io_vars.o) \
                       $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file. 
