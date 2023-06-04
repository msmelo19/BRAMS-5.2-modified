include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(cancap.o)     $(ARC)(fcdch.o)      $(ARC)(microbe.o)    $(ARC)(qsat_wat.o)    $(ARC)(sf_aero.o)    $(ARC)(sf_flux.o)    $(ARC)(sf_orog.o)    $(ARC)(sice_htf.o)   $(ARC)(vgrav.o) \
      $(ARC)(dewpnt.o)     $(ARC)(im_sf_pt.o)   $(ARC)(phi_m_h.o)    $(ARC)(phi_m_h_vol.o) $(ARC)(raero.o)      $(ARC)(sf_evap.o)    $(ARC)(sf_impl.o)    $(ARC)(sf_resist.o)  $(ARC)(smc_ext.o) \
      $(ARC)(dustresb.o)   $(ARC)(leaf.o)       $(ARC)(physiol.o)    $(ARC)(root_frac.o)   $(ARC)(sf_expl.o)    $(ARC)(sf_melt.o)    $(ARC)(sf_stom.o)    $(ARC)(stdev1.o) \
      $(ARC)(elevate.o)    $(ARC)(leaf_lit.o)   $(ARC)(qsat.o)       $(ARC)(screen_tq.o)   $(ARC)(sf_exch.o)    $(ARC)(sfl_int.o)    $(ARC)(sf_rib.o)     $(ARC)(soil_evap.o) \
      $(ARC)(ex_sf_ex.o)   $(ARC)(leaf_limits.o)  $(ARC)(sf_orog_gb.o)  $(ARC)(qsat_mix.o) $(ARC)(generate_anthrop_heat.o)  $(ARC)(urbanz0.o)  $(ARC)(get_us.o)

all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(cancap.o): $(ARC)(pftparm.o) $(ARC)(surf_param.o) $(ARC)(trif.o) $(THISFILE)
$(ARC)(dewpnt.o): $(ARC)(c_epslon.o) $(ARC)(c_r_cp.o) $(ARC)(c_lheat.o) \
                  $(ARC)(c_0_dg_c.o) $(THISFILE)
$(ARC)(dustresb.o): $(ARC)(c_sulchm.o) $(ARC)(dust_param.o) $(ARC)(dust_param.o)  \
                    $(ARC)(c_r_cp.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_pi.o)  \
                    $(ARC)(c_g.o) $(ARC)(dust_param.o) $(THISFILE)
$(ARC)(elevate.o): $(ARC)(c_r_cp.o) $(ARC)(c_lheat.o) $(ARC)(c_epslon.o) \
                   $(ARC)(c_g.o) $(THISFILE)
$(ARC)(ex_sf_ex.o): $(ARC)(dust_param.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_r_cp.o) \
                    $(ARC)(c_g.o) $(ARC)(csigma.o) $(ARC)(soil_param.o)       \
                    $(ARC)(nstypes.o) $(ARC)(fldtype.o) $(ARC)(veg_param.o)   \
                    $(THISFILE)
$(ARC)(fcdch.o): $(ARC)(c_vkman.o) $(ARC)(surf_param.o) $(ARC)(blopt8a.o) $(THISFILE)
$(ARC)(generate_anthrop_heat.o): $(ARC)(nstypes.o) $(ARC)(urban_param.o) \
                                 $(ARC)(switches_urban.o) $(ARC)(switches.o) \
                                 $(THISFILE)
$(ARC)(im_sf_pt.o): $(ARC)(c_r_cp.o) $(ARC)(c_lheat.o) $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(leaf.o): $(ARC)(pftparm.o) $(ARC)(surf_param.o) $(ARC)(c_rmol.o) $(THISFILE)
$(ARC)(leaf_limits.o): $(ARC)(c_0_dg_c.o) $(ARC)(pftparm.o) $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(leaf_lit.o): $(ARC)(pftparm.o) $(THISFILE)
$(ARC)(microbe.o): $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(phi_m_h.o): $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(phi_m_h_vol.o): $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(physiol.o): $(ARC)(c_densty.o) $(ARC)(nstypes.o) $(ARC)(nvegparm.o)  \
                   $(ARC)(soil_param.o) $(ARC)(pftparm.o) $(ARC)(surf_param.o) \
                   $(ARC)(urban_param.o) $(ARC)(switches_urban.o) \
                   $(ARC)(ancil_info.o) $(THISFILE)
$(ARC)(qsat.o): $(ARC)(c_0_dg_c.o) $(ARC)(c_epslon.o) $(THISFILE)
$(ARC)(qsat_mix.o): $(ARC)(c_0_dg_c.o) $(ARC)(c_epslon.o) $(THISFILE)
$(ARC)(qsat_wat.o): $(ARC)(c_epslon.o) $(ARC)(c_0_dg_c.o) $(THISFILE)
$(ARC)(raero.o): $(ARC)(c_vkman.o) $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(root_frac.o): $(THISFILE)
$(ARC)(screen_tq.o): $(ARC)(c_ht_m.o) $(ARC)(surf_param.o) $(THISFILE)
$(ARC)(sf_aero.o): $(ARC)(dust_param.o) $(THISFILE)
$(ARC)(sf_evap.o): $(ARC)(c_r_cp.o) $(ARC)(c_lheat.o) $(ARC)(c_0_dg_c.o) $(THISFILE)
$(ARC)(sf_exch.o): $(ARC)(c_z0h_z0m.o) $(ARC)(c_vkman.o) $(ARC)(c_rough.o)  \
                   $(ARC)(c_densty.o) $(ARC)(blend_h.o) $(ARC)(c_perma.o)  \
                   $(ARC)(c_r_cp.o) $(ARC)(c_g.o) $(ARC)(c_0_dg_c.o)       \
                   $(ARC)(dust_param.o) $(ARC)(c_kappai.o) $(ARC)(c_lheat.o) \
                   $(ARC)(blopt8a.o) $(ARC)(nstypes.o) $(ARC)(snow_param.o) \
                   $(ARC)(surf_param.o) $(ARC)(jules_mod.o) $(THISFILE)
$(ARC)(sf_expl.o): $(ARC)(dust_param.o) $(ARC)(c_0_dg_c.o) $(ARC)(c_r_cp.o) \
                   $(ARC)(c_g.o) $(ARC)(csigma.o) $(ARC)(c_mdi.o) $(ARC)(soil_param.o) \
                   $(ARC)(snow_param.o) $(ARC)(nstypes.o) $(ARC)(fldtype.o)  \
                   $(ARC)(veg_param.o) $(ARC)(switches.o) $(ARC)(ancil_info.o)  \
                   $(ARC)(prognostics.o) $(ARC)(c_elevate.o) $(ARC)(blopt8a.o) \
                   $(ARC)(solinc_data.o) $(THISFILE)
$(ARC)(sf_flux.o): $(ARC)(csigma.o) $(ARC)(c_r_cp.o) $(ARC)(c_lheat.o) $(ARC)(c_g.o) \
                   $(ARC)(c_0_dg_c.o) $(ARC)(c_epslon.o) $(ARC)(snow_param.o)   \
                   $(ARC)(surf_param.o) $(ARC)(urban_param.o) \
                   $(ARC)(nstypes.o) $(ARC)(switches_urban.o) \
                   $(ARC)(switches.o) $(THISFILE)
$(ARC)(sf_impl.o): $(ARC)(c_mdi.o) $(ARC)(caoptr.o) $(ARC)(csigma.o)   \
                   $(ARC)(c_lheat.o) $(ARC)(c_kappai.o) $(ARC)(c_r_cp.o) \
                   $(ARC)(c_0_dg_c.o) $(ARC)(surf_param.o) $(ARC)(ancil_info.o) \
                   $(ARC)(prognostics.o) $(ARC)(solinc_data.o) \
                   $(ARC)(jules_mod.o) $(THISFILE)
$(ARC)(sf_melt.o): $(ARC)(c_r_cp.o) $(ARC)(c_lheat.o) $(ARC)(c_0_dg_c.o)  \
                   $(ARC)(switches.o) $(ARC)(rad_param.o) $(THISFILE)
$(ARC)(sf_orog.o): $(ARC)(c_surf.o) $(ARC)(c_vkman.o) $(ARC)(c_mdi.o)  \
                   $(ARC)(surf_param.o) $(ARC)(blopt8a.o) $(THISFILE)
$(ARC)(sf_orog_gb.o): $(ARC)(c_surf.o) $(ARC)(c_vkman.o) $(ARC)(c_mdi.o)  \
                      $(ARC)(surf_param.o) $(ARC)(blopt8a.o) $(THISFILE)
$(ARC)(sf_resist.o): $(ARC)(switches.o) $(ARC)(rad_param.o) $(THISFILE)
$(ARC)(sf_rib.o): $(ARC)(c_r_cp.o) $(ARC)(c_g.o) $(THISFILE)
$(ARC)(sf_stom.o): $(ARC)(ccarbon.o) $(ARC)(pftparm.o) $(ARC)(surf_param.o)  \
                   $(ARC)(c_rmol.o) $(THISFILE)
$(ARC)(sfl_int.o): $(ARC)(c_vkman.o)     $(ARC)(surf_param.o) \
                   $(ARC)(urban_param.o) $(ARC)(switches_urban.o) \
                   $(ARC)(nstypes.o)     $(THISFILE)
$(ARC)(sice_htf.o): $(ARC)(c_mdi.o) $(ARC)(c_sicehc.o) $(ARC)(c_kappai.o)  \
                    $(ARC)(c_0_dg_c.o) $(THISFILE)
$(ARC)(smc_ext.o): $(THISFILE)
$(ARC)(soil_evap.o): $(THISFILE)
$(ARC)(stdev1.o): $(ARC)(c_g.o) $(THISFILE)
$(ARC)(vgrav.o): $(ARC)(dust_param.o) $(ARC)(c_sulchm.o) $(ARC)(c_0_dg_c.o)  \
                 $(ARC)(c_g.o) $(THISFILE)
$(ARC)(urbanz0.o) : $(ARC)(urban_param.o) $(ARC)(um_types.o) $(ARC)(nstypes.o) \
                    $(THISFILE)
$(ARC)(get_us.o)  : $(ARC)(urban_param.o) $(THISFILE)

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.

