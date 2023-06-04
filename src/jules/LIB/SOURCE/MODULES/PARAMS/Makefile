include ../../../Makefile.common.mk
THISFILE = ./Makefile

OBJ = $(ARC)(blend_h.o)     $(ARC)(c_ch4.o)      $(ARC)(c_epslon.o)   $(ARC)(c_mdi.o)       \
      $(ARC)(c_r_cp.o)      $(ARC)(c_z0h_z0m.o)  $(ARC)(caoptr.o)     $(ARC)(c_charnk.o)    \
      $(ARC)(c_g.o)         $(ARC)(c_pdm.o)      $(ARC)(c_sicehc.o)   $(ARC)(descent.o)     \
      $(ARC)(sigm.o)        $(ARC)(ccarbon.o)    $(ARC)(c_densty.o)   $(ARC)(c_gamma.o)     \
      $(ARC)(c_perma.o)     $(ARC)(c_sulchm.o)   $(ARC)(fldtype.o)    $(ARC)(soil_param.o)  \
      $(ARC)(csigma.o)      $(ARC)(c_ht_m.o)     $(ARC)(c_pi.o)        \
      $(ARC)(c_surf.o)      $(ARC)(nvegparm.o)   $(ARC)(surf_param.o) $(ARC)(csmin.o)       \
      $(ARC)(c_kappai.o)   $(ARC)(c_rhowat.o)   $(ARC)(c_topog.o)     \
      $(ARC)(pftparm.o)     $(ARC)(veg_param.o)  $(ARC)(c_0_dg_c.o)   \
      $(ARC)(c_lheat.o)     $(ARC)(c_rough.o)    $(ARC)(c_vkman.o)    $(ARC)(rad_param.o)   \
      $(ARC)(seed.o)        $(ARC)(trif.o)       $(ARC)(nstypes.o)    $(ARC)(timeconst.o)   \
      $(ARC)(c_rmol.o)      $(ARC)(max_dimensions.o)                  $(ARC)(nvegparm_io.o) \
      $(ARC)(pftparm_io.o)  $(ARC)(snow_param.o) $(ARC)(trif_io.o)    $(ARC)(urban_param.o) \
      $(ARC)(atmos_constants_mod.o)              $(ARC)(conversions_mod.o)                  \
      $(ARC)(earth_constants_mod.o)              $(ARC)(water_constants_mod.o)              \
      $(ARC)(dust_param.o)                       $(ARC)(dust_parameters_mod.o)
all : $(OBJ) $(THISFILE)

##############################
## Dependencies             ##
##############################
$(ARC)(blend_h.o): $(THISFILE)
$(ARC)(c_0_dg_c.o): $(ARC)(conversions_mod.o) $(ARC)(water_constants_mod.o) $(THISFILE)
$(ARC)(c_ch4.o): $(THISFILE)
$(ARC)(c_charnk.o): $(THISFILE)
$(ARC)(c_densty.o): $(ARC)(water_constants_mod.o) $(THISFILE)
$(ARC)(c_epslon.o): $(ARC)(atmos_constants_mod.o) $(THISFILE)
$(ARC)(c_g.o): $(ARC)(earth_constants_mod.o) $(THISFILE)
$(ARC)(c_gamma.o): $(THISFILE)
$(ARC)(c_ht_m.o):  $(THISFILE)
$(ARC)(c_kappai.o):  $(THISFILE)
$(ARC)(c_kinds.o):  $(THISFILE)
$(ARC)(c_lheat.o): $(ARC)(water_constants_mod.o) $(THISFILE)
$(ARC)(c_mdi.o): $(THISFILE)
$(ARC)(c_pdm.o):  $(THISFILE)
$(ARC)(c_perma.o):  $(THISFILE)
$(ARC)(c_pi.o): $(ARC)(conversions_mod.o) $(THISFILE)
$(ARC)(c_r_cp.o):  $(ARC)(atmos_constants_mod.o) $(THISFILE)
$(ARC)(c_rhowat.o): $(ARC)(water_constants_mod.o) $(THISFILE)
$(ARC)(c_rmol.o):  $(THISFILE)
$(ARC)(c_rough.o): $(THISFILE)
$(ARC)(c_sicehc.o): $(THISFILE)
$(ARC)(c_sulchm.o): $(THISFILE)
$(ARC)(c_surf.o): $(THISFILE)
$(ARC)(c_topog.o): $(THISFILE)
$(ARC)(c_vkman.o): $(ARC)(atmos_constants_mod.o) $(THISFILE)
$(ARC)(c_z0h_z0m.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(caoptr.o): $(THISFILE)
$(ARC)(ccarbon.o): $(THISFILE)
$(ARC)(csigma.o): $(THISFILE)
$(ARC)(csmin.o): $(THISFILE)
$(ARC)(descent.o):  $(THISFILE)
$(ARC)(dust_param.o): $(ARC)(dust_parameters_mod.o) $(THISFILE)
$(ARC)(dust_parameters_mod.o):  $(THISFILE)
$(ARC)$(ARC)(fldtype.o):  $(THISFILE)
$(ARC)(max_dimensions.o): $(THISFILE)
$(ARC)(nstypes.o):  $(THISFILE)
$(ARC)(nvegparm.o):  $(THISFILE)
$(ARC)(nvegparm_io.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(pftparm.o): $(THISFILE)
$(ARC)(pftparm_io.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(rad_param.o): $(THISFILE)
$(ARC)(seed.o):  $(THISFILE)
$(ARC)(sigm.o):  $(THISFILE)
$(ARC)(snow_param.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(soil_param.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(surf_param.o): $(ARC)(c_epslon.o) $(ARC)(c_lheat.o) $(ARC)(c_r_cp.o) $(ARC)(c_g.o) $(THISFILE)
$(ARC)(timeconst.o) : $(THISFILE)
$(ARC)(trif.o): $(THISFILE)
$(ARC)(trif_io.o): $(ARC)(max_dimensions.o) $(THISFILE)
$(ARC)(urban_param.o): $(THISFILE)
$(ARC)(veg_param.o):  $(THISFILE)
$(ARC)(atmos_constants_mod.o): $(THISFILE)
$(ARC)(conversions_mod.o): $(THISFILE)
$(ARC)(earth_constants_mod.o): $(THISFILE)
$(ARC)(water_constants_mod.o): $(FILE) 

.PHONY : clean
clean :
	$(RM) $(wildcard *.o *$(MOD_FSUF))

## End of file.
