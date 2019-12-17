#
# Makefile.am for active/ISM
#

#
# binaries
#

EXTRA_DIST = social.hh isi.hh 3dvecs.hh vici.hh ovicsek.hh

bin_PROGRAMS = ism ism-pl ism_create_conf vicsek ovicsek ovicsek_st traj212col examine	\
               ck cknew ckiso ckt cktiso cko polarization tcorr1 ism_potential

ism_SOURCES = ism.cc isi.cc social.cc

ism_pl_SOURCES = ism-pl.cc isi.cc social.cc

ism_create_conf_SOURCES = ism_create_conf.cc

ism_potential_SOURCES = ism_potential.cc isi_potential.cc social_potential.cc

vicsek_SOURCES = vicsek.cc vici.cc social.cc

ovicsek_SOURCES = ovicsek.cc ovicsek_env.cc

ovicsek_st_SOURCES = ovicsek_ST.cc ovicsek_env.cc

traj212col_SOURCES = traj212col.cc

examine_SOURCES = examinetraj.cc

tcorr1_SOURCES = tcorr1.cc

ck_SOURCES = ck.cc

cknew_SOURCES = cknew.cc

ckiso_SOURCES = ckiso.cc
ckiso_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)

ckt_SOURCES = ckt.cc

cktiso_SOURCES = cktiso.cc
cktiso_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)

cko_SOURCES = cko.cc

#clambdat_SOURCES = clambdat.cc
#clambdat_LDFLAGS = -llapacke -lglsim_an

polarization_SOURCES = polarization.cc
