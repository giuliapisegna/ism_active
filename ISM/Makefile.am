#
# Makefile.am for active/ISM
#

#
# binaries
#

EXTRA_DIST = social.hh isi.hh 3dvecs.hh vici.hh ovicsek.hh

bin_PROGRAMS = ism ism_create_conf vicsek ovicsek traj212col examine	\
               ck cknew ckiso ckt cktiso cko polarization tcorr1

ism_SOURCES = ism.cc isi.cc social.cc

ism_create_conf_SOURCES = ism_create_conf.cc

vicsek_SOURCES = vicsek.cc vici.cc social.cc

ovicsek_SOURCES = ovicsek.cc

traj212col_SOURCES = traj212col.cc

examine_SOURCES = examinetraj.cc

tcorr1_SOURCES = tcorr1.cc
tcorr1_LDADD = -lglsim_an $(AM_LDADD)

ck_SOURCES = ck.cc

cknew_SOURCES = cknew.cc
cknew_LDADD = -lglsim_an $(AM_LDADD)

ckiso_SOURCES = ckiso.cc
ckiso_LDADD = -lglsim_an $(AM_LDADD)
ckiso_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)

ckt_SOURCES = ckt.cc
ckt_LDADD = -lglsim_an $(AM_LDADD)

cktiso_SOURCES = cktiso.cc
cktiso_CXXFLAGS=$(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)

cko_SOURCES = cko.cc
cko_LDADD = -lglsim_an $(AM_LDADD)

#clambdat_SOURCES = clambdat.cc
#clambdat_LDFLAGS = -llapacke -lglsim_an

polarization_SOURCES = polarization.cc
