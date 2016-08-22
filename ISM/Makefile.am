#
# Makefile.am for active/ISM
#

#
# binaries
#

bin_PROGRAMS = ism ism_create_conf vicsek traj212col examine ck ckiso ckt	\
               cktiso cko polarization tcorr1

ism_SOURCES = ism.cc isi.cc social.cc

ism_create_conf_SOURCES = ism_create_conf.cc

vicsek_SOURCES = vicsek.cc vici.cc social.cc

traj212col_SOURCES = traj212col.cc

examine_SOURCES = examinetraj.cc

tcorr1_SOURCES = tcorr1.cc
tcorr1_LDADD = -lglsim_an $(AM_LDADD)

ck_SOURCES = ck.cc

ckiso_SOURCES = ckiso.cc

ckt_SOURCES = ckt.cc
ckt_LDADD = -lglsim_an $(AM_LDADD)

cktiso_SOURCES = cktiso.cc

cko_SOURCES = cko.cc
cko_LDADD = -lglsim_an $(AM_LDADD)

#clambdat_SOURCES = clambdat.cc
#clambdat_LDFLAGS = -llapacke -lglsim_an

polarization_SOURCES = polarization.cc
