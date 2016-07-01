#
# Simulations of the ISM model
#

local=$(HOME)/local

HDF5_inc = /usr/include/hdf5/serial

# Debug
# CPPOPTFLAGS=-DDEBUG
# CXXOPTFLAGS=-g
# LDOPTFLAGS=-g
# Profile
#CPPOPTFLAGS=
#CXXOPTFLAGS=-g -pg
#LDOPTFLAGS=-g -pg
#Release
CPPOPTFLAGS=
CXXOPTFLAGS=-O3
LDOPTFLAGS=
# CPPOPTFLAGS=
# CXXOPTFLAGS=-O3 -fopenmp
# LDOPTFLAGS=-lgomp

CPPFLAGS=-I$(local)/include -I$(HDF5_inc) $(CPPOPTFLAGS)
CXXFLAGS=--std=c++11 $(CXXOPTFLAGS)
LDFLAGS=-L$(local)/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial $(LDOPTFLAGS)
LDLIBS=-lglsim_ol -lglsim -lhdf5 -lboost_serialization -lboost_timer -lboost_system -lboost_program_options -lgsl  -Wl,-rpath -Wl,$(HOME)/local/lib
#LDLIBS=-static -lglsim_ol -lglsim -lhdf5_hl -lhdf5 -lboost_serialization -lboost_timer -lboost_system -lboost_program_options -lgsl -lz -lboost_chrono -lpthread -ldl -Wl,-rpath -Wl,$(HOME)/local/lib

all: ism ism_create_conf vicsek traj212col examine tcorr1 ck ckt cko polarization

ism: ism.o isi.o social.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

ism.o: ism.cc isi.hh social.hh

vicsek: vicsek.o vici.o social.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

vici.o: vici.cc vici.hh social.hh

traj212col: traj212col.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

trak212col.o: traj212col.cc

examine: examinetraj.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

examinetraj.o: examinetraj.cc

tcorr1: tcorr1.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS) -lglsim_an

tcorr1.o: tcorr1.cc

ck: ck.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS) -lglsim_an

ck.o: ck.cc

ckt: ckt.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS) -lglsim_an

ckt.o: ckt.cc

cko: cko.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS) -lglsim_an

cko.o: cko.cc

clambdat.o: clambdat.cc

clambdat: clambdat.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS) -llapacke -lglsim_an


polarization: polarization.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS) -lglsim_an

polarization.o: polarization.cc
