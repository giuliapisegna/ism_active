#
# Simulations of the ISM model
#

local=$(HOME)/local

HDF5_inc = /usr/include/hdf5/serial

# Debug
#CPPOPTFLAGS=-DDEBUG
#CXXOPTFLAGS=-g
#LDOPTFLAGS=-g
# Profile
CPPOPTFLAGS=
CXXOPTFLAGS=-g -pg
LDOPTFLAGS=-g -pg
#Release
#CPPOPTFLAGS=
#CXXOPTFLAGS=-O3
#LDOPTFLAGS=

CPPFLAGS=-I$(local)/include -I$(HDF5_inc) $(CPPOPTFLAGS)
CXXFLAGS=--std=c++11 $(CXXOPTFLAGS)
LDFLAGS=-L$(local)/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial $(LDOPTFLAGS)
LDLIBS=-lglsim_ol -lglsim -lhdf5 -lboost_serialization -lboost_timer -lboost_system -lboost_program_options -lgsl -Wl,-rpath -Wl,$(HOME)/local/lib

all: ism ism_create_conf vicsek traj212col

ism: ism.o isi.o social.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

ism.o: ism.cc isi.hh social.hh

vicsek: vicsek.o vici.o social.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

vici.o: vici.cc vici.hh social.hh

traj212col: traj212col.o
	g++ -o $@ $(LDFLAGS) $^ $(LDLIBS)

trak212col.o: traj212col.cc
