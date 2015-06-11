
local=$(HOME)/local

HDF5_inc = /usr/include/hdf5/serial

CPPFLAGS=-I$(local)/include -I$(HDF5_inc)
CXXFLAGS=--std=c++11 -g
LDFLAGS=-L$(local)/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -g
LDLIBS=-lglsim_ol -lglsim -lhdf5 -lboost_serialization -lboost_system -lboost_program_options -Wl,-rpath -Wl,$(HOME)/local/lib

all: coulomb coulomb-ham

coulomb: coulomb.cc

coulomb-ham: coulomb-ham.cc

