
local=$(HOME)/local

HDF5_inc = /usr/include/hdf5/serial

CPPFLAGS=-I$(local)/include -I$(HDF5_inc)
CXXFLAGS=--std=c++11 -g
LDFLAGS=-L$(local)/lib -g
LDLIBS=-lglsim_ol -lglsim -lboost_serialization -lboost_system

all: coulomb

coulomb: coulomb.cc


