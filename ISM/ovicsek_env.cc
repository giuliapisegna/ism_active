/*
 * ovicsek_env.cc -- Environment for original Vicsek model
 *
 */

/*****************************************************************************
 *
 * VicsekEnvironment
 * 
 */

#include "ovicsek_env.hh"

OVicsekParameters::OVicsekParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("OVicsek.steps",po::value<int>()->required(),"Steps to run")
    ("OVicsek.fixed_graph",po::bool_switch()->required(),"False if birds are flying")
    ("OVicsek.v0",po::value<double>()->required(),"Speed")
    ("OVicsek.cutoff",po::value<double>()->required(),"Interaction cutoff")
    ("OVicsek.eta",po::value<double>()->required(),"Noise (range 0--1)")
    ;
}

OVicsekEnvironment::OVicsekEnvironment(const char* scope) :
  SimEnvironment(scope),
  VSsteps(0),
  fixed_graph(false),
  v0(1.),
  cutoff(1.),
  eta(0.1),
  total_number(0),
  polarization(0),
  SE(scope),
  par(scope)
{}

void OVicsekEnvironment::common_init()
{
  VSsteps=par.value("OVicsek.steps").as<int>();
  fixed_graph=par.value("OVicsek.fixed_graph").as<bool>();
  v0=par.value("OVicsek.v0").as<double>();
  cutoff=par.value("OVicsek.cutoff").as<double>();
  eta=par.value("OVicsek.eta").as<double>();
}

template <typename Archive>
inline void OVicsekEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("VicsekEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & VSsteps & fixed_graph & cutoff & eta;
}

