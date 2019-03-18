/*
 * ovicsek_env.hh -- Environment for original Vicsek model
 *
 */

#ifndef OVICSEK_ENV_HH
#define OVICSEK_ENV_HH

#include "glsim/random.hh"
#include "glsim/stochastic.hh"
#include "glsim/simulation.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/observable.hh"
#include "glsim/md.hh"
#include "glsim/avevar.hh"

/*****************************************************************************
 *
 * Environment
 * 
 */

class OVicsekParameters : public glsim::Parameters {
public:
  OVicsekParameters(const char* scope);
} ;

class OVicsekEnvironment : public glsim::SimEnvironment {
public:
  OVicsekEnvironment(const char* scope=glsim::Parameters::default_scope);

  long    VSsteps;
  bool    fixed_graph;
  double  cutoff,rcsq;
  double  eta;
  double  v0;

  // Observables
  int     total_number;
  double  polarization,v0sqave;
  double  Vcm[3];

protected:
  void    init_local() {SimEnvironment::init_local(); common_init();}
  void    warm_init_local() {SimEnvironment::warm_init_local(); common_init();}
  void    update_observables();

private:
  glsim::StochasticEnvironment SE;
  OVicsekParameters par;

  void common_init();
  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const unsigned int class_version=1;
} ;

BOOST_CLASS_VERSION(OVicsekEnvironment,OVicsekEnvironment::class_version);

#endif /* OVICSEK_ENV_HH */
