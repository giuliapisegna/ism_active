/*
 * ovicsek.hh -- Simulation class for Vicsek model (as originally formulated
 *
 */

#ifndef OVICSEK_HH
#define OVICSEK_HH

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


/*****************************************************************************
 *
 * OVicsekSimulation
 * 
 */

class OVicsekSimulation : public glsim::Simulation {
public:
  OVicsekSimulation(OVicsekEnvironment& e,glsim::OLconfiguration &c);
  ~OVicsekSimulation();
  const char* name() const {return "Vicsek's model (original)";}

  void step();
  void log(),log_start_sim();

protected:
  void update_observables();

  OVicsekEnvironment&         env;
  glsim::OLconfiguration& conf;

private:
  void vnoise(double*);
  void update_velocities();

  glsim::NeighbourList_subcells *NN;
  double  (*confb)[3];
  double  rcsq,v0sq;

  glsim::Uniform_real *ranz,*ranphi;
} ;


/******************************************************************************
 *
 * Observable
 *
 */

class OVicsekObservable_parameters : public glsim::Parameters {
public:
  OVicsekObservable_parameters(const char* scope);
} ;

class OVicsekObservable :  public glsim::SBObservable {
public:
  OVicsekObservable(OVicsekEnvironment&,glsim::OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

private:
  OVicsekEnvironment  &env;
  glsim::OLconfiguration &conf;
  OVicsekObservable_parameters par;

  void update();
} ;

inline OVicsekObservable::OVicsekObservable(OVicsekEnvironment& e,glsim::OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}


#endif /* OVICSEK_HH */
