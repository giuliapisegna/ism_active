/*
 * ovicsek_ST.hh -- Simulation class for the original Vicsek model
 * with self tuning
 *
 */

#ifndef OVICSEK_ST_HH
#define OVICSEK_ST_HH

#include <queue>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/queue.hpp>

#include "glsim/random.hh"
#include "glsim/stochastic.hh"
#include "glsim/simulation.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/observable.hh"
#include "glsim/md.hh"
#include "glsim/avevar.hh"

#include "ovicsek_env.hh"

/*****************************************************************************
 *
 * Environment
 * 
 */

class OVicsek_STParameters : public glsim::Parameters {
public:
  OVicsek_STParameters(const char* scope);
} ;

class OVicsek_STEnvironment : public OVicsekEnvironment {
public:
  OVicsek_STEnvironment(const char* scope=glsim::Parameters::default_scope);

  bool    tune;
  bool    can_tune;
  int     tune_time;
  int     tune_step;
  short   dsign;
  long    last_tuning;
  double  tuned_eta;    // because OVicsek_Environment::eta is always re-read from .ini
  double  tune_factor;  // kappa in Dante's manuscript
  double  polarizationAveSQ,polarizationVar;
  std::queue<double>  polarization_prev;
  double  polprev;
  double  AC1;                           // Self-correlation at time 1 (running)
  double  AC1_prev;

protected:
  void    init_local() {OVicsekEnvironment::init_local(); tuned_eta=eta; common_init();}
  void    warm_init_local() {OVicsekEnvironment::warm_init_local(); eta=tuned_eta; common_init();}

private:
  OVicsek_STParameters par;

  void common_init();
  void vserial(oarchive_t &ar) {ar << *this;}
  void vserial(iarchive_t &ar) {ar >> *this;}
  template <typename Archive>
  void serialize(Archive &ar,const unsigned int version);
  friend class boost::serialization::access;

public:
  static const unsigned int class_version=2;
} ;

BOOST_CLASS_VERSION(OVicsek_STEnvironment,OVicsek_STEnvironment::class_version);


/*****************************************************************************
 *
 * OVicsek_STSimulation
 * 
 */

class OVicsek_STSimulation : public glsim::Simulation {
public:
  OVicsek_STSimulation(OVicsek_STEnvironment& e,glsim::OLconfiguration &c);
  ~OVicsek_STSimulation();
  const char* name() const {return "Original Vicsek's model with self-tuning)";}

  void step();
  void log(),log_start_sim();

protected:
  void update_observables();

  OVicsek_STEnvironment&         env;
  glsim::OLconfiguration& conf;

private:
  glsim::AveVar<false> polAV,polSQAV,AC1AV;

  void vnoise(double*);
  void update_velocities(),update_velocities_small_system();
  void tune_eta();

  glsim::NeighbourList_subcells *NN;
  double  (*confb)[3];
  double  rcsq,v0sq;

  glsim::Uniform_real *ranu,*ranphi;
} ;

/******************************************************************************
 *
 * Observable
 *
 */

class OVicsek_STObservable_parameters : public glsim::Parameters {
public:
  OVicsek_STObservable_parameters(const char* scope);
} ;

class OVicsek_STObservable :  public glsim::SBObservable {
public:
  OVicsek_STObservable(OVicsek_STEnvironment&,glsim::OLconfiguration&);

  void interval_and_file();
  void write_header();
  void observe();

private:
  OVicsek_STEnvironment  &env;
  glsim::OLconfiguration &conf;
  OVicsek_STObservable_parameters par;

  void update();
} ;

inline OVicsek_STObservable::OVicsek_STObservable(OVicsek_STEnvironment& e,glsim::OLconfiguration &c) :
  SBObservable(e),
  env(e),
  conf(c),
  par(e.scope())
{}

#endif /* OVICSEK_ST_HH */


