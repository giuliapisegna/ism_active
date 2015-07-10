/*
 * isi.cc -- The inertial spin integrator
 *
 * Integrates the equation of motion for an inertial spin model (with
 * Brownian noise), without specifiying the velocity-velocity
 * interaction.
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#

#include "3dvecs.hh"
#include "isi.hh"

/*****************************************************************************
 *
 * ISMEnvironment
 * 
 */

ISMParameters::ISMParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("ISM.steps",po::value<int>()->required(),"Steps to run")
    ("ISM.time_step",po::value<int>()->required(),"Delta t")
    ("ISM.fixed_graph",po::bool_switch()->required(),"False if birds are flying")
    ("ISM.v0",po::value<double>()->required(),"Speed")
    ("ISM.temperature",po::value<double>()->required(),"Temperature for dv/dt friction")
    ("ISM.eta",po::value<double>()->required(),"Friction coefficient")
    ;
}

ISMEnvironment::ISMEnvironment(const char* scope) :
  SimEnvironment(scope),
  ISsteps(0),
  time_step(1e-5),
  fixed_graph(false),
  v0(1.),
  temperature(1.),
  eta(1.),
  total_number(0),
  total_social_mass(0),
  polarization(0),
  SE(scope),
  par(scope)
{}

void ISMEnvironment::common_init()
{
  ISsteps=par.value("ISM.steps").as<int>();
  fixed_graph=par.value("ISM.fixed_graph").as<bool>();
  time_step=par.value("ISM.time_step").as<double>();
  v0=par.value("ISM.v0").as<double>();
  temperature=par.value("ISM.temperature").as<double>();
  eta=par.value("ISM.eta").as<double>();
}

template <typename Archive>
inline void ISMEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & ISsteps & fixed_graph & time_step & v0 & temperature & eta;
}

/*****************************************************************************
 *
 * ISMSimulation
 * 
 */

ISMSimulation::ISMSimulation(ISMEnvironment& e,glsim::OLconfiguration &c,SocialInteractions *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i),
  confb(0)
{
  env.total_number=conf.N;
  env.total_social_mass=0;
  for (int i=0; i<conf.N; ++i) {
    env.total_social_mass+=inter->social_mass(conf.type[i]);
    // env.Ptot[0]+=conf.v[i][0]*inter->mass(conf.type[i]);
    // env.Ptot[1]+=conf.v[i][1]*inter->mass(conf.type[i]);
    // env.Ptot[2]+=conf.v[i][2]*inter->mass(conf.type[i]);
  }

  if (conf.a==0) {
    conf.a=new double[conf.N][3];
    memset(conf.a,0,conf.N*3*sizeof(double));
  }
  confb=new double[conf.N][3];
  env.social_potential_energy=inter->social_potential_energy_and_acceleration(conf,confb);
  // env.Epot=inter->potential_energy(conf);

  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  //update_observables();

  // Constants for Langevin integration
  Dt=env.time_step;
  double xDt = env.fixed_graph ? 0 : Dt;
  double mass=inter->social_mass(0);
  double xi=env.eta/mass;
  double xidt=xi*Dt;
  double c0l,c1,c2,sx,sv,rho;
  double exi=exp(-xidt);
  if (xidt<1e-3) {
    c0l=1 - xidt + xidt*xidt/2 - xidt*xidt*xidt/6;
    c1=1 - xidt/2 + xidt*xidt/6 - xidt*xidt*xidt/24;
    c2=0.5 - xidt/6 + xidt*xidt/24;
    rho=sqrt(3.)*(0.5-xidt/16.-(17./1280.)*xidt*xidt
		  +(17./6144)*xidt*xidt*xidt);
  } else {
    c0l=exi;
    c1=(1-c0l)/xidt;
    c2=(1-c1)/xidt;
    rho=(1-exi)*(1-exi)/sqrt( (1-exi*exi)*(3*xidt-3+4*exi-exi*exi) );
  }
  sa=(env.temperature/mass)*(1-exi*exi);
  sa=sqrt(sa);
  if (env.eta<1e-3) {
    sv=env.temperature*Dt*Dt*Dt*env.eta*(2./3.-0.5*xidt)/(mass*mass);
  } else {
    sv=(env.temperature/env.eta)*(2*Dt-(3-4*exi+exi*exi)/xi);
  }
  sv=sqrt(sv);
  noise=new glsim::BivariateGaussian_distribution(sv,sa,rho);
  
  c0=c0l;
  c1dt=Dt*c1;
  c2dt=Dt*c2;
  c1mc2=c1-c2;

  v0sq=0;
}

ISMSimulation::~ISMSimulation()
{
  delete noise;
}

template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void ISMSimulation::step()
{
  double w,xiv[3],xia[3],deltav[3];
  
  for (int i=0; i<conf.N; ++i) {

    // Get noise*3
    (*noise)(xiv[0],xia[0]);
    (*noise)(xiv[1],xia[1]);
    (*noise)(xiv[2],xia[2]);

    // 1. Update x
    conf.r[i][0] += conf.v[i][0]*xDt;
    conf.r[i][1] += conf.v[i][1]*xDt;
    conf.r[i][2] += conf.v[i][2]*xDt;

    // 2. Delta v w/o rattle corrections
    deltav[0] = c1dt*conf.a[i][0] + c2dt*Dt*confb[i][0] + xiv[0];
    deltav[1] = c1dt*conf.a[i][1] + c2dt*Dt*confb[i][1] + xiv[0];
    deltav[2] = c1dt*conf.a[i][2] + c2dt*Dt*confb[i][2] + xiv[0];

    // 3. Solve for lagrange multiplier (first part of rattle)
    double a = v0sq;
    double b = 2*dotp(conf.v[i],deltav);
    double c = modsq(deltav)-v0sq;
    if (b==0)
      w = sqrt(-c/a);
    else {
      double q = -0.5 * (b+sgn(b)*sqrt(b*b-4*a*c));
      w = b>0 ? c/q : q/a;
    }
    double hlambda = (w -1)/c2dt;

    // 4. Partial update of a (missing second rattle correction and
    // acceleration at next step)
    conf.a[i][0] = c0*conf.a[i][0] + c1mc2*(Dt*confb[i][0] + hlambda*conf.v[i][0]) + xia[0];
    conf.a[i][1] = c0*conf.a[i][1] + c1mc2*(Dt*confb[i][1] + hlambda*conf.v[i][1]) + xia[1];
    conf.a[i][2] = c0*conf.a[i][2] + c1mc2*(Dt*confb[i][2] + hlambda*conf.v[i][2]) + xia[2];

    // 5. Update v
    conf.v[i][0] = (1 + c2dt*hlambda)*conf.v[i][0] + deltav[0];
    conf.v[i][1] = (1 + c2dt*hlambda)*conf.v[i][1] + deltav[1];
    conf.v[i][2] = (1 + c2dt*hlambda)*conf.v[i][2] + deltav[2];

  }

  // 6. Compute forces with new positions and velocities
  env.social_potential_energy=inter->social_potential_energy_and_acceleration(conf,confb);

  for (int i=0; i<conf.N; ++i) {

    // 7. Further update of a
    conf.a[i][0] += c2dt*confb[i][0];
    conf.a[i][1] += c2dt*confb[i][1];
    conf.a[i][2] += c2dt*confb[i][2];

    // 8. Solve for second lagrange multiplier (last step of rattle)
    double c2hgamma = - dotp(conf.v[i],conf.a[i])/v0sq;

    // 9. Final correction for a
    conf.a[i][0] += c2hgamma * conf.v[i][0];
    conf.a[i][1] += c2hgamma * conf.v[i][1];
    conf.a[i][2] += c2hgamma * conf.v[i][2];

  }
  inter->fold_coordinates(conf);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.ISsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

/*
 * Observables and log
 *
 */

void ISMSimulation::update_observables()
{
  double V[3];
  env.social_kinetic_energy=0;
  
  for (int i=0; i<conf.N; ++i) {
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
    env.social_kinetic_energy+=inter->social_mass(conf.type[i])*modsq(conf.a[i]);
  }
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.social_kinetic_energy*=0.5;
}

void ISMSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  glsim::logs(glsim::info) << "    Step       Time    SocEpot    SocEkin    SocEtot         Phi\n";

  sprintf(buff," Initial            %10.3e %10.3e %10.3e %10.3e\n",
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,
	  env.social_total_energy/conf.N,env.polarization/conf.N);
  glsim::logs(glsim::info) << buff;
}

void ISMSimulation::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,
	  env.social_total_energy/conf.N,env.polarization/conf.N);
  glsim::logs(glsim::info) << buff;
}

