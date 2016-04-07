/*
 * vici.cc -- Vicsek model integrator
 *
 * Integrates the equation of motion for the Vicsek model, understood
 * as the overdamped inertial spin model (with Brownian noise),
 * without specifiying the velocity-velocity interaction.  This means
 * that the model is not integrated in the usual way but with
 * overdamped Brownian dynamics ("position Langevin") plus a
 * constraint.
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, September 2015
 *
 */

#include "3dvecs.hh"
#include "vici.hh"

/*****************************************************************************
 *
 * VicsekEnvironment
 * 
 */

VicsekIntegratorParameters::VicsekIntegratorParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("VicsekI.steps",po::value<int>()->required(),"Steps to run")
    ("VicsekI.time_step",po::value<double>()->required(),"Delta t")
    ("VicsekI.fixed_graph",po::bool_switch()->required(),"False if birds are flying")
    ("VicsekI.temperature",po::value<double>()->required(),"Temperature for dv/dt friction")
    ("VicsekI.rescale_v0",po::bool_switch()->default_value(false),"If true, rescale speed to v0 before first step")
    ("VicsekI.planar_noise",po::bool_switch()->default_value(false),"If true, the random force will lie in the XY plane.  With appropriate initial conditions, this can be used to simulate the 2-D Vicsek model")
    ;
}

VicsekEnvironment::VicsekEnvironment(const char* scope) :
  SimEnvironment(scope),
  VSsteps(0),
  time_step(1e-5),
  fixed_graph(false),
  temperature(1.),
  total_number(0),
  total_social_mass(0),
  polarization(0),
  constraint_fails(0),
  SE(scope),
  par(scope)
{}

void VicsekEnvironment::common_init()
{
  VSsteps=par.value("VicsekI.steps").as<int>();
  fixed_graph=par.value("VicsekI.fixed_graph").as<bool>();
  time_step=par.value("VicsekI.time_step").as<double>();
  temperature=par.value("VicsekI.temperature").as<double>();
  rescale_v0=par.value("VicsekI.rescale_v0").as<bool>();
  planar_noise=par.value("VicsekI.planar_noise").as<bool>();
  constraint_fails=0;
}

template <typename Archive>
inline void VicsekEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("VicsekEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & VSsteps & fixed_graph & time_step & temperature & rescale_v0 & planar_noise;
  ar & constraint_fails;
}

/*****************************************************************************
 *
 * VicsekSimulation
 * 
 */

VicsekSimulation::VicsekSimulation(VicsekEnvironment& e,
				   glsim::OLconfiguration &c,
				   VicsekInteraction *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i),
  confb(0)
{
  env.total_number=conf.N;
  env.total_social_mass=0;
  if (!conf.type) {
    conf.type=new short[conf.N];
    memset(conf.type,0,conf.N*sizeof(short));
  }
  env.v0=inter->speed();
  v0sq=env.v0*env.v0;
  env.social_mass[0]=inter->social_mass(0);
  for (int i=0; i<conf.N; ++i)
    env.total_social_mass+=inter->social_mass(conf.type[i]);

  if (conf.v==0) {
    conf.v=new double[conf.N][3];
    glsim::Spherical3d_distribution sr;
    double u[3];
    for (int i=0; i<conf.N; ++i) {
      sr(u);
      conf.v[i][0]=env.v0*u[0];
      conf.v[i][1]=env.v0*u[1];
      conf.v[i][2]=env.v0*u[2];
    }
    glsim::logs(glsim::info) << "Initialized velocities to random directions\n";
  }

  if (conf.a==0) {
    conf.a=new double[conf.N][3];
    memset(conf.a,0,conf.N*3*sizeof(double));
  }

  // If asked, rescale v0
  if (env.initialization_kind()!=glsim::Environment::load && env.rescale_v0) {
    for (int i=0; i<conf.N; ++i) {
      double v0=sqrt(modsq(conf.v[i]));
      conf.v[i][0]*=env.v0/v0;
      conf.v[i][1]*=env.v0/v0;
      conf.v[i][2]*=env.v0/v0;
    }
  }

  env.social_potential_energy=inter->social_potential_energy_and_force(conf,conf.a);

  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();

  // Constants for Langevin integration
  Dt=env.time_step;
  xDt = env.fixed_graph ? 0 : Dt;
  double etasv0=1./v0sq;   // Because in the Vicsek friction is actually rotational friction
                           // and here eta=1 (overdamped, time rescaling)
  double sigma=sqrt(2*env.temperature*etasv0*Dt);
  noisexy=new glsim::Gaussian_distribution(sigma,0);
  if (env.planar_noise)
    noisez=new glsim::Gaussian_distribution(0,0);
  else
    noisez=noisexy;
}

VicsekSimulation::~VicsekSimulation()
{
  delete noisexy;
  if (env.planar_noise) delete noisez;
}

template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void VicsekSimulation::step()
{
  double w,deltav[3];
  
  for (int i=0; i<conf.N; ++i) {

    // 1. Update x
    conf.r[i][0] += conf.v[i][0]*xDt;
    conf.r[i][1] += conf.v[i][1]*xDt;
    conf.r[i][2] += conf.v[i][2]*xDt;

    // 2. Delta v w/o shake corrections
    deltav[0] = v0sq*(Dt*conf.a[i][0] + (*noisexy)());
    deltav[1] = v0sq*(Dt*conf.a[i][1] + (*noisexy)());
    deltav[2] = v0sq*(Dt*conf.a[i][2] + (*noisez)());

    // 3. Solve for lagrange multiplier
    double a = v0sq;
    double b = 2*dotp(conf.v[i],deltav);
    double c = modsq(deltav)-v0sq;
    if (b==0)
      w = sqrt(-c/a);
    else {
      double q = -0.5 * (b+sgn(b)*sqrt(b*b-4*a*c));
      w = b>0 ? c/q : q/a;
    }

    // 4. Update v
    if (isnan(w)) {
      // Problem: constraint cannot be enforced with a lagrange multiplier (centripetal force)
      // Do it the Vicsek way (just renormalizing v)

      conf.v[i][0]+=deltav[0];
      conf.v[i][1]+=deltav[1];
      conf.v[i][2]+=deltav[2];
      w=sqrt(v0sq/modsq(conf.v[i]));
      conf.v[i][0]*=w;
      conf.v[i][1]*=w;
      conf.v[i][2]*=w;

      env.constraint_fails++;

    } else {

      conf.v[i][0] = w*conf.v[i][0] + deltav[0];
      conf.v[i][1] = w*conf.v[i][1] + deltav[1];
      conf.v[i][2] = w*conf.v[i][2] + deltav[2];

    }

  }

  // 5. Compute forces with new positions and velocities (in
  //    overdamped Langevin, velocity, rather than acceleration, is
  //    proportional to the force, here conf.a is the "velocity's velocity")
  env.social_potential_energy=inter->social_potential_energy_and_force(conf,conf.a);

  inter->fold_coordinates(conf);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.VSsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

/*
 * Observables and log
 *
 */

void VicsekSimulation::update_observables()
{
  double V[3];

  env.social_kinetic_energy=0;
  env.v0sqave=0;
  memset(V,0,3*sizeof(double));
    
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
    env.social_kinetic_energy+=inter->social_mass(conf.type[i])*modsq(conf.a[i]);
  }
  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.social_kinetic_energy*=0.5;
  env.social_total_energy=env.social_kinetic_energy+env.social_potential_energy;
}

void VicsekSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  glsim::logs(glsim::info) << "    Step       Time    SocEtot      <vsq>        Phi  Fails/N/step\n";

  double Ntimessteps=conf.N*env.steps_completed;
  sprintf(buff," Initial            %10.3e %10.3e %10.3e %12.5e\n",
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  env.constraint_fails/Ntimessteps);
  glsim::logs(glsim::info) << buff;
}

void VicsekSimulation::log()
{
  update_observables();
  static char buff[301];
  snprintf(buff,300,"%8ld %10.3e %10.3e %10.3e %10.3e %12.5e\n",
	  env.steps_completed,env.time_completed,
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  (double)env.constraint_fails/(conf.N*env.steps_completed));
  glsim::logs(glsim::info) << buff;
}

/******************************************************************************
 *
 * Observable
 *
 */

VicsekObservable_parameters::VicsekObservable_parameters(const char* scope) :
  glsim::Parameters(scope)
{
  parameter_file_options().add_options()
    ("VicsekI.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("VicsekI.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

void VicsekObservable::interval_and_file()
{
  obs_interval=par.value("VicsekI.obs_interval").as<int>();
  obs_file_prefix=par.value("VicsekI.obs_file_prefix").as<std::string>();
}

void VicsekObservable::write_header()
{
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)| |     (8)| |     (9)| |    (10)| |       (11)|\n");
  fprintf(of,"#- Step and time -| |------- Social energy --------| | Av v^2 | |--- Center of mass velocity --| |Polariz.| |  Constraint\n");
  fprintf(of,"#   Step       Time  Potential    Kinetic      Total  <|v_i|^2>       VCMx       VCMy       VXMz        Phi  fails/N/step\n");
}

void VicsekObservable::observe()
{
  update();
  double Ntimessteps=conf.N*env.steps_completed;
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %12.5e\n",
	  env.steps_completed,env.time_completed,
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,env.social_total_energy/conf.N,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization,env.constraint_fails/Ntimessteps);
}

void VicsekObservable::update()
{
  double V[3];

  env.social_kinetic_energy=0;
  env.v0sqave=0;
  memset(V,0,3*sizeof(double));
  
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
    env.social_kinetic_energy+=env.social_mass[conf.type[i]]*modsq(conf.a[i]);
  }
  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.Vcm[0]=V[0]/conf.N;
  env.Vcm[1]=V[1]/conf.N;
  env.Vcm[2]=V[2]/conf.N;
  env.social_kinetic_energy*=0.5;
  env.social_total_energy=env.social_kinetic_energy+env.social_potential_energy;
}
