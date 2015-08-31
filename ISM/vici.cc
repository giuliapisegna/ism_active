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
    ("VicsekI.eta",po::value<double>()->required(),"Friction coefficient")
    ("VicsekI.rescale_v0",po::bool_switch()->default_value(false),"If true, rescale speed to v0 before first step")
    ;
}

VicsekEnvironment::VicsekEnvironment(const char* scope) :
  SimEnvironment(scope),
  VSsteps(0),
  time_step(1e-5),
  fixed_graph(false),
  temperature(1.),
  eta(1.),
  total_number(0),
  total_social_mass(0),
  polarization(0),
  SE(scope),
  par(scope)
{}

void VicsekEnvironment::common_init()
{
  VSsteps=par.value("VicsekI.steps").as<int>();
  fixed_graph=par.value("VicsekI.fixed_graph").as<bool>();
  time_step=par.value("VicsekI.time_step").as<double>();
  temperature=par.value("VicsekI.temperature").as<double>();
  eta=par.value("VicsekI.eta").as<double>();
  rescale_v0=par.value("VicsekI.rescale_v0").as<bool>();
}

template <typename Archive>
inline void VicsekEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & VSsteps & fixed_graph & time_step & temperature & eta;
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
  if (env.rescale_v0) {
    for (int i=0; i<conf.N; ++i) {
      double v0=sqrt(modsq(conf.v[i]));
      conf.v[i][0]*=env.v0/v0;
      conf.v[i][1]*=env.v0/v0;
      conf.v[i][2]*=env.v0/v0;
    }
  }

  env.social_potential_energy=inter->social_potential_energy_and_acceleration(conf,conf.a);

  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();

  // Constants for Langevin integration
  Dt=env.time_step;
  xDt = env.fixed_graph ? 0 : Dt;
  mass=inter->social_mass(0);
  double etasv0=env.eta/v0sq;   // Because in the Vicsek friction is actually rotational friction
  double xi=etasv0/mass;
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
  if (etasv0<1e-3) {
    sv=env.temperature*Dt*Dt*Dt*etasv0*(2./3.-0.5*xidt)/(mass*mass);
  } else {
    sv=(env.temperature/etasv0)*(2*Dt-(3-4*exi+exi*exi)/xi);
  }
  sv=sqrt(sv);
  noise=new glsim::Gaussian_distribution(sigma,0);
}

VicsekSimulation::~VicsekSimulation()
{
  delete noise;
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
    deltav[0] = Dt*conf.a[i][0] + (*noise)();
    deltav[1] = Dt*conf.a[i][1] + (*noise)();
    deltav[2] = Dt*conf.a[i][2] + (*noise)();

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
    conf.v[i][0] = w*conf.v[i][0] + deltav[0];
    conf.v[i][1] = w*conf.v[i][1] + deltav[1];
    conf.v[i][2] = w*conf.v[i][2] + deltav[2];

  }

  // 5. Compute forces with new positions and velocities (in
  //    overdamped Langevin, velocity, rather than acceleration, is
  //    proportional to the force, here conf.a is the "velocity's velocity")
  env.social_potential_energy=inter->social_potential_energy_and_acceleration(conf,conf.a);

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
  double V[3],S[3];

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
    S[0]*=inter->social_mass(0);
    S[1]*=inter->social_mass(0);
    S[2]*=inter->social_mass(0);
    env.total_spin[0]+=S[0];
    env.total_spin[1]+=S[1];
    env.total_spin[2]+=S[2];
  }
  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.social_kinetic_energy*=0.5;
  env.social_total_energy=env.social_kinetic_energy+env.social_potential_energy;
  env.total_spinsq=modsq(env.total_spin);
}

void VicsekSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  glsim::logs(glsim::info) << "    Step       Time    SocEtot      <vsq>        Phi   <Stot^2>\n";

  sprintf(buff," Initial            %10.3e %10.3e %10.3e %10.3e\n",
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  env.total_spinsq);
  glsim::logs(glsim::info) << buff;
}

void VicsekSimulation::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  env.total_spinsq);
  glsim::logs(glsim::info) << buff;
}

/******************************************************************************
 *
 * Observable
 *
 */

ISMObservable_parameters::ISMObservable_parameters(const char* scope) :
  glsim::Parameters(scope)
{
  parameter_file_options().add_options()
    ("ISM.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("ISM.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

void ISMObservable::interval_and_file()
{
  obs_interval=par.value("ISM.obs_interval").as<int>();
  obs_file_prefix=par.value("ISM.obs_file_prefix").as<std::string>();
}

void ISMObservable::write_header()
{
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)| |     (8)| |     (9)| |    (10)| |    (11)| |    (12)| |    (13)| |    (14)| |    (15)| |    (16)| |    (17)|\n");
  fprintf(of,"#- Step and time -| |------- Social energy --------| | Av v^2 | |--- Center of mass velocity --| |Polariz.| |---------- Total spin --------|  |---------- Spin (single conf) ----------|\n");
  fprintf(of,"#   Step       Time  Potential    Kinetic      Total  <|v_i|^2>       VCMx       VCMy       VXMz        Phi         Sx         Sy         Sz    Average   Variance        Min        Max\n");

}

void ISMObservable::observe()
{
  update();
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,env.social_total_energy/conf.N,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization,
	  env.total_spin[0],env.total_spin[1],env.total_spin[2],
	  env.spinsqavar.ave(),env.spinsqavar.var(),env.spinsqavar.min(),env.spinsqavar.max());
}

void ISMObservable::update()
{
  double V[3],S[3];

  env.social_kinetic_energy=0;
  env.v0sqave=0;
  memset(env.total_spin,0,3*sizeof(double));
  memset(V,0,3*sizeof(double));
  env.spinsqavar.clear();
  
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
    env.social_kinetic_energy+=env.social_mass[conf.type[i]]*modsq(conf.a[i]);
    vprod(S,conf.v[i],conf.a[i]);
    S[0]*=env.social_mass[0];
    S[1]*=env.social_mass[0];
    S[2]*=env.social_mass[0];
    env.total_spin[0]+=S[0];
    env.total_spin[1]+=S[1];
    env.total_spin[2]+=S[2];
    env.spinsqavar.push(modsq(S));
  }
  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.Vcm[0]=V[0]/conf.N;
  env.Vcm[1]=V[1]/conf.N;
  env.Vcm[2]=V[2]/conf.N;
  env.social_kinetic_energy*=0.5;
  env.social_total_energy=env.social_kinetic_energy+env.social_potential_energy;
  env.total_spinsq=modsq(env.total_spin);
}
