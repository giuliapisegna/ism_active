/*
 * isi_potential.cc -- The inertial spin integrator with external positional potential
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

#include <cmath>

#include "3dvecs.hh"
#include "isi_potential.hh"

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
    ("ISM.time_step",po::value<double>()->required(),"Delta t")
    ("ISM.fixed_graph",po::bool_switch()->required(),"False if birds are flying")
    ("ISM.temperature",po::value<double>()->required(),"Temperature for dv/dt friction")
    ("ISM.eta",po::value<double>()->required(),"Friction coefficient")
    ("ISM.rescale_v0",po::bool_switch()->default_value(false),"If true, rescale speed to v0 before first step")
    ("ISM.A", po::value<double>()-> required(),"Amplitude of oscillation of the potential")
    ("ISM.omega", po::value<double>()-> required(), "Frequency of oscillation of the potential")
    ;
}

ISMEnvironment::ISMEnvironment(const char* scope) :
  SimEnvironment(scope),
  ISsteps(0),
  time_step(1e-5),
  fixed_graph(false),
  temperature(1.),
  eta(1.),
  A(0),
  omega(0),
  total_number(0),
  total_social_mass(0),
  polarization(0),
  drl_x(0),
  SE(scope),
  par(scope)
{ rl[0]=0;
  rl[1]=0;
  rl[2]=0;
  Rcm[0]=0;
  Rcm[1]=0;
  Rcm[2]=0;
   }

void ISMEnvironment::common_init()
{
  ISsteps=par.value("ISM.steps").as<int>();
  fixed_graph=par.value("ISM.fixed_graph").as<bool>();
  time_step=par.value("ISM.time_step").as<double>();
  temperature=par.value("ISM.temperature").as<double>();
  eta=par.value("ISM.eta").as<double>();
  rescale_v0=par.value("ISM.rescale_v0").as<bool>();
  A=par.value("ISM.A").as<double>();
  omega=par.value("ISM.omega").as<double>();
}

template <typename Archive>
inline void ISMEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version)
    throw glsim::Environment_wrong_version("MDEnvironment",version,class_version);
  ar & boost::serialization::base_object<SimEnvironment>(*this);
  ar & ISsteps & fixed_graph & time_step & temperature & eta & A & omega;
}

/*****************************************************************************
 *
 * ISMSimulation
 * 
 */

ISMSimulation::ISMSimulation(ISMEnvironment& e,glsim::OLconfiguration &c,VicsekInteraction *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i),
  confb(0)
{

//proviamo a inizializzare il landmark sul centro di massa del sistema!!
/*
double R[3]={0.};
for(int i=0; i < conf.N; i++){
	R[0] += conf.r[i][0];
	R[1] += conf.r[i][1];
	R[2] += conf.r[i][2];
}

	R[0] /= conf.N; 
	R[1] /= conf.N; 
	R[2] /= conf.N; 

	env.rl[0] =  R[0]; 
	env.rl[1] =  R[1]; 
	env.rl[2] =  R[2]; 
        env.R0 =  R[0]; 
/*
//pongo il landmark al centro della scatola*/

  env.rl[0] = conf.box_length[0]/2;
  env.rl[1] = conf.box_length[1]/2;
  env.rl[2] = conf.box_length[2]/2;
  env.drl_x = env.A*env.omega;


//dobbiamo inizializzare posizioni e accelerazioni, mi serve chi e k che sono nei vicsek parameters, per ora li metto anche nell'environment, poi vediamo meglio li chiamo k2 e chi2
  env.v0=inter->speed();
  double k = inter -> stiffness(); 
//fprintf(stderr, "k %lf \n", k);
  double m = inter -> social_mass(0);

/*
  for(int i=0; i < conf.N; i++){
//fprintf(stderr, "vo %lf \n", env.v0);
    conf.a[i][0] = -(k/m) *(conf.r[i][0] - env.rl[0]);
    conf.a[i][1] = -(k/m) *(conf.r[i][1] - env.rl[1]);
    conf.a[i][2] = -(k/m) *(conf.r[i][2] - env.rl[2]);
  }

//scelgo un vettore a caso
  double n[3]={1,1,1};
  //lo normalizzo
  n[0] /= sqrt( n[0]*n[0] +n[1]*n[1] +n[2]*n[2] );
  n[1] /= sqrt( n[0]*n[0] +n[1]*n[1] +n[2]*n[2] );
  n[2] /= sqrt( n[0]*n[0] +n[1]*n[1] +n[2]*n[2] );

// faccio si' che la velocita'sia perpendicolare alla sua derivata = perpendicolare ad a = perpendicolare ad r)
 for(int i=0; i < conf.N; i++) {
 
  vprod(conf.v[i], n, conf.a[i]);
//fprintf(stderr, "velocita %lf  %lf %lf\n", conf.v[i][0] , conf.v[i][1] , conf.v[i][2] );
//fprintf(stderr, "accc %lf  %lf %lf\n", conf.a[i][0] , conf.a[i][1] , conf.a[i][2] );
  conf.v[i][0] /= sqrt( conf.v[i][0]*conf.v[i][0]  + conf.v[i][1]*conf.v[i][1]  + conf.v[i][2]*conf.v[i][2]);
  conf.v[i][1] /= sqrt( conf.v[i][0]*conf.v[i][0]  + conf.v[i][1]*conf.v[i][1]  + conf.v[i][2]*conf.v[i][2]);
  conf.v[i][2] /= sqrt( conf.v[i][0]*conf.v[i][0]  + conf.v[i][1]*conf.v[i][1]  + conf.v[i][2]*conf.v[i][2]);
 }


//fine della inizializzazione di accelerazione e velocita'*/ 

/*
//UN'ALTRA INIZIALIZZAZIONE DELLE Ai


double p[3];

for(int i = 0; i < conf.N; i++){
	
	vprod(p,conf.v[i],conf.r[i]);
        vprod(conf.a[i],p,conf.v[i]);
	conf.a[i][0] = conf.a[i][0]*k/(env.v0*env.v0);
        conf.a[i][1] = conf.a[i][1]*k/(env.v0*env.v0);
	conf.a[i][2] = conf.a[i][2]*k/(env.v0*env.v0);
}*/

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

  // If asked, rescale v0 conserving angular velocity
  if (env.initialization_kind()!=glsim::Environment::load && env.rescale_v0) {
    for (int i=0; i<conf.N; ++i) {
      double v0=sqrt(modsq(conf.v[i]));
      conf.v[i][0]*=env.v0/v0;
      conf.v[i][1]*=env.v0/v0;
      conf.v[i][2]*=env.v0/v0;
      conf.a[i][0]*=env.v0/v0;
      conf.a[i][1]*=env.v0/v0;
      conf.a[i][2]*=env.v0/v0;
    }
  }

  confb=new double[conf.N][3];
  env.social_potential_energy=inter->social_potential_energy_and_acceleration(conf,confb, env.rl, env.drl_x);

  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();

  // Constants for Langevin integration
  Dt=env.time_step;
  xDt = env.fixed_graph ? 0 : Dt;
  mass=inter->social_mass(0);
  double etasv0=env.eta/v0sq;   // Because in the ISM friction is actually rotational friction
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
  if (sa<=0) throw glsim::Runtime_error("Error computing variance of noise (got negative)",HERE);
  sa=sqrt(sa);
  if (etasv0<1e-3) {
    sv=env.temperature*Dt*Dt*Dt*etasv0*(2./3.-0.5*xidt)/(mass*mass);
  } else {
    sv=(env.temperature/etasv0)*(2*Dt-(3-4*exi+exi*exi)/xi);
  }
  if (sv<=0) throw glsim::Runtime_error("Error computing variance of noise (got negative)",HERE);
  sv=sqrt(sv);
  noise=new glsim::BivariateGaussian_distribution(sv,sa,rho);
  
  c0=c0l;
  c1dt=Dt*c1;
  c2dt=Dt*c2;
  c1mc2=c1-c2;
}

ISMSimulation::~ISMSimulation()
{
//I_AM_HERE;
  delete noise;
//I_AM_HERE;
}

template <typename T> T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void ISMSimulation::step()
{
  double w,xiv[3],xia[3],deltav[3];

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.ISsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
 
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
    deltav[1] = c1dt*conf.a[i][1] + c2dt*Dt*confb[i][1] + xiv[1];
    deltav[2] = c1dt*conf.a[i][2] + c2dt*Dt*confb[i][2] + xiv[2];

    // 3. Solve for lagrange multiplier (first part of rattle)
    double a = v0sq;
    double b = 2*dotp(conf.v[i],deltav);
    double c = modsq(deltav)-v0sq;
    if (b==0)
      w = sqrt(-c/a);
    else {
      double q = -0.5 * (b+sgn(b)*sqrt(b*b-4*a*c));
      w = b>0 ? c/q : q/a;
      if (std::isnan(w)) {
	glsim::logs(glsim::error) << "RATTLE failed: cannot enforce constraint\n";
	env.run_completed=true;
	return;
      }
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

  //5.bis update della coordinata del landmark

    env.rl[0] = conf.box_length[0]/2 + env.A*sin(env.omega*(env.time_completed));
    env.drl_x = env.A*env.omega*cos(env.omega*(env.time_completed));

  // 6. Compute forces with new positions and velocities
  env.social_potential_energy=inter->social_potential_energy_and_acceleration(conf,confb, env.rl, env.drl_x);

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


    double s[3];
    vprod(s, conf.v[i], conf.a[i]);
  //  fprintf(stderr," sx %lf sy %lf sz %lf \n", s[0],s[1],s[2]);

  }


  inter->rebuildcell(conf);

}

/*
 * Observables and log
 *
 */

void ISMSimulation::update_observables()
{

  double V[3],S[3],R[3];

  env.social_kinetic_energy=0;
  env.v0sqave=0;
  memset(env.total_spin,0,3*sizeof(double));
  memset(V,0,3*sizeof(double));
  memset(R,0,3*sizeof(double));
    
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
    R[0] += conf.r[i][0];
    R[1] += conf.r[i][1];
    R[2] += conf.r[i][2];
    env.social_kinetic_energy+=inter->social_mass(conf.type[i])*modsq(conf.a[i]);
    vprod(S,conf.v[i],conf.a[i]);
    S[0]*=inter->social_mass(0);
    S[1]*=inter->social_mass(0);
    S[2]*=inter->social_mass(0);
    env.total_spin[0]+=S[0];
    env.total_spin[1]+=S[1];
    env.total_spin[2]+=S[2];
  }

  env.Rcm[0]=R[0]/conf.N; 
  env.Rcm[1]=R[1]/conf.N; 
  env.Rcm[2]=R[2]/conf.N; 

  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.social_kinetic_energy*=0.5;
  env.social_total_energy=env.social_kinetic_energy+env.social_potential_energy;
  env.total_spinsq=modsq(env.total_spin);
}

void ISMSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();

  glsim::logs(glsim::info) << "    Step       Time    SocEtot      <vsq>        Phi   <Stot^2>   Kin	Pot	RCMx\n";

/*
  sprintf(buff," Initial            %10.3e %10.3e %10.3e %10.3e\n",
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  env.total_spinsq);*/


  sprintf(buff," Initial            %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  env.total_spinsq, env.social_kinetic_energy, env.social_potential_energy, env.Rcm[0]);

  glsim::logs(glsim::info) << buff;
}

void ISMSimulation::log()
{
  update_observables();
  static char buff[300];
  sprintf(buff,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e  %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.social_total_energy/conf.N,env.v0sqave,env.polarization,
	  env.total_spinsq, env.social_kinetic_energy, env.social_potential_energy, env.Rcm[0]);
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

/*
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)| |     (8)| |     (9)| |     (10)| |    (11)| |    (12)| |    (13)| |    (14)| |    (15)| |    (16)| |    (17)||    (18)| |    (19)| |    (20)|\n");
  fprintf(of,"#- Step and time -| |------- Social energy --------| | Av v^2 | |--- Center of mass velocity --| | Polariz.| |---------- Total spin --------|  |---------- Spin (single conf) ----------||--- Center of mass position --|\n");
  fprintf(of,"#   Step       Time  Potential    Kinetic      Total  <|v_i|^2>       VCMx       VCMy       VXMz         Phi         Sx         Sy         Sz    Average   Variance        Min        Max    	RCMx	  RCMy	    RCMz \n"); */


  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)| |     (8)| |     (9)| |     (10)| |    (11)| |    (12)| |    (13)| |    (14)| |    (15)| |    (16)| |    (17)||    (18)| |    (19)| |    (20)|\n");
  fprintf(of,"#- Step and time -| |------- Social energy --------| | Av v^2 | |--- Center of mass velocity --| | Polariz.| |---------- Landmark position --------|  |---------- Spin (single conf) ----------||--- Center of mass position --|\n");
  fprintf(of,"#   Step       Time  Potential    Kinetic      Total  <|v_i|^2>       VCMx       VCMy       VXMz         Phi         Rlx         Rly         Rlz    Average   Variance        Min        Max    	RCMx	  RCMy	    RCMz \n"); 


}

void ISMObservable::observe()
{
  update();
/*
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %11.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
	  env.steps_completed,env.time_completed,
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,env.social_total_energy/conf.N,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization,
	  env.total_spin[0],env.total_spin[1],env.total_spin[2],
	  env.spinsqavar.ave(),env.spinsqavar.var(),env.spinsqavar.min(),env.spinsqavar.max(), env.Rcm[0], env.Rcm[1], env.Rcm[2]);*/

  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %11.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e \n",
	  env.steps_completed,env.time_completed,
	  env.social_potential_energy/conf.N,env.social_kinetic_energy/conf.N,env.social_total_energy/conf.N,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization,
	  env.rl[0],env.rl[1],env.rl[2],
	  env.spinsqavar.ave(),env.spinsqavar.var(),env.spinsqavar.min(),env.spinsqavar.max(), env.Rcm[0], env.Rcm[1], env.Rcm[2]);


}

void ISMObservable::update()
{
  double V[3],S[3], R[3];

  env.social_kinetic_energy=0;
  env.v0sqave=0;
  memset(env.total_spin,0,3*sizeof(double));
  memset(V,0,3*sizeof(double));
  memset(R,0,3*sizeof(double));
  env.spinsqavar.clear();
  
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
    R[0] += conf.r[i][0];
    R[1] += conf.r[i][1];
    R[2] += conf.r[i][2];
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
  env.Rcm[0]=R[0]/conf.N; 
  env.Rcm[1]=R[1]/conf.N; 
  env.Rcm[2]=R[2]/conf.N; 
  env.social_kinetic_energy*=0.5;
  env.social_total_energy=env.social_kinetic_energy+env.social_potential_energy;
  env.total_spinsq=modsq(env.total_spin);
}
