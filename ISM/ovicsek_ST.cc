/*
 * ovicsek_ST.cc -- Simulation of the original Vicsek model with self-tuning
 *
 */

#include "3dvecs.hh"
#include "ovicsek_ST.hh"
#include "glsim/offlattice.hh"
#include "glsim/avevar.hh"

/*****************************************************************************
 *
 * VicsekEnvironment
 * 
 */

OVicsek_STParameters::OVicsek_STParameters(const char *scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("OVicsek_ST.tune",po::value<bool>()->required(),"Whether to do the tuning")
    ("OVicsek_ST.tune_time",po::value<int>()->required(),"Time to use for tuning with autocorrelation, i.e. use AC(arg)")
    ("OVicsek_ST.tune_step",po::value<int>()->required(),"Steps between tuning (means arg values of AC will be averaged before tuning")
    ("OVicsek_ST.tune_factor",po::value<double>()->required(),"Tuning factor")
    ;
}

OVicsek_STEnvironment::OVicsek_STEnvironment(const char* scope) :
  OVicsekEnvironment(scope),
  tune(false),
  can_tune(false),
  tune_time(1),
  tune_step(100),
  dsign(1),
  last_tuning(0),
  tune_factor(0.1),
  polarizationAveSQ(0.),
  polarizationVar(0.),
  polarization_prev(),
  polprev(0.),
  AC1(0.),
  AC1_prev(0.),
  par(scope)
{}

void OVicsek_STEnvironment::common_init()
{
  tune=par.value("OVicsek_ST.tune").as<bool>();
  tune_time=par.value("OVicsek_ST.tune_time").as<int>();
  tune_step=par.value("OVicsek_ST.tune_step").as<int>();
  tune_factor=par.value("OVicsek_ST.tune_factor").as<double>();
}

template <typename Archive>
inline void OVicsek_STEnvironment::serialize(Archive &ar,const unsigned int version)
{
  if (version!=class_version && version!=1)
    throw glsim::Environment_wrong_version("Vicsek_STEnvironment",version,class_version);
  ar & boost::serialization::base_object<OVicsekEnvironment>(*this);
  if (version==1)  {
    ar & tune & tune_step & dsign & last_tuning & tuned_eta & tune_factor
      & polarizationAveSQ & polarizationVar & polprev & AC1 & AC1_prev;
  } else {
    ar & tune & can_tune & tune_time & tune_step & dsign & last_tuning & tuned_eta & tune_factor
      & polarizationAveSQ & polarizationVar & polarization_prev & polprev & AC1 & AC1_prev;
  }
}

/*****************************************************************************
 *
 * OVicsek_STSimulation
 * 
 */

OVicsek_STSimulation::OVicsek_STSimulation(OVicsek_STEnvironment& e,
					   glsim::OLconfiguration &c) :
  Simulation(e,c),
  env(e),
  conf(c),
  confb(0)
{
  NN=new glsim::NeighbourList_subcells(env.cutoff);
  try {
    NN->rebuild(conf);
  } catch (const glsim::System_too_small &e) {
    glsim::logs(glsim::warn) << "Small system, using naive nearest neighbours\n";
    delete NN;
    NN=0;
  }

  rcsq=env.cutoff*env.cutoff;
  env.total_number=conf.N;
  if (!conf.type) {
    conf.type=new short[conf.N];
    memset(conf.type,0,conf.N*sizeof(short));
  }
  v0sq=env.v0*env.v0;

  confb=new double[conf.N][3];
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
  } else if (env.initialization_kind()!=glsim::Environment::load) {
    for (int i=0; i<conf.N; ++i) {
      double v0=sqrt(modsq(conf.v[i]));
      conf.v[i][0]*=env.v0/v0;
      conf.v[i][1]*=env.v0/v0;
      conf.v[i][2]*=env.v0/v0;
    }
  }

  ranu=new glsim::Uniform_real(0.,1.);
  ranphi=new glsim::Uniform_real(0,2*M_PI);
  
  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();
}

OVicsek_STSimulation::~OVicsek_STSimulation()
{
  delete ranu;
  delete ranphi;
  delete[] confb;
  delete NN;
}

// Modified sign function that returns 1 for val=0
template <typename T> T ssg(T val) {
    return (T(0) <= val) - (val < T(0));
}

void OVicsek_STSimulation::vnoise(double v[])
{
  double u1[3],u2[3],a,range;

  u1[0]=-ssg(v[0])*v[2];
  u1[1]=-ssg(v[1])*v[2];
  u1[2]=ssg(v[0])*v[0] + ssg(v[1])*v[1];
  normalize(v);
  normalize(u1);
  vprod(u2,v,u1);
  a=1-2*env.eta;
  range=1-a;
  double z=a+range*(*ranu)();
  double phi=(*ranphi)();
  v[0]=env.v0*(sqrt(1-z*z)*cos(phi)*u1[0] + sqrt(1-z*z)*sin(phi)*u2[0] + z*v[0]);
  v[1]=env.v0*(sqrt(1-z*z)*cos(phi)*u1[1] + sqrt(1-z*z)*sin(phi)*u2[1] + z*v[1]);
  v[2]=env.v0*(sqrt(1-z*z)*cos(phi)*u1[2] + sqrt(1-z*z)*sin(phi)*u2[2] + z*v[2]);
}
  

void OVicsek_STSimulation::update_velocities()
{
  memset(confb,0,conf.N*3*sizeof(double));

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {

    double dsq=conf.distancesq(conf.r[p->first],conf.r[p->second]);
    if (dsq>rcsq) continue;

    confb[p->first][0] += conf.v[p->second][0];
    confb[p->first][1] += conf.v[p->second][1];
    confb[p->first][2] += conf.v[p->second][2];
    confb[p->second][0] += conf.v[p->first][0];
    confb[p->second][1] += conf.v[p->first][1];
    confb[p->second][2] += conf.v[p->first][2];

  }

  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]+=confb[i][0];
    conf.v[i][1]+=confb[i][1];
    conf.v[i][2]+=confb[i][2];
    vnoise(conf.v[i]); // Add noise to the resulting velocity
  }
}

void OVicsek_STSimulation::update_velocities_small_system()
{
  memset(confb,0,conf.N*3*sizeof(double));

  for (int i=0; i<conf.N-1; ++i)
    for (int j=i+1; j<conf.N; ++j) {
      double dsq=conf.distancesq(conf.r[i],conf.r[j]);
      if (dsq>rcsq) continue;

      confb[i][0] += conf.v[j][0];
      confb[i][1] += conf.v[j][1];
      confb[i][2] += conf.v[j][2];
      confb[j][0] += conf.v[i][0];
      confb[j][1] += conf.v[i][1];
      confb[j][2] += conf.v[i][2];
  }

  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]+=confb[i][0];
    conf.v[i][1]+=confb[i][1];
    conf.v[i][2]+=confb[i][2];
    vnoise(conf.v[i]); // Add noise to the resulting velocity
  }
}

void OVicsek_STSimulation::step()
{
  if (NN) update_velocities();
  else update_velocities_small_system();
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0] += conf.v[i][0];
    conf.r[i][1] += conf.v[i][1];
    conf.r[i][2] += conf.v[i][2];
  }

  conf.fold_coordinates();
  if (NN) NN->update(env.v0);
  update_observables();
  if (env.steps_completed % env.tune_step == 0) tune_eta();
  env.time_completed+=1;
  env.time_in_run+=1;
  env.run_completed = env.steps_in_run>=env.VSsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

void OVicsek_STSimulation::tune_eta()
{
  env.last_tuning=env.steps_completed;
  if (env.tune && env.AC1_prev!=0) {
    double s=ssg(env.AC1-env.AC1_prev);
    env.dsign*=s;
    double delta=1-env.AC1;
    delta*=delta;
    env.eta+=delta*env.dsign*env.tune_factor;
    if (env.eta<0) env.eta=1e-4;
    if (env.eta>1) env.eta=1;
    env.tuned_eta=env.eta;
  }
  env.AC1_prev=env.AC1;
  AC1AV.clear();
  polAV.clear();
  polSQAV.clear();
}

/*
 * Observables and log
 *
 */

void OVicsek_STSimulation::update_observables()
{
  double V[3];

  env.v0sqave=0;
  memset(V,0,3*sizeof(double));

  env.polarization_prev.push(env.polarization);
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
  }
  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
  env.Vcm[0]=V[0]/conf.N;
  env.Vcm[1]=V[1]/conf.N;
  env.Vcm[2]=V[2]/conf.N;

  polAV.push(env.polarization);
  env.polarizationAveSQ=polAV.ave()*polAV.ave();
  env.polarizationVar=polAV.var();
  if (env.polarization_prev.size()<env.tune_time) return;
  env.can_tune=true;
  env.polprev=env.polarization_prev.front();
  env.polarization_prev.pop();
  AC1AV.push(env.polarization*env.polprev);
  env.AC1=(AC1AV.ave()-env.polarizationAveSQ)/env.polarizationVar;
}

void OVicsek_STSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  glsim::logs(glsim::info) << "    Step       Time    <vsq>        Phi        eta\n";

  double Ntimessteps=conf.N*env.steps_completed;
  sprintf(buff," Initial            %10.3e %10.3e\n",
	  env.v0sqave,env.polarization);
  glsim::logs(glsim::info) << buff;
}

void OVicsek_STSimulation::log()
{
  static char buff[301];
  snprintf(buff,300,"%8ld %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	   env.v0sqave,env.polarization,env.eta);
  glsim::logs(glsim::info) << buff;
}

/******************************************************************************
 *
 * Observable
 *
 */

OVicsek_STObservable_parameters::OVicsek_STObservable_parameters(const char* scope) :
  glsim::Parameters(scope)
{
  parameter_file_options().add_options()
    ("OVicsek_ST.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("OVicsek_ST.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

void OVicsek_STObservable::interval_and_file()
{
  obs_interval=par.value("OVicsek_ST.obs_interval").as<int>();
  obs_file_prefix=par.value("OVicsek_ST.obs_file_prefix").as<std::string>();
}

void OVicsek_STObservable::write_header()
{
  fprintf(of,"#   Tuning every %d steps using AC(n) with n=%d\n",env.tune_step,env.tune_time);
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)|       (8)|       (9)|      (10)|     (11)|      (12)|      (13)|\n");
  fprintf(of,"#- Step and time -| | Av v^2 | |--- Center of mass velocity --| |---------------- Polarizazion  ------------------- |          | Elapsed  |\n");
  fprintf(of,"#   Step       Time  <|v_i|^2>       VCMx       VCMy       VXMz        Phi PreviousPh    PhiAv^2     PhiVar     AC(n)        eta since tune\n");
}

void OVicsek_STObservable::observe()
{
  update();
  double Ntimessteps=conf.N*env.steps_completed;
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10d\n",
	  env.steps_completed,env.time_completed,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization,env.polprev,env.polarizationAveSQ,env.polarizationVar,env.AC1,env.eta,
	  env.steps_completed-env.last_tuning);
}

void OVicsek_STObservable::update()
{
}

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  OVicsek_STEnvironment  env;
  glsim::OLconfiguration conf;
  OVicsek_STObservable obs(env,conf);
  glsim::Trajectory traj(env,conf);

  glsim::SimulationCL CL("ovicsek_ST (Original Vicsek's model with self-tuning)","(C) 2015-2019 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  OVicsek_STSimulation sim(env,conf);

  traj.observe_first();
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
