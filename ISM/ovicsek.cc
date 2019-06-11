/*
 * ovicsek.cc -- Simulation of the original Vicsek model
 *
 */

#include "3dvecs.hh"
#include "ovicsek_env.hh"
#include "ovicsek.hh"
#include "glsim/offlattice.hh"

/*****************************************************************************
 *
 * OVicsekSimulation
 * 
 */

OVicsekSimulation::OVicsekSimulation(OVicsekEnvironment& e,
				     glsim::OLconfiguration &c) :
  Simulation(e,c),
  env(e),
  conf(c),
  confb(0)
{
  NN=new glsim::NeighbourList_subcells(env.cutoff);
  NN->rebuild(conf);
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

  ranz=new glsim::Uniform_real(1-2*env.eta,1);
  ranphi=new glsim::Uniform_real(0,2*M_PI);
  
  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();
}

OVicsekSimulation::~OVicsekSimulation()
{
  delete ranz,ranphi;
  delete[] confb;
}

// Modified sign function that returns 1 for val=0
template <typename T> T ssg(T val) {
    return (T(0) <= val) - (val < T(0));
}

void OVicsekSimulation::vnoise(double v[])
{
  double u1[3],u2[3];

  u1[0]=-ssg(v[0])*v[2];
  u1[1]=-ssg(v[1])*v[2];
  u1[2]=ssg(v[0])*v[0] + ssg(v[1])*v[1];
  normalize(v);
  normalize(u1);
  vprod(u2,v,u1);
  double z=(*ranz)();
  double phi=(*ranphi)();
  v[0]=env.v0*(sqrt(1-z*z)*cos(phi)*u1[0] + sqrt(1-z*z)*sin(phi)*u2[0] + z*v[0]);
  v[1]=env.v0*(sqrt(1-z*z)*cos(phi)*u1[1] + sqrt(1-z*z)*sin(phi)*u2[1] + z*v[1]);
  v[2]=env.v0*(sqrt(1-z*z)*cos(phi)*u1[2] + sqrt(1-z*z)*sin(phi)*u2[2] + z*v[2]);
}
  

void OVicsekSimulation::update_velocities()
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

void OVicsekSimulation::step()
{
  update_velocities();
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0] += conf.v[i][0];
    conf.r[i][1] += conf.v[i][1];
    conf.r[i][2] += conf.v[i][2];
  }

  conf.fold_coordinates();
  NN->update(env.v0);
  env.time_completed+=1;
  env.time_in_run+=1;
  env.run_completed = env.steps_in_run>=env.VSsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}

/*
 * Observables and log
 *
 */

void OVicsekSimulation::update_observables()
{
  double V[3];

  env.v0sqave=0;
  memset(V,0,3*sizeof(double));
    
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    env.v0sqave+=vs;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
  }
  env.v0sqave/=conf.N;
  env.polarization=sqrt(modsq(V))/(conf.N*env.v0);
}

void OVicsekSimulation::log_start_sim()
{
  char buff[300];
  
  Simulation::log_start_sim();
  glsim::logs(glsim::info) << "    Step       Time    <vsq>        Phi\n";

  double Ntimessteps=conf.N*env.steps_completed;
  sprintf(buff," Initial            %10.3e %10.3e\n",
	  env.v0sqave,env.polarization);
  glsim::logs(glsim::info) << buff;
}

void OVicsekSimulation::log()
{
  update_observables();
  static char buff[301];
  snprintf(buff,300,"%8ld %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	   env.v0sqave,env.polarization);
  glsim::logs(glsim::info) << buff;
}

/******************************************************************************
 *
 * Observable
 *
 */

OVicsekObservable_parameters::OVicsekObservable_parameters(const char* scope) :
  glsim::Parameters(scope)
{
  parameter_file_options().add_options()
    ("OVicsek.obs_interval",po::value<int>()->default_value(0),
     "Interval for standard observation, 0=skip")
    ("OVicsek.obs_file_prefix",po::value<std::string>(),"Observation file prefix")
    ;
}

void OVicsekObservable::interval_and_file()
{
  obs_interval=par.value("OVicsek.obs_interval").as<int>();
  obs_file_prefix=par.value("OVicsek.obs_file_prefix").as<std::string>();
}

void OVicsekObservable::write_header()
{
  fprintf(of,"#   (1)| |     (2)| |     (3)| |     (4)| |     (5)| |     (6)| |     (7)||\n");
  fprintf(of,"#- Step and time -| | Av v^2 | |--- Center of mass velocity --| |Polariz.|\n");
  fprintf(of,"#   Step       Time  <|v_i|^2>       VCMx       VCMy       VXMz        Phi\n");
}

void OVicsekObservable::observe()
{
  update();
  double Ntimessteps=conf.N*env.steps_completed;
  fprintf(of,"%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	  env.steps_completed,env.time_completed,
	  env.v0sqave,env.Vcm[0],env.Vcm[1],env.Vcm[2],
	  env.polarization);
}

void OVicsekObservable::update()
{
  double V[3];

  env.v0sqave=0;
  memset(V,0,3*sizeof(double));
  
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
}

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  OVicsekEnvironment  env;
  // OVicsekParameters   VP;
  glsim::OLconfiguration conf;
  OVicsekObservable obs(env,conf);
  glsim::Trajectory traj(env,conf);

  glsim::SimulationCL CL("ovicsek (Original Vicsek's model)","(C) 2015-2016 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  OVicsekSimulation sim(env,conf);

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
