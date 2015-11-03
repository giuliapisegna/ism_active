/*
 * vicsek.cc -- The Vicsek model
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, September 2015
 *
 */

#include "glsim/histogram.hh"
#include "vici.hh"


/*****************************************************************************/

class Energy_histogram_par :  public glsim::Parameters {
public:
  Energy_histogram_par(const char *scope);
} ;

Energy_histogram_par::Energy_histogram_par(const char *scope)
{
  parameter_file_options().add_options()
    ("Ehistogram.file_prefix",po::value<std::string>()->required(),"file")
    ("Ehistogram.Nbins",po::value<int>()->default_value(100),"Number of bins")
    ("Ehistogram.emin",po::value<double>()->required(),"Minimum energy to record")
    ("Ehistogram.emax",po::value<double>()->required(),"Maximum energy to record")
    ;
}

class Energy_histogram :  public glsim::SBObservable {
public:
  Energy_histogram(VicsekEnvironment&);
  ~Energy_histogram();
  void interval_and_file();
  void write_header() {}
  void observe();

private:
  Energy_histogram_par par;
  VicsekEnvironment    &env;
  glsim::Histogram     *histo;
} ;

Energy_histogram::Energy_histogram(VicsekEnvironment& e) :
  SBObservable(e),
  par(e.scope()),
  env(e),
  histo(0)
{}

void Energy_histogram::interval_and_file()
{
  obs_interval=1;
  obs_file_prefix=par.value("Ehistogram.file_prefix").as<std::string>();
  histo=new glsim::Histogram(par.value("Ehistogram.Nbins").as<int>(),
				       par.value("Ehistogram.emin").as<double>(),
				       par.value("Ehistogram.emax").as<double>());
}

void Energy_histogram::observe()
{
  histo->push(env.social_potential_energy);
}

Energy_histogram::~Energy_histogram()
{
  if (!of) return;
  fprintf(of,"# Energy histogram (%d outliers)\n",histo->outliers());
  std::ostringstream o;
  histo->probability_output();
  o << *histo;
  fputs(o.str().c_str(),of);
  delete histo;
}

/*****************************************************************************/


void wmain(int argc, char *argv[])
{
  VicsekEnvironment  env;
  VicsekParameters VP;
  glsim::OLconfiguration conf;
  VicsekObservable obs(env,conf);
  Energy_histogram ehis(env);
  // glsim::Trajectory traj(env,conf,
  // 			 glsim::OLconfig_file::options().r_frame());

  glsim::SimulationCL CL("vicsek (Vicsek's model)","(C) 2015 Tomas S. Grigera",env.scope());
  CL.parse_command_line(argc,argv);
  glsim::prepare(CL,env,conf);

  MetricVicsekInteraction inter(VP,conf);
  VicsekSimulation sim(env,conf,&inter);

// traj.observe_first();
  ehis.observe_first();
  obs.observe_first();
  sim.run();
  env.save();
  conf.save(env.configuration_file_fin);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
