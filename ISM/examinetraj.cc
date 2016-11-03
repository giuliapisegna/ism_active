/*
 * examinetraj.cc
 *
 */

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"

#include "3dvecs.hh"

using namespace glsim;

/*****************************************************************************
 *
 * options
 *
 */

static struct ooptions {
  std::vector<std::string>   ifiles;

} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : UtilityCL("examinetraj")
{
  hidden_command_line_options().add_options()
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "\nusage: " << progname << " [options] ifile [..]\n\n"
    << "Reads the given glsim trajectory/configuration files does things\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
  std::cerr  << "\n";
  show_parameters(std::cerr);
  std::cerr  << "\n";
}

void polarization(glsim::OLconfiguration& conf,double &v0sqave,double Vcm[],double &P)
{
  double V[3];

  v0sqave=0;
  memset(V,0,3*sizeof(double));
  memset(Vcm,0,3*sizeof(double));
  
  for (int i=0; i<conf.N; ++i) {
    double vs=modsq(conf.v[i]);
    v0sqave+=vs;
    vs=sqrt(vs);
    V[0]+=conf.v[i][0]/vs;
    V[1]+=conf.v[i][1]/vs;
    V[2]+=conf.v[i][2]/vs;
    Vcm[0]+=conf.v[i][0];
    Vcm[1]+=conf.v[i][1];
    Vcm[2]+=conf.v[i][2];
  }
  v0sqave/=conf.N;
  P=sqrt(modsq(V))/conf.N;
  Vcm[0]/=conf.N;
  Vcm[1]/=conf.N;
  Vcm[2]/=conf.N;
}

/*****************************************************************************
 *
 * main
 *
 */

void wmain(int argc,char *argv[])
{
  CLoptions    opt;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(&conf);
  H5_multi_file ifs(options.ifiles,cfile);

  glsim::logs.set_stream(std::cout,glsim::error);
  ifs.read();

  double v0sqave,Vcm[3],P;
  do {
    polarization(conf,v0sqave,Vcm,P);
    printf("%8ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	   conf.step,conf.time,
	   v0sqave,Vcm[0],Vcm[1],Vcm[2],
	   P);
  } while (ifs.read());
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
