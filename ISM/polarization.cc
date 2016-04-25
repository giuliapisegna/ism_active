/*
 * polarization.cc
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
    << "Reads the given glsim trajectory/configuration files and computes polarization.\n"
    << "\nOptions:\n";
  show_command_line_options(std::cerr);
  std::cerr  << "\n";
  show_parameters(std::cerr);
  std::cerr  << "\n";
}

/*****************************************************************************
 *
 * main
 *
 */

double polarization(glsim::OLconfiguration& conf)
{
  double V[3];
  memset(V,0,3*sizeof(double));
  for (int i=0; i<conf.N; ++i) {
    double v0=sqrt(modsq(conf.v[i]));
    V[0]+=conf.v[i][0]/v0;
    V[1]+=conf.v[i][1]/v0;
    V[2]+=conf.v[i][2]/v0;
  }
  return sqrt(modsq(V))/conf.N;
}

void wmain(int argc,char *argv[])
{
  CLoptions    opt;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(&conf);
  H5_multi_file ifs(options.ifiles,cfile);

  glsim::logs.set_stream(std::cout,glsim::error);
  ifs.read();

  std::cout << "#Step  Time  Polarization\n";
  do {
    double pol=polarization(conf);
    std::cout << conf.step << "  " << conf.time << "  " << pol << '\n';
  } while (ifs.read());
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
