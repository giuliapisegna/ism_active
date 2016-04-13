/*
 * examinetraj.cc
 *
 */

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"

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

  do {
    printf("%7d anom %14.7e %14.7e %14.7e norm %14.7e %14.7e %14.7e\n",conf.step,
	   conf.v[55][0],conf.v[55][1],conf.v[55][2],
	   conf.v[10][0],conf.v[10][1],conf.v[10][2]);
    printf("        acce %14.7e %14.7e %14.7e norm %14.7e %14.7e %14.7e\n",conf.step,
	   conf.a[55][0],conf.a[55][1],conf.a[55][2],
	   conf.a[10][0],conf.a[10][1],conf.a[10][2]);
  } while (ifs.read());
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
