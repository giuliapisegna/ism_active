/*
 * traj212col.cc -- convert glsim's trajectory files to COBBS 12 column
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

CLoptions::CLoptions() : UtilityCL("tarj212col")
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
    << "Reads the given glsim trajectory/configuration files and writes them in ascii\n"
    << "COBBS 12-column format.  Writes to stdout.\n"
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

  for (int npart=0; npart<conf.N; ++npart) {
    do {
      printf("TRK_%04d %8ld %14.7e %14.7e %14.7e ",
	     npart,conf.step,conf.r[npart][0],conf.r[npart][1],conf.r[npart][2]);
      if (conf.v)
	printf("%14.7e %14.7e %14.7e ",
	       conf.v[npart][0],conf.v[npart][1],conf.v[npart][2]);
      else
	printf("%14.7e %14.7e %14.7e ",0.,0.,0.);
      if (conf.a)
	printf("%14.7e %14.7e %14.7e 0\n",
	       conf.a[npart][0],conf.a[npart][1],conf.a[npart][2]);
      else
	printf("%14.7e %14.7e %14.7e 0\n",0.,0.,0.);
    } while (ifs.read());
    printf("\n");
    ifs.rewind();
    ifs.read();
  }
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
