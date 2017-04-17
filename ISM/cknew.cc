/******************************************************************************
 *
 * Options and main
 *
 */

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/Grk.hh"

#include "3dvecs.hh"

struct optlst {
public:
  int  nk,kdir;
  bool connect;
  std::vector<std::string> ifiles;
} options;


class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("ck")
{
  hidden_command_line_options().add_options()
    ("nk",po::value<int>(&options.nk)->required(),"number of wavevectors to compute")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("connect,C",po::bool_switch(&options.connect)->default_value(false),
     "Computed connected C(k)")
     ;

  positional_options().add("nk",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] nk ifile [ifile ....]\n\n"
    << "Computes C(k) for nk wavevectors kn (Delta k computed automatically).\n"
    << "Directions must be specified with option -d.\n"
    << "This computes\n\n"
    << "              C(k) = (1/N) \\sum_{ij} v_i v_i exp(i k r_{ij}) / v0^2\n\n"
    << "or, if -C (connect) is given, the connected C(k):\n\n"
    << "      C(k) = (1/N) \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i exp(i k r_{ij})\n\n"
    << "where \\delta v_i = v_i/v_0 - (1/N) \\sum_i v_i/v_0 and \\delta\\hat v_i is\n"
    << "normalized.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main and normalize
 *
 */

void normalize_vel(glsim::OLconfiguration& conf)
{
  /* Compute polarization */
  double V[3];
  memset(V,0,3*sizeof(double));
  for (int i=0; i<conf.N; ++i) {
    double v0=sqrt(modsq(conf.v[i]));
    conf.v[i][0]/=v0;
    conf.v[i][1]/=v0;
    conf.v[i][2]/=v0;
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
  }
  V[0]/=conf.N;
  V[1]/=conf.N;
  V[2]/=conf.N;

  /* Substract polarization */
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]-=V[0];
    conf.v[i][1]-=V[1];
    conf.v[i][2]-=V[2];
  }

  /* Normalize */
  double fn=0;
  for (int i=0; i<conf.N; ++i)
    fn+=modsq(conf.v[i]);
  fn= fn==0 ? 1 : sqrt(fn/conf.N);
  for (int i=0; i<conf.N; ++i) {
    conf.v[i][0]/=fn;
    conf.v[i][1]/=fn;
    conf.v[i][2]/=fn;
  }
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf,conft;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  ifs.read();
  glsim::Gk G(conf.box_length,options.nk);
  ifs.rewind();
  while (ifs.read()) {
    if (options.connect) normalize_vel(conf);
    G.push(conf,conf.v);
  }
  std::cout << G;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
