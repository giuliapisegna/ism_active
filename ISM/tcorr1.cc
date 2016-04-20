/*
 * tcorr1.cc -- compute correlation def1
 *
 */

#define HAVE_LIBFFTW3 1

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/fft.hh"
#include "glsim/timecorr.hh"

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

CLoptions::CLoptions() : UtilityCL("tcorr1")
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
    << "Computes the time correlation, def 1\n"
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

std::vector<double>& operator+=(std::vector<double>& a,const std::vector<double>&b)
{
  for (int i=0; i<a.size(); ++i)
    a[i]+=b[i];
  return a;
}

void wmain(int argc,char *argv[])
{
  CLoptions    opt;
  opt.parse_command_line(argc,argv);

  OLconfiguration conf;
  OLconfig_file cfile(&conf);
  H5_multi_file ifs(options.ifiles,cfile);

  glsim::logs.set_stream(std::cout,glsim::error);

  // Compute the time series (vectorial polarization)

  std::vector<double> phix,phiy,phiz;

  double deltat;
  while (ifs.read()) {
    double phiv[3];
    phiv[0]=phiv[1]=phiv[2]=0;
    for (int i=0; i<conf.N; ++i) {
      double v0=sqrt(modsq(conf.v[i]));
      for (int j=0; j<3; ++j)
	phiv[j]+=conf.v[i][j]/v0;
    }
    phix.push_back(phiv[0]/conf.N);
    phiy.push_back(phiv[1]/conf.N);
    phiz.push_back(phiv[2]/conf.N);
    if (phix.size()==1)
      deltat=conf.time;
    else if (phix.size()==2)
      deltat=conf.time-deltat;
  }

  // Compute self-correlation
  typedef glsim::RealFFTW     rFFT;

  rFFT realft(glsim::FFT::in_place);
  std::vector<double> tcorr;
  correlation_connected_1d_tti_fft(phix,realft,phix.size(),false);
  tcorr=realft.tdata();
  correlation_connected_1d_tti_fft(phiy,realft,phiy.size(),false);
  tcorr+=realft.tdata();
  correlation_connected_1d_tti_fft(phiz,realft,phiz.size(),false);
  tcorr+=realft.tdata();
  /* Integral */
  double sum=0;
  for (int i=0; i<tcorr.size(); i++)
    sum+=tcorr[i];
  sum*=deltat;


  // Print
  printf("# Time   corr(def1)\n#\n");
  printf("#\n# Integral = %g\n#\n",sum);
  for (int i=0; i<tcorr.size(); i++)
      printf("%g %g\n",i*deltat,tcorr[i]);
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
