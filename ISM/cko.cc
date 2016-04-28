/*
 * cko.cc -- Def2 in k,omega space
 *
 *
 */

#define HAVE_LIBFFTW3

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/timecorr.hh"
#include "glsim/fft.hh"

using glsim::vcomplex;
using glsim::dcomplex;

#ifdef HAVE_LIBFFTW3
typedef glsim::RealFFTW     rFFT;
typedef glsim::ComplexFFTW  cFFT;
#else
typedef glsim::RealFFT_gsl_2n     rFFT;
typedef glsim::ComplexFFT_gsl_2n  cFFT;
#endif /* HAVE_LIBFFTW3 */


/*****************************************************************************
 *
 *  Class for C(k).  I'll choose one direction of k, parallel to the
 *  cooordinate axes (set on construction), this way PBCs are
 *  correctly enforced.
 *
 */

class Cko {
public:
  Cko(double box_length[],double deltat_,int kn_,int kdir_);

  Cko& push_config(glsim::OLconfiguration &conf);
  Cko& compute_Cko();
  const vcomplex& Cko_data() const {return Cko_;} 

private:
  int                               Npart,kn,kdir;
  double                            k[3],deltak_[3];
  double                            deltat,delta_omega;
  vcomplex                          jkx_,jky_,jkz_;
  vcomplex                          Cko_;
  
  vcomplex j_k(glsim::OLconfiguration &);

  friend std::ostream& operator<<(std::ostream&,const Cko&); 
} ;

Cko::Cko(double box_length[],double deltat_,int kn_,int kdir_) :
  deltat(deltat_),
  kn(kn_), kdir(kdir_),
  jkx_(3),jky_(3),jkz_(3)
{
  // choose deltak
  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];

  k[0]=k[1]=k[2]=0.;
  k[kdir]=deltak_[kdir]*kn;
}

vcomplex Cko::j_k(glsim::OLconfiguration &conf)
{
  std::vector<dcomplex> rk(3,dcomplex(0));

  for (int i=0; i<conf.N; i++) {
    double kr= k[0]*conf.r[i][0] + k[1]*conf.r[i][1] + k[2]*conf.r[i][2];
    rk[0]+=conf.v[i][0]*exp(dcomplex(0,-1)*kr);
    rk[1]+=conf.v[i][1]*exp(dcomplex(0,-1)*kr);
    rk[2]+=conf.v[i][2]*exp(dcomplex(0,-1)*kr);
  }
  return rk;
}

Cko& Cko::push_config(glsim::OLconfiguration &conf)
{
  Npart=conf.N;
  std::vector<dcomplex> jk_;
  jk_=j_k(conf);
  jkx_.push_back(jk_[0]);
  jky_.push_back(jk_[1]);
  jkz_.push_back(jk_[2]);

  return *this;
}

//
// Although the correlation function transforms jk, computes squared
// modulus and transforms back, I don't simply skip the back transform
// part to get the correlation in omega space, because the
// normalization of the time correlation is done in the time space
Cko& Cko::compute_Cko()
{
  dcomplex fac=1./dcomplex(Npart,0);
  int Ckolen=jkx_.size();
  Cko_.clear();
  Cko_.resize(Ckolen,dcomplex(0,0));
  cFFT FF(glsim::FFT::in_place);

  correlation_1d_tti_fft(jkx_,FF,Ckolen);
  FF.forward();
  for (int j=0; j<Ckolen; j++) Cko_[j]+=fac*FF.fdatum(j);
  correlation_1d_tti_fft(jky_,FF,Ckolen);
  FF.forward();
  for (int j=0; j<Ckolen; j++) Cko_[j]+=fac*FF.fdatum(j);
  correlation_1d_tti_fft(jkz_,FF,Ckolen);
  FF.forward();
  for (int j=0; j<Ckolen; j++) Cko_[j]+=fac*FF.fdatum(j);

  delta_omega=2*M_PI/(deltat*Ckolen);

  return *this;
}

std::ostream& operator<<(std::ostream& o,const Cko& Cko_)
{
  o << "# omega   Cko'   Cko''  |Ck|\n";
  for (int i=0; i<Cko_.Cko_.size(); i++)
    o << i*Cko_.delta_omega << "  " << Cko_.Cko_[i].real()
      << "  " << Cko_.Cko_[i].imag() << " " << std::norm(Cko_.Cko_[i]) << '\n';
}

/******************************************************************************
 *
 * Options and main
 *
 */

struct optlst {
public:
  int  kn,kdir;
  std::vector<std::string> ifiles;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("cko")
{
  hidden_command_line_options().add_options()
    ("kn",po::value<int>(&options.kn)->required(),"wavevector")
    ("kdir",po::value<int>(&options.kdir)->required(),"k direction")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  // command_line_options().add_options()
  //   ("nave,n",po::value<int>(&options.nave)->default_value(10),
  //    "do arg averages over random directions of the wavevector")
  //   ("seed,S",po::value<long>(&options.seed)->default_value(1L),
  //    "random number seed (with -s used to generate only one random direction)")
  //   ;

  positional_options().add("kn",1).add("kdir",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] kn kdir ifile [ifile ....]\n\n"
    << "Computes C(k,omega) (def2) at wave one vector in one of the x,y,z directions\n"
    << "(specified by kdir=0,1,2) and kn (the desired multiple of Delta k).\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main and deltat
 *
 */

double get_deltat(glsim::H5_multi_file &ifs,glsim::OLconfiguration &conf)
{
  ifs.read();
  double t0=conf.time;
  ifs.read();
  ifs.rewind();
  return conf.time-t0;
}

void substract_cm(glsim::OLconfiguration &conf)
{
  std::vector<double> CM=conf.center_of_mass();
  for (int i=0; i<conf.N; ++i) {
    conf.r[i][0]-=CM[0];
    conf.r[i][1]-=CM[1];
    conf.r[i][2]-=CM[2];
  }
}

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  if (options.kdir<0 || options.kdir>2)
    throw glsim::Runtime_error("kdir must be 0,1,2");

  double deltat=get_deltat(ifs,conf);
  Cko C(conf.box_length,deltat,options.kn,options.kdir);
  while (ifs.read()) {
    conf.unfold_coordinates();
    substract_cm(conf);
    C.push_config(conf);
  }
  C.compute_Cko();
  std::cout << C;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
