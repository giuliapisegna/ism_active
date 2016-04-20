/*
 * ck.cc -- Def2 in k space
 *
 *
 */

#include "parameters.hh"
#include "olconfiguration.hh"
#include "timecorr.hh"

#define HAVE_LIBFFTW2
#ifdef HAVE_LIBFFTW3
typedef glsim::RealFFTW     rFFT;
typedef glsim::ComplexFFTW  cFFT;
#else
typedef glsim::RealFFT_gsl_2n     rFFT;
typedef glsim::ComplexFFT_gsl_2n  cFFT;
#endif /* HAVE_LIBFFTW3 */

class Ck {
public:
  Ck(double k,double deltat,int Nav);

  Ck& push_config(double r[][3],int N);
  Ck& compute_Fk();
  const vcomplex& Ck_data() const {return Ck_;} 

private:
  int                               Nav,Npart;
  double                            k,deltat;
  std::vector<std::vector<double> > kr;
  std::vector<vcomplex>             jkx_,jky_,jkz_;
  vcomplex                          Ck_;
  
  dcomplex rho_k(double r[][3],double k[]);

  friend std::ostream& operator<<(std::ostream&,const Fk&); 
} ;

Ck::Ck(double box_length,double deltat_,int Nav_) :
  k(k_),
  deltat(deltat_),
  Nav(Nav_),
  kr(Nav_),
  rhok_(Nav_),
  Npart(0)
{
  // choose deltak
  deltak_[0]=2*M_PI/box_length[0];
  deltak_[1]=2*M_PI/box_length[1];
  deltak_[2]=2*M_PI/box_length[2];
}

std::vector<dcomplex> Ck::j_k(double k[],double r[][3],double v[][3])
{
  std::vector<dcomplex> rk(3,dcomplex(0));

  for (int i=0; i<Npart; i++) {
    double kr= k[0]*r[i][0] + k[1]*r[i][1] + k[2]*r[i][2];
    rk[0]+=v[i][0]*exp(dcomplex(0,-1)*kr);
    rk[1]+=v[i][1]*exp(dcomplex(0,-1)*kr);
    rk[2]+=v[i][2]*exp(dcomplex(0,-1)*kr);
  }
  return rk;
}

Ck& Ck::push_config(glsim::OLconfiguration &conf)
{
  Npart=N;

  std::vector<dcomplex> jk_;
  for (int i=0; i<Nav; i++) {
    jk__=j_k(&(kr[i][0]),conf.r,conf.v);
    jkx_[i].push_back(jk_[0]);
    jky_[i].push_back(jk_[1]);
    jkz_[i].push_back(jk_[2]);
  }

  return *this;
}

Ck& Ck::compute_Ck()
{
  // compute Ck for each direction and average
  dcomplex fac=1./dcomplex(Npart*Nav,0);
  int Cklen=rhok_[0].size()/2;
  Ck_.clear();
  Ck_.resize(Cklen,dcomplex(0,0));
  cFFT FF(FFT::in_place);
  for (int i=0; i<Nav; i++) {
    correlation_1d_tti_fft(jkx_[i],FF,Cklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Ck_[j]+=fac*FF.tdatum(j);
    correlation_1d_tti_fft(jky_[i],FF,Cklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Ck_[j]+=fac*FF.tdatum(j);
    correlation_1d_tti_fft(jkz_[i],FF,Cklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Ck_[j]+=fac*FF.tdatum(j);
  }
  return *this;
}

std::ostream& operator<<(std::ostream& o,const Fk& Fk_)
{
  o << "# time   Ck'   Ck''\n";
  for (int i=0; i<Fk_.Fk_.size(); i++)
    o << i*Fk_.deltat << "  " << Fk_.Fk_[i].real()
      << "  " << Fk_.Fk_[i].imag() << '\n';
}

/*****************************************************************************
 *
 * F_s(k,t)
 *
 */

Fsk::Fsk(double k_,double deltat_) :
  k(k_),
  deltat(deltat_),
  Npart(0)
{
  // choose a random dir
  Spherical3d_distribution sr;
  kr.resize(3);
  sr(&(kr[0]));
  for (int j=0; j<3; j++) kr[j]*=k;
}

Fsk& Fsk::push_config(double r[][3],int N)
{
  Npart=N;
  expkr.resize(N);
  
  for (int i=0; i<N; i++)
    expkr[i].push_back( exp( dcomplex(0,-1)*
			     (kr[0]*r[i][0] + kr[1]*r[i][1] + 
			      kr[2]*r[i][2]) ) );
  
  return *this;
}

Fsk& Fsk::compute_Fsk()
{
  dcomplex fac=dcomplex(1./Npart,0);
  int Fklen=expkr[0].size()/2;
  Fsk_.clear();
  Fsk_.resize(Fklen);
  cFFT FF(FFT::in_place);
  for (int i=0; i<Npart; i++) {
    correlation_1d_tti_fft(expkr[i],FF,Fklen);
    for (int j=0; j<FF.tdata_rw().size(); j++) Fsk_[j]+=fac*FF.tdatum(j);
  }
  
}

std::ostream& operator<<(std::ostream& o,const Fsk& Fsk_)
{
  o << "# time   F_s'(k)   F_s''(k)\n";
  for (int i=0; i<Fsk_.Fsk_.size(); i++)
    o << i*Fsk_.deltat << "  " << Fsk_.Fsk_[i].real()
      << "  " << Fsk_.Fsk_[i].imag() << '\n';
}

} /* namespace */




struct optlst {
public:
  double k;
  std::vector<std::string> ifiles;
  bool   self_part;
  bool   substract_cm;
  int    nave;
  long   seed;
} options;

class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("GS_olconf_dump")
{
  hidden_command_line_options().add_options()
    ("k",po::value<double>(&options.k)->required(),"wavevector")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("substract-cm,c",po::bool_switch(&options.substract_cm)->default_value(false),
     "substract center of mass at each time (assumes all particles have the same mass)")
    ("self,s",po::bool_switch(&options.self_part)->default_value(false),
     "compute only the self part F_s(k) [default is the full F(k,t)]")
    ("nave,n",po::value<int>(&options.nave)->default_value(10),
     "do arg averages over random directions of the wavevector")
    ("seed,S",po::value<long>(&options.seed)->default_value(1L),
     "random number seed (with -s used to generate only one random direction)")
    ;

  positional_options().add("k",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] k ifile [ifile ....]\n\n"
    << "Computes the (self part of) the intermediate scattering function,"
    << "at the given wavevector and from the given trajectory files.\n"
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
  glsim::Random_number_generator RNG(glsim::gsl_rng_mt19937,options.seed);

  double deltat=get_deltat(ifs,conf);
  glsim::Fk F(options.k,deltat,options.nave);
  while (ifs.read()) {
    // conf.unfold_coordinates();
    // if (options.substract_cm) substract_cm(conf);
    F.push_config(conf.r,conf.N);
  }
  F.compute_Fk();
  std::cout << F;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
