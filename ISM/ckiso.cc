/*
 * ckiso.cc -- Compute the isotropic (averaged over k) C(k)
 *
 */

#include <valarray>
#include <gsl/gsl_sf_trig.h>

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/timecorr.hh"
#include "glsim/binvec.hh"

#include "3dvecs.hh"


/******************************************************************************
 *
 * This computes C(r)
 *
 */
class Cr {
public:
  Cr(const glsim::OLconfiguration& c,double rnn);
  ~Cr();

  /// Number of bins
  int    size() const {return vec->nbins();}
  /// \f$\Delta r\f$
  double deltar() const {return vec->delta();}
  /// Radius corresponding to ith bin
  double r(int i) const {return vec->binc(i);}
  /// 4\pi\rho r^2 C(r) at ith bin, excluding the self (i=j) contribution
  double corr_nonorm(std::size_t i) const
  {return (*vec)[i]/nconf;}

  void push(const glsim::OLconfiguration&);
  void push(const Cr&);

private:
  double Cr0;
  glsim::Binned_vector<double> *vec;
  double rnn;
  int    nconf;
  double rho;

  friend std::ostream& operator<<(std::ostream& o,const Cr& corr);
} ;

Cr::Cr(const glsim::OLconfiguration &c,double rnn_) :
  rnn(rnn_), nconf(0), Cr0(0.)
{
  double rmax=sqrt(c.box_length[0]*c.box_length[0] + c.box_length[1]*c.box_length[1] +
		   c.box_length[2]*c.box_length[2])/1.99;
  int nbins=(int) ceil(rmax/(0.1*rnn));
  vec=new glsim::Binned_vector<double>(nbins,0,rmax);
  for (int i=0; i<nbins; (*vec)[i++]=0.) ;
  rho=c.number_density();
}

Cr::~Cr()
{
  delete vec;
}

void Cr::push(const glsim::OLconfiguration& conf)
{
  for (int i=0; i<conf.N; ++i) {
    Cr0+=dotp(conf.v[i],conf.v[i])/conf.N;
  }
  for (int i=0; i<conf.N-1; ++i) {
    for (int j=i+1; j<conf.N; ++j) {
      double r=sqrt(conf.distancesq(i,j));
      (*vec)[r]+=2*dotp(conf.v[i],conf.v[j])/conf.N;
    }
  }
  nconf++;
}

void Cr::push(const Cr& C)
{
  Cr0+=C.Cr0;
  for (int i=0; i<vec->nbins(); ++i)
    (*vec)[i]+=(*(C.vec))[i];
  nconf+=C.nconf;
}
 
std::ostream& operator<<(std::ostream& o,const Cr& corr)
{
  double f=1./(4*M_PI*corr.rho*corr.nconf*corr.vec->delta());
  o << 0. << " " << corr.Cr0/corr.nconf << '\n';
  for (int i=0; i < corr.vec->nbins(); ++i) {
    double r=corr.vec->binc(i);
    o << r << "  " << f*(*corr.vec)[i]/(r*r) << '\n';
  }
}

/*****************************************************************************
 *
 * C(k)
 *
 */

class Ckiso {
public:
  Ckiso(const glsim::OLconfiguration&,const Cr&,int);
  std::vector<double> Ck;
  double              deltak;
} ;

Ckiso::Ckiso(const glsim::OLconfiguration &conf,const Cr& C,int nk)
{
  deltak=0.2*M_PI/conf.box_length[0];

  Ck.resize(nk);
  for (int ik=0; ik<nk; ++ik) {
    double k=ik*deltak/M_PI; // Divide by PI because of GSL's definition 
                             // of the sinc function
  // We start with 1 instead of 0 because corr_nonorm, which is used
  // in the loop below, does note include the self (i=j) contribution
    Ck[ik]=1;
    for (int i=0; i<C.size(); ++i) {
      double kr=k*C.r(i);
      Ck[ik]+=C.corr_nonorm(i)*gsl_sf_sinc(kr);
    }
  }
}

std::ostream& operator<<(std::ostream& o,const Ckiso& ck)
{
  for (int i=0; i < ck.Ck.size(); ++i) {
    double k=i*ck.deltak;
    o << k << "  " << ck.Ck[i] << '\n';
  }
}


/******************************************************************************
 *
 * Options and main
 *
 */

struct optlst {
public:
  int  nk;
  double rnn;
  std::vector<std::string> ifiles;
  std::string c_of_r_file;

} options;


class CLoptions : public glsim::UtilityCL {
public:
  CLoptions();
  void show_usage() const;
} ;

CLoptions::CLoptions() : glsim::UtilityCL("ckiso")
{
  hidden_command_line_options().add_options()
    ("nk",po::value<int>(&options.nk)->required(),"number of wavevectors to compute")
    ("rnn",po::value<double>(&options.rnn)->required(),"mean NN distance")
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("c-of-r,r",po::value<std::string>(&options.c_of_r_file),
     "also write C(r) to given file")
     ;

  positional_options().add("nk",1).add("rnn",1).add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] nk rnn ifile [ifile ....]\n\n"
    << "nk is the number of k values, rnn is the mean nearest-neighbour distance.\n"
    << "This computes the isotropic C(k) for nk wavevectors kn (Delta k computed automatically).\n"
    << "This computes\n\n"
    << "      C(k) = (1/N) \\sum_{ij} \\delta \\hat v_i \\delta \\hat v_i \\sinc( k r_{ij})\n\n"
    << "where \\delta v_i = v_i/v_0 - (1/N) \\sum_i v_i/v_0 and \\delta\\hat v_i is\n"
    << "normalized.\n"
    << "\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}


/*****************************************************************************
 *
 * main and normalize_vel
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

  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  ifs.read();
  Cr C(conf,options.rnn);
  ifs.rewind();

  #pragma omp parallel
  {
    glsim::OLconfiguration confloc=conf;
    Cr Cloc(confloc,options.rnn);

    #pragma omp for schedule(static) nowait
    for (hsize_t rec=0; rec<ifs.size(); ++rec) {
      #pragma omp critical
      {
	ifs.read();
	confloc=conf;
      }
      normalize_vel(confloc);
      Cloc.push(confloc);
    }

    #pragma omp critical
    C.push(Cloc);
  }

  if (options.c_of_r_file.length()>0) {
    std::ofstream rfile(options.c_of_r_file);
    int ixi=0;
    while (C.corr_nonorm(ixi)>0 && ixi <C.size()) ++ixi;
    rfile << "# Velocity correlation (real space)\n";
    rfile << "#\n# r0 (crossing zero) = " << C.r(ixi) << '\n';
    rfile << "#\n# r    C(r)\n";
    rfile << C;
  }
  Ckiso Ck(conf,C,options.nk);
  std::cout << "# k   C(k)\n";
  std::cout << Ck;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
