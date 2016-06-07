/*
 * clambdat.cc -- Compute C_\lambda(t)
 *
 */

#define HAVE_LIBFFTW3

#include "glsim/blib.hh"
#include "glsim/offlattice.hh"
#include "glsim/timecorr.hh"

#include <lapacke.h>
#ifdef HAVE_LIBFFTW3
typedef glsim::RealFFTW     rFFT;
typedef glsim::ComplexFFTW  cFFT;
#else
typedef glsim::RealFFT_gsl_2n     rFFT;
typedef glsim::ComplexFFT_gsl_2n  cFFT;
#endif /* HAVE_LIBFFTW3 */

#include "3dvecs.hh"
#include "cmarray.hh"

typedef column_major_array<double>  LaplacianT,EvecT;

void laplacian_matrix(glsim::OLconfiguration &conf,double rc,LaplacianT &L)
{
  glsim::NeighbourList_subcells NN(rc);
  NN.rebuild(conf,rc);

  memset(L.data(),0,L.nrows()*L.ncols());

  for (auto p = NN.pairs_begin(); p!=NN.pairs_end(); ++p) {
    L(p->first,p->second)--;
    L(p->second,p->first)--;
    L(p->first,p->first)++;
    L(p->second,p->second)++;
  }
  // Test laplacian: linear periodic chain
  // L(0,0)=2;
  // L(0,1)=-1;
  // L(0,L.ncols()-1)=-1;
  // for (int i=1; i<L.nrows()-1; ++i) {
  //   L(i,i)=2;
  //   L(i,i+1)=L(i,i-1)=-1;
  // }
  // L(L.nrows()-1,0)=-1;
  // L(L.nrows()-1,L.nrows()-2)=-1;
  // L(L.nrows()-1,L.nrows()-1)=2;
}

void check_laplacian_tralation_invariance(LaplacianT &L,int dim)
{
  double sum;
  for (int row=0; row<L.nrows(); ++row) {
    sum=0;
    for (int col=0; col<L.ncols(); col+=dim)
      sum+=L(row,col);
    std::cout << "Sum of row " << row << " stride " << dim <<
      " = " << sum;
  }
}

void compute_evecs(LaplacianT &L,double *evals,EvecT &evecs)
{
  int* ifail=new lapack_int[L.ncols()];
  // lapack_int LAPACKE_dsyevx( int matrix_order, char jobz, char range, char uplo,
  //                            lapack_int n, double* a, lapack_int lda, double vl,
  //                            double vu, lapack_int il, lapack_int iu,
  //                            double abstol, lapack_int* m, double* w, double* z,
  //                            lapack_int ldz, lapack_int* ifail );

  int M;
  // To request only evals
  // int info=LAPACKE_dsyevx(LAPACK_COL_MAJOR,'N','I','U',
  // 			  L.nrows(),L.data(),L.nrows(),0,
  // 			  0,1,evecs.ncols(),
  // 			  0,&M,evals,evecs.data(),
  // 			  L.nrows(),ifail);
  // To compute evals and evectors
  int info=LAPACKE_dsyevx(LAPACK_COL_MAJOR,'V','I','U',
  			  L.nrows(),L.data(),L.nrows(),0,
  			  0,1,evecs.ncols(),
  			  0,&M,evals,evecs.data(),
  			  L.nrows(),ifail);
  assert(info==0);
  delete[] ifail;
}

/*****************************************************************************
 *
 * options
 *
 */

struct optlst {
public:
  double rc;
  bool normalize;
  int  Neval;
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
    ("ifiles",po::value<std::vector<std::string> >(&options.ifiles)->required(),"input files")
    ;
  command_line_options().add_options()
    ("cut-off,c,d",po::value<double>(&options.rc)->required(),
     "cut-off radius (metric)")
    ("neigen,n",po::value<int>(&options.Neval)->default_value(10),
     "compute correlation for the first arg eigenvalues/eigenvectors")
    ("normalize,N",po::bool_switch(&options.normalize)->default_value(false),
     "Normalize correlation to 1 at t=0")
     ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] ifile [ifile ....]\n"
    << "MUST give -c\n"
    << " Options:\n";
  show_command_line_options(std::cerr);
}

/*****************************************************************************/

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

void substract_polarization(glsim::OLconfiguration& conf)
{
  double V[3]={0,0,0};
  for (int i=0; i<conf.N; ++i) {
    V[0]+=conf.v[i][0];
    V[1]+=conf.v[i][1];
    V[2]+=conf.v[i][2];
  }
  V[0]/=conf.N;
  V[1]/=conf.N;
  V[2]/=conf.N;
  for (int i=0; i<conf.N; ++i) {
    double v0=sqrt(modsq(conf.v[i]));
    conf.v[i][0]=(conf.v[i][0]-V[0])/v0;
    conf.v[i][1]=(conf.v[i][1]-V[1])/v0;
    conf.v[i][2]=(conf.v[i][2]-V[2])/v0;
  }
}

void wmain(int argc,char *argv[])
{
  glsim::logs.set_stream(std::cerr,glsim::warn);
  
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  double deltat=get_deltat(ifs,conf);
  // Ckt C(conf.box_length,deltat,options.kn,options.kdir,options.connect);

  LaplacianT L(conf.N,conf.N);
  double     *evals=new double[options.Neval];
  EvecT      evecs(conf.N,options.Neval);

  ifs.read();
  laplacian_matrix(conf,options.rc,L);
  compute_evecs(L,evals,evecs);

  // for (int n=0; n<conf.N; ++n) {
  //   std::cout << "eval " << n << ": " << evals[n] << "\n";
  // }
  // std::cout << "\n\n\n";
  // for (int n=0; n<3; ++n) {
  //   std::cout << "eval " << n << ": " <<  evals[n] << "\n\n";
  //   std::cout << "evec " << n << '\n';
  //   for (int i=0; i<conf.N; ++i)
  //     std::cout << "component " << i << ": " << evecs(i,n) << '\n';
  //   std::cout << "\n\n\n";
  // }
  // exit(1);


  double V[3];
  typedef std::vector<double> seriesT;
  std::vector<double> deltaVx(conf.N),deltaVy(conf.N),deltaVz(conf.N);
  std::vector<seriesT> dvx(options.Neval),dvy(options.Neval),dvz(options.Neval),
    Corr(options.Neval);

  do {
    conf.unfold_coordinates();
    substract_cm(conf);

    // Compute \deltav_i (site base)
    substract_polarization(conf);
    // Project \delta v onto eigenvector and store
    for (int n=0; n<options.Neval; ++n) {
      double dx,dy,dz;
      dx=dy=dz=0;
      for (int i=0; i<conf.N; ++i) {
	double v0=sqrt(modsq(conf.v[i]));
	dx+=evecs(i,n)*conf.v[i][0];
	dy+=evecs(i,n)*conf.v[i][1];
	dz+=evecs(i,n)*conf.v[i][2];
      }
      dvx[n].push_back(dx);
      dvy[n].push_back(dy);
      dvz[n].push_back(dz);
    }

  }  while (ifs.read());

  // Compute autocorrelations
  rFFT FF(glsim::FFT::in_place);
  for (int n=0; n<options.Neval; ++n) {
    Corr[n].resize(dvx[n].size()/2);
    correlation_1d_tti_fft(dvx[n],FF,dvx[n].size()/2);
    for (int j=0; j<FF.tdata_rw().size(); j++)
      Corr[n][j]+=FF.tdatum(j);
    correlation_1d_tti_fft(dvy[n],FF,dvx[n].size()/2);
    for (int j=0; j<FF.tdata_rw().size(); j++)
      Corr[n][j]+=FF.tdatum(j);
    correlation_1d_tti_fft(dvz[n],FF,dvx[n].size()/2);
    for (int j=0; j<FF.tdata_rw().size(); j++)
      Corr[n][j]+=FF.tdatum(j);
  }

  // Print
  std::cout << "# Eigenvalues\n\n";
  for (int n=0; n<options.Neval; ++n)
    std::cout << "# lambda" << n << "= " << evals[n] << '\n';
  std::cout << "#   C_lambda0  C_lambda1 ...\n";
  std::vector<double> fac(options.Neval,1.);
  if (options.normalize)
    for (int n=0; n<options.Neval; ++n)
      fac[n]=1./Corr[n][0];
  for (int i=0; i<Corr[0].size(); i++) {
    std::cout << i*deltat;
    for (int n=0; n<options.Neval; ++n)
      std::cout << "  " << fac[n]*Corr[n][i];
    std::cout << '\n';
  }

  delete[] evals;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
