/*
 * clambdat.cc -- Compute C_\lambda(t)
 *
 */

#include "cmarray.hh"

typedef column_major_array<int>  LaplacianT;

void laplacian_matrix(glsim::OLconfiguration &conf,double rc,LaplacianT &L)
{
  glsim::NeighbourList_subcells NN(rc);
  NN.rebuild(conf,rc);

  memcpy(L.data(),0,L.nx()*L.ny());

  for (auto p = NN->pairs_begin(); p!=NN->pairs_end(); ++p) {
    L(p->first,p->second)--;
    L(p->second,p->first)--;
    L(p->first,p->first)++;
    L(p->second,p->second)++;
  }
}

void compute_avecs(LaplacianT &L,double *avals,AvecT &avecs)
{
  ifail=new lapack_int[L.nx()];
  // lapack_int LAPACKE_dsyevx( int matrix_order, char jobz, char range, char uplo,
  //                            lapack_int n, double* a, lapack_int lda, double vl,
  //                            double vu, lapack_int il, lapack_int iu,
  //                            double abstol, lapack_int* m, double* w, double* z,
  //                            lapack_int ldz, lapack_int* ifail );

  // double avals[N];
  // column_major_array<double> avecs(conf.N,desired_avecs);

  info=LAPACKE_dsyevx(LAPACK_COL_MAJOR,'V','I','U'
		      L.nx(),L.data(),L.nx(),0,
		      0,1,10,
		      0,&M,avals,avecs,
		      conf.N,ifail);
}

/*****************************************************************************
 *
 * options
 *
 */

struct optlst {
public:
  double rc;
  bool normalize,connect;
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
    ("connect,C",po::bool_switch(&options.connect)->default_value(false),
     "Computed connected C(k,t) (Roman style)")
    ("normalize,N",po::bool_switch(&options.normalize)->default_value(false),
     "Normalize correlation to 1 at t=0")
     ;

  positional_options().add("ifiles",-1);
}

void CLoptions::show_usage() const
{
  std::cerr
    << "usage: " << progname << "[options] ifile [ifile ....]\n\n"
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

void wmain(int argc,char *argv[])
{
  CLoptions o;
  o.parse_command_line(argc,argv);
  
  glsim::OLconfiguration conf;
  glsim::OLconfig_file   cfile(&conf);
  glsim::H5_multi_file   ifs(options.ifiles,cfile);

  double deltat=get_deltat(ifs,conf);
  Ckt C(conf.box_length,deltat,options.kn,options.kdir,options.connect);

  column_major_array L(conf.N,conf.N);

  ifs.read();
  laplacian_matrix(conf,options.rc,L);
  do {
    conf.unfold_coordinates();
    substract_cm(conf);
    // C.push_config(conf);
  }  while (ifs.read());
  C.compute_Ckt();
  std::cout << C;
}

int main(int argc, char *argv[])
{
  return glsim::StandardEC(argc,argv,wmain);
}
