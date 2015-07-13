/*
 * social.cc -- Social interactions for inertial spin and related integrators
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#include "social.hh"

#include "3dvecs.hh"

VicsekParameters::VicsekParameters(const char* scope) :
  Parameters(scope)
{
  parameter_file_options().add_options()
    ("Vicsek.mass",po::value<double>()->required())
    ("Vicsek.Joverv0sq",po::value<double>()->required(),"Vicseck coupling constant J/v0^2")
    ("Vicsek.metric",po::bool_switch()->required(),"True for metric, false for topological interactions")
    ("Vicsek.cutoff",po::value<double>()->required(),"Cutoff for metric interactions")
    ;
}

VicsekInteraction::VicsekInteraction(VicsekParameters &p,glsim::OLconfiguration &c,
				     glsim::NearestNeighbours *n) :
  SocialInteractions(),
  par(p)
{
  mass=par.value("Vicsek.mass").as<double>();
  Jsv0sq=par.value("Vicsek.Joverv0sq").as<double>();
  metric=par.value("Vicsek.metric").as<bool>();
  if (!metric) throw glsim::Unimplemented("Topological interactions");
  rc=par.value("Vicsek.cutoff").as<double>();
  rcsq=rc*rc;

  if (n) {
    NN=n;
    own_NN=false;
  } else {
    NN=new glsim::NeighbourList_naive(rc);
    own_NN=true;
  }
  NN->rebuild(c,rc);
}

/*
   This is the social force, as is only good for metric interactions
*/
double VicsekInteraction::social_potential_energy_and_acceleration(glsim::OLconfiguration &conf,
								   double b[][3])
{
  double E=0;
  memset(b,0,conf.N*3*sizeof(double));

  double efac=-Jsv0sq;
  double ffac=Jsv0sq/mass;

  for (auto p = NN->pair_begin(); p!=NN->pair_end(); ++p) {

    double rxmn=conf.ddiff(conf.r[p->first][0],conf.r[p->second][0],conf.box_length[0]);
    double rymn=conf.ddiff(conf.r[p->first][1],conf.r[p->second][1],conf.box_length[1]);
    double rzmn=conf.ddiff(conf.r[p->first][2],conf.r[p->second][2],conf.box_length[2]);
    double dsq=rxmn*rxmn+rymn*rymn+rzmn*rzmn;

    if (dsq>rcsq) continue;

    E += efac*dotp(conf.v[p->first],conf.v[p->second]);

    b[p->first][0] += ffac * conf.v[p->second][0];
    b[p->first][1] += ffac * conf.v[p->second][1];
    b[p->first][2] += ffac * conf.v[p->second][2];
    b[p->second][0] += ffac * conf.v[p->first][0];
    b[p->second][1] += ffac * conf.v[p->first][1];
    b[p->second][2] += ffac * conf.v[p->first][2];

  }
  return E;
}
