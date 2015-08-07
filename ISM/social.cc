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
    ("Vicsek.v0",po::value<double>()->required(),"Modulus of velocity")
    ("Vicsek.chi",po::value<double>()->required(),"Social moment of inertia (m*v0sq), actually in the inertial model rather than Vicsek")
    ("Vicsek.J",po::value<double>()->required(),"Vicseck coupling constant")
    ("Vicsek.metric",po::bool_switch()->required(),"True for metric, false for topological interactions")
    ("Vicsek.cutoff",po::value<double>()->required(),"Cutoff for metric interactions")
    ;
}

VicsekInteraction::VicsekInteraction(VicsekParameters &p,glsim::OLconfiguration &c,
				     glsim::NearestNeighbours *n) :
  SocialInteractions(),
  par(p)
{
  v0=par.value("Vicsek.v0").as<double>();
  v0sq=v0*v0;
  chi=par.value("Vicsek.chi").as<double>();
  J=par.value("Vicsek.J").as<double>();
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
   This is computes the social interactions, but only for the metric case.

   This routine assumes that the input velocities have modulus v0 (the
   modulus of the velocity is fixed in the Vicsek model).

   Note that the acceleration can be computed as F/m or as v0sq*F/chi.

*/
double VicsekInteraction::social_potential_energy_and_acceleration(glsim::OLconfiguration &conf,
								   double b[][3])
{
  double E=0;
  memset(b,0,conf.N*3*sizeof(double));

  double efac=-J/v0sq;
  double ffac=J/chi;

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
