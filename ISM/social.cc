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

/*****************************************************************************
 * 
 * Vicsek
 *
 */

/*
 * Parameters
 *
 */

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

/*
 * Generic Vicsek
 *
 */

VicsekInteraction::VicsekInteraction(VicsekParameters &p,
				     glsim::OLconfiguration &c) :
  SocialInteractions(),
  par(p)
{
  v0=par.value("Vicsek.v0").as<double>();
  v0sq=v0*v0;
  chi=par.value("Vicsek.chi").as<double>();
  J=par.value("Vicsek.J").as<double>();
  metric=par.value("Vicsek.metric").as<bool>();
  rc=par.value("Vicsek.cutoff").as<double>();
}

/*
 * Topological Vicsek
 *
 */

TopologicalVicsekInteraction::TopologicalVicsekInteraction(VicsekParameters &p,
						glsim::OLconfiguration &c,
						glsim::TopologicalNearestNeighbours *n) :
  VicsekInteraction(p,c)
{
  if (metric) throw glsim::Runtime_error("You asked for metric interactions but created the topological object");
  NNeighbours=(int) rc;

  if (n) {
    NN=n;
    own_NN=false;
  } else {
    NN=new glsim::TopologicalNeighbours_naive(NNeighbours);
    own_NN=true;
  }
  NN->rebuild(c,rc);
}

/*
   This is computes the social interactions, but in the topological case.

   This routine assumes that the input velocities have modulus v0 (the
   modulus of the velocity is fixed in the Vicsek model).

   Note that the acceleration can be computed as F/m or as v0sq*F/chi.

*/
double TopologicalVicsekInteraction::social_potential_energy_and_acceleration(glsim::OLconfiguration &conf,
									      double b[][3])
{
  double E=0;
  memset(b,0,conf.N*3*sizeof(double));

  double efac=-J/v0sq;
  double ffac=J/chi;

  for (int i=0; i<conf.N; i++) {
    for (auto p = NN->neighbours_begin(i); p!=NN->neighbours_end(i); ++p) {

      E += efac*dotp(conf.v[i],conf.v[*p]);

      b[i][0] += ffac * conf.v[*p][0];
      b[i][1] += ffac * conf.v[*p][1];
      b[i][2] += ffac * conf.v[*p][2];

    }
  }
  return E;
}
