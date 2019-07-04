/*
 * social.cc -- Social interactions for inertial spin and related integrators
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#include "social_potential.hh"

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
 // I_AM_HERE;
  parameter_file_options().add_options()
    ("Vicsek.v0",po::value<double>()->required(),"Modulus of velocity")
    ("Vicsek.chi",po::value<double>()->required(),"Social moment of inertia (m*v0sq), actually in the inertial model rather than Vicsek")
    ("Vicsek.J",po::value<double>()->required(),"Vicseck coupling constant")
    ("Vicsek.metric",po::bool_switch()->required(),"True for metric, false for topological interactions")
    ("Vicsek.cutoff",po::value<double>()->required(),"Cutoff for metric interactions")
    ("Vicsek.k", po::value<double>() -> default_value(0.), "Stiffness of the harmonic potential")
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
  k=par.value("Vicsek.k").as<double>();
}


