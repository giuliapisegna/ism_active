/*
 * social.hh -- Header for social interactions
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */

#ifndef SOCIAL_HH
#define SOCIAL_HH

#include "glsim/parameters.hh"
#include "glsim/olconfiguration.hh"
#include "glsim/nneighbours.hh"

/*****************************************************************************
 * 
 * SocialInteracions
 * 
 */

class SocialInteractions {
public:
  SocialInteractions(glsim::OLconfiguration&) {}
  virtual ~SocialInteractions() {}

  virtual double social_mass(short type) const=0;
  virtual double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3])=0;
  virtual double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3])=0;
  virtual void fold_coordinates(glsim::OLconfiguration&,double maxdisp=-1);
} ;

inline void SocialInteractions::fold_coordinates(glsim::OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
}

/*****************************************************************************
 *
 * Vicsek
 * 
 */

class VicsekParameters : public glsim::Parameters {
public:
  VicsekParameters(const char *scope=glsim::Parameters::default_scope);
} ;

class VicsekInteraction : public SocialInteractions {
public:
  VicsekInteraction(glsim::OLconfiguration&,glsim::NearestNeighbours *NN=0);
  double social_mass(short type) const {return mass;}
  double social_potential_energy_and_force(glsim::OLconfiguration&,double b[][3]) {}
  double social_potential_energy_and_acceleration(glsim::OLconfiguration&,double b[][3]);
  void   fold_coordinates(glsim::OLconfiguration&,double maxdisp=-1);

private:
  VicsekParameters  par;
  double            mass,Jsv0sq,rc,rcsq;
  bool              metric;

  bool                     own_NN;
  glsim::NearestNeighbours *NN;
} ;

inline void VicsekInteraction::fold_coordinates(glsim::OLconfiguration& conf,double maxdisp)
{
  conf.fold_coordinates();
  if (maxdisp<0)
    NN->rebuild(conf,rc);
  else
    NN->update(conf,maxdisp);
}

#endif /* SOCIAL_HH */
