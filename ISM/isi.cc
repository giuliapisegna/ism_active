/*
 * isi.cc -- The inertial spin integrator
 *
 * Integrates the equation of motion for an inertial spin model (with
 * Brownian noise), without specifiying the velocity-velocity
 * interaction.
 *
 * Tomás S. Grigera <tgrigera@iflysib.unlp.edu.ar>
 *
 * La Plata, July 2015
 *
 */


ISMSimulation::ISMSimulation(ISMEnvironment& e,OLConfiguration &c,glsim::Interactions *i) :
  Simulation(e,c),
  env(e),
  conf(c),
  inter(i)
{
  if (conf.a==0) {
    conf.a=new double[conf.N][3];
    memset(conf.a,0,conf.N*3*sizeof(double));
  }
  if (conf.v==0) {
    conf.v=new double[conf.N][3];
    memset(conf.v,0,conf.N*sizeof(double));
  }

  env.total_number=conf.N;
  env.total_mass=0;
  env.inter=inter;
  for (int i=0; i<conf.N; ++i) {
    env.total_mass+=inter->mass(conf.type[i]);
    env.Ptot[0]+=conf.v[i][0]*inter->mass(conf.type[i]);
    env.Ptot[1]+=conf.v[i][1]*inter->mass(conf.type[i]);
    env.Ptot[2]+=conf.v[i][2]*inter->mass(conf.type[i]);
  }
  // Substract Vcm (if P is conserved)
  if (inter->conserve_P())
    for (int i=0; i<conf.N; ++i) {
      conf.v[i][0]-=env.Ptot[0]/env.total_mass;
      conf.v[i][1]-=env.Ptot[1]/env.total_mass;
      conf.v[i][2]-=env.Ptot[2]/env.total_mass;
    }

  env.Epot=inter->potential_energy(conf);
  conf.step=env.steps_completed;
  conf.time=env.time_completed;
  update_observables();

  // Constants for Langevin integration
  Dt=env.time_step;
  double chi=inter->chi(conf.type[i]);
  double xi=env.eta/chi;
  double xidt=xi*Dt;
  double c0l,c1,c2,sx,sv,rho;
  double exi=exp(-xidt);
  if (xidt<1e-3) {
    c0l=1 - xidt + xidt*xidt/2 - xidt*xidt*xidt/6;
    c1=1 - xidt/2 + xidt*xidt/6 - xidt*xidt*xidt/24;
    c2=0.5 - xidt/6 + xidt*xidt/24;
    rho=sqrt(3.)*(0.5-xidt/16.-(17./1280.)*xidt*xidt
		  +(17./6144)*xidt*xidt*xidt);
  } else {
    c0l=exi;
    c1=(1-c0l)/xidt;
    c2=(1-c1)/xidt;
    rho=(1-exi)*(1-exi)/sqrt( (1-exi*exi)*(3*xidt-3+4*exi-exi*exi) );
  }
  sa=(env.temperature/chi)*(1-exi*exi);
  sa=sqrt(sa);
  if (env.eta<1e-3) {
    sv=env.temperature*Dt*Dt*Dt*env.eta*(2./3.-0.5*xidt)/(chi*chi);
  } else {
    sv=(env.temperature/env.eta)*(2*Dt-(3-4*exi+exi*exi)/xi);
  }
  sv=sqrt(sv);
  noise=new BivariateGaussian_distribution(sv,sa,rho);
  
  c0=c0l;
  c1dt=Dt*c1;
  c2dt=Dt*c2;
  c1mc2=c1-c2;

  v0sq=0;
}

ISMSimulation::~ISMSimulation()
{
  delete noise;
}

void ISMSimulation::step()
{
  double w,xiv[3],xia[3],deltav[3];
  
  for (int i=0; i<conf.N; ++i) {

    // Get noise*3
    noise(xiv[0],xia[0]);
    noise(xiv[1],xia[1]);
    noise(xiv[2],xia[2]);

    // 1. Update x
    conf.r[i][0] += conf.v[i][0]*Dt;
    conf.r[i][1] += conf.v[i][1]*Dt;
    conf.r[i][2] += conf.v[i][2]*Dt;

    // 2. Delta v w/o rattle corrections
    deltav[0] = c1dt*conf.a[i][0] + c2dt*Dt*b[i][0] + xiv[0];
    deltav[1] = c1dt*conf.a[i][1] + c2dt*Dt*b[i][1] + xiv[0];
    deltav[2] = c1dt*conf.a[i][2] + c2dt*Dt*b[i][2] + xiv[0];

    // 3. Solve for lagrange multiplier (first part of rattle)
    double a = v0sq;
    double b = 2*dotp(conf.v[i],deltav);
    double c = modsq(deltav)-v0sq;
    if (b==0)
      w = sqrt(-c/a);
    else {
      double q = -0.5 * (b+sgn(b)*sqrt(b*b-4*a*c));
      w = b>0 ? c/q : q/a;
    }
    double hlambda = (w -1)/c2Dt;

    // 4. Partial update of a (missing second rattle correction and
    // acceleration at next step)
    conf.a[i][0] = c0*conf.a[i][0] + c1mc2*(Dt*b[i][0] + hlambda*conf.v[i][0]) + xia[0];
    conf.a[i][1] = c0*conf.a[i][1] + c1mc2*(Dt*b[i][1] + hlambda*conf.v[i][1]) + xia[1];
    conf.a[i][2] = c0*conf.a[i][2] + c1mc2*(Dt*b[i][2] + hlambda*conf.v[i][2]) + xia[2];

    // 5. Update v
    conf.v[i][0] = (1 + c2dt*hlambda)*conf.v[i][0] + deltav[0];
    conf.v[i][1] = (1 + c2dt*hlambda)*conf.v[i][1] + deltav[1];
    conf.v[i][2] = (1 + c2dt*hlambda)*conf.v[i][2] + deltav[2];

  }

  // 6. Compute forces with new positions and velocities
  env.Epot=inter->acceleration_and_potential_energy(conf);

  for (int i=0; i<conf.N; ++i) {

    // 7. Further update of a
    conf.a[i][0] += c2dt*b[i][0];
    conf.a[i][1] += c2dt*b[i][1];
    conf.a[i][2] += c2dt*b[i][2];

    // 8. Solve for second lagrange multiplier (last step of rattle)
    double c2hgamma = - dotp(conf.v[i],conf.a[i])/v0sq;

    // 9. Final correction for a
    conf.a[i][0] += c2hgamma * conf.v[i][0];
    conf.a[i][1] += c2hgamma * conf.v[i][1];
    conf.a[i][2] += c2hgamma * conf.v[i][2];

  }
  inter->fold_coordinates(conf);

  env.time_completed+=Dt;
  env.time_in_run+=Dt;
  env.run_completed = env.steps_in_run>=env.MDsteps;
  conf.time=env.time_completed;
  conf.step=env.steps_completed;
}
