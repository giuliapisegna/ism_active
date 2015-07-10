

#ifndef T3DVECS_HH
#define T3DVECS_HH

inline double dotp(double x[3],double y[3])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline double modsq(double x[3])
{
  return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

#endif /* T3DVECS_HH */

