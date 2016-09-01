#include <pebble.h>
#include "SmallMaths.h"

const float PI     = 3.14159265;
const float TWOPI  = 6.28318531;
const float HALFPI = 1.57079633;
const float LN2    = 0.693147181;
const float RAD2DEG= 57.2957795;
const float DEG2RAD= 0.0174532925;

// --------------------------------------------------------------------------------------
#define SQRT_MAGIC_F 0x5f3759df 
float sm_sqrt(const float x)
{
  const float xhalf = 0.5f*x;
 
  union // get bits for floating value
  {
    float x;
    int i;
  } u;
  u.x = x;
  u.i = SQRT_MAGIC_F - (u.i >> 1);  // gives initial guess y0
  
  u.x = u.x*(1.5f - xhalf*u.x*u.x);   // This can be removed for increasing speed
  u.x = u.x*(1.5f - xhalf*u.x*u.x);   // This can be removed for increasing speed
  return x*u.x*(1.5f - xhalf*u.x*u.x);// Newton step, repeating increases accuracy 
}

// --------------------------------------------------------------------------------------
float sm_exp(float x) {
  int LIMIT = 100;
  int INVERSE_LIMIT = 100000; // Inverse fractional limit 100000 = 0.00001 == 0.001%
  double x1 = 1;// running x, x^2, x^3 etc
  double y = 1; // result
  double yl = 1; // last result
  double yd = 1; // delta result
  double d = 1; // denominator
  for (int i=1; i<=LIMIT; i++) {
    yl = y;
    x1 *= x;
    d *= i;
    y += x1 / d;
    if (y == yl) {
//      printf("sm_exp iterations: %i", i);
      break;
    }
    yd = y / (y - yl);
    if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
//      printf("sm_exp iterations: %i", i);
      break;
    }
  }
  return y;
}

// --------------------------------------------------------------------------------------
float sm_agm(float x1, float x2) {
// Arithmetic-Geometric Mean  
  int LIMIT = 100;
  int MIN_LIMIT = 6;
  int INVERSE_LIMIT = 100000; // Inverse fractional limit 100000 = 0.00001 == 0.001%
  float a1 = 1;
  float g1 = 1000;
  float a2 = x1;
  float g2 = x2;
  
  for (int n=1; n<=LIMIT; n++) {
    a1 = a2;
    g1 = g2;
    a2 = (a1 + g1)/2.0;
    g2 = sm_sqrt(a1 * g1);
    float yd = 1.0f;
   
    if (n >= MIN_LIMIT) {
      if (a2 == g2) {
//        printf("sm_agm iterations: %i", n);
        break;
      }
      yd = g2 / (g2 - g1);
      if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
//        printf("sm_agm iterations: %i", n);
        break;
      }
    }
  }
  
  return a2;
}

// --------------------------------------------------------------------------------------
float sm_powint(float x, int y) {
  float powint = 1.0f;
  if (y > 0) {
    for (int i=0; i<y; i++) {
      powint *= x;
    }
  } else {
    for (int i=0; i<y; i++) {
      powint /= x;    
    }
  }
  return powint;
}

// --------------------------------------------------------------------------------------
float sm_pow(float x, float y) {
  // relies on exp and ln, so you're only going to get 6 sig figs out of this.
  //printf("sm_pow: sm_ln(x)=%i/1000", (int)(1000*sm_ln(x)));
  //printf("sm_pow: y * sm_ln(x)=%i/1000", (int)(1000*(y * sm_ln(x))));
  if (x <= 0.0) {
    return 0.0;
  }
  if (y == 0.0) {
    return 1.0;
  }
  
  if (y == (int)y) {
    return sm_powint(x, (int)y);
  }
  
  float res = sm_exp(y * sm_ln(x));
  return res;
}

// --------------------------------------------------------------------------------------
float sm_ln(float x) {
  const int p = 18; // bits of precision
  int m = 9;
  int twoPowerM = sm_powint(2, m);  
  float s = x * twoPowerM;
  
  while (s < sm_powint(2, p/2)) {
    m+=1;
    twoPowerM = sm_powint(2, m);
    s = x * twoPowerM;
  }
  
  float part1 = 4.0f / s;
  float part2 = sm_agm(1.0f, part1) * 2;
  float part3 = PI / part2;
  float part4 = LN2 * m;
  float y = part3 - part4;
  return y;
}

// --------------------------------------------------------------------------------------
float sm_sind(float angleDegrees) {
  return sm_sin(angleDegrees * DEG2RAD);
}

// --------------------------------------------------------------------------------------
float sm_cosd(float angleDegrees) {
  return sm_cos(angleDegrees * DEG2RAD);
}

// --------------------------------------------------------------------------------------
float sm_tand(float angleDegrees) {
  float radians = angleDegrees * DEG2RAD;
  return sm_sin(radians) / sm_cos(radians);
}

// --------------------------------------------------------------------------------------
float sm_asind(float x) {
  return sm_atand(x/sm_sqrt(1.0 - x * x));
}

// --------------------------------------------------------------------------------------
float sm_acosd(float x) {
  return 90.0 - sm_asind(x);
}

// --------------------------------------------------------------------------------------
float sm_atand(float x) {
  return sm_atan(x) * RAD2DEG;
}

// --------------------------------------------------------------------------------------
float sm_sin(float x) {
  int LIMIT = 100;
  int MIN_LIMIT = 6;
  int INVERSE_LIMIT = 100000; // Inverse fractional limit 100000 = 0.00001 == 0.001%
  bool negative = false;
  float xr = x;
  if (x < 0) {
    negative = true;
    xr = -xr;
  }
  // Reduce x to between 0 and pi
  while (xr > TWOPI) {
    xr -= TWOPI;
  }
  if (xr > PI) {
    xr -= PI;
    negative = !negative;
  }
  
  float x1 = xr;// running x, x^3, x^5 etc
  float y = xr; // result
  float yl = 1; // last result
  float yd = 1; // delta result
  float d = 1; // denominator
  bool negIter = false;
  for (int i=1; i<=LIMIT; i++) {
    negIter = !negIter;
    yl = y;
    x1 *= xr*xr;
    d *= i*2;
    d *= i*2+1;
    if (negIter) {
      y -= x1 / d;
    } else {
      y += x1 / d;     
    }
    if (i >= MIN_LIMIT) {
      if (y == yl) {
//        printf("sm_sin iterations: %i", i);
        break;
      }
      yd = y / (y - yl);
      if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
//        printf("sm_sin iterations: %i", i);
        break;
      }
    }
  }
  if (negative) {
    return -y;
  }
  return y;
}

// --------------------------------------------------------------------------------------
float sm_cos(float x) {
  return sm_sin(x + HALFPI);
}

/*
// --------------------------------------------------------------------------------------
float sm_asin(float x) {
  return sm_atan(x/sm_sqrt(1.0 - x * x));
}
*/
// --------------------------------------------------------------------------------------
float sm_asin(float x) {
  int LIMIT = 100;
  int MIN_LIMIT = 6;
  int INVERSE_LIMIT = 100000; // Inverse fractional limit 100000 = 0.00001 == 0.001%
  float xr = x;
  
  float n1 = 1; // 1st numerator
  float x1 = xr;// running x, x^3, x^5 etc (2nd numerator)
  float y = xr; // result
  float yl = 1; // last result
  float yd = 1; // delta result
  float d1 = 1; // 1st denominator
  float d2 = 1; // 2nd denominator
  for (int i=1; i<=LIMIT; i++) {
    yl = y;
    n1 *= i*2 - 1;
    x1 *= xr*xr;
    d1 *= i*2;
    d2 = i*2 + 1;
    y += n1 * x1 / (d1 * d2);     
    if (i >= MIN_LIMIT) {
      if (y == yl) {
//        printf("sm_asin iterations: %i", i);
        break;
      }
      yd = y / (y - yl);
      if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
//        printf("sm_asin iterations: %i", i);
        break;
      }
    }
  }
  return y;
}

// --------------------------------------------------------------------------------------
float sm_acos(float x) {
  return HALFPI - sm_asin(x);
}

// --------------------------------------------------------------------------------------
float sm_atan(float x) {
  
  int LIMIT = 100;
  int MIN_LIMIT = 6;
  int INVERSE_LIMIT = 100000; // Inverse fractional limit 100000 = 0.00001 == 0.001%

/*  
//  Maclaurin series
//
//                x^3   x^5   x^7   x^9
// tan-1(x) = x - ___ + ___ - ___ + ___ - ...
//                 3     5     7     9

// !!! Only works for x >= -1 && x <= 1
    
  float xr = x;
    
  float x1 = xr;// running x, x^3, x^5 etc
  float y = xr; // result
  float yl = 1; // last result
  float yd = 1; // inverse delta result
  float d = 1; // denominator
  bool negIter = false;
  for (int i=1; i<=LIMIT; i++) {
    negIter = !negIter;
    yl = y;
    x1 *= xr*xr;
    d = i*2 + 1;
    if (negIter) {
      y -= x1 / d;
    } else {
      y += x1 / d;     
    }
    if (i >= MIN_LIMIT) {
      if (y == yl) {
        printf("iterations: %i", i);
        break;
      }
      yd = y / (y - yl);
      if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
        printf("iterations: %i", i);
        break;
      }
    }
  }
  return y;
*/
  
/*
//  Euler:
//             __inf  2^2n (n!)^2      x^(2n+1)
//  tan-1(x) = \      ___________  _______________
//             /       (2n + 1)!   (1 + x^2)^(n+1)
//             ¬¬n=0

// Problem is that the numerator and denominator start getting very large - possibly overflow

  float oneplusxsquared = 1 + x*x;
  float y = x / oneplusxsquared;

  double twotothe2n = 1;
  double nfactorial = 1;
  double twonplus1factorial = 1;

  float xtothe2nplus1 = x;
  float oneplusxsquaredtonplus1 = oneplusxsquared;
  float ydelta = 0;
  float yd = 1; // inverse delta result

  for (int n=1; n<=LIMIT; n++) {
    twotothe2n *= 4;
    nfactorial *= n;
    twonplus1factorial *= 2*n;
    twonplus1factorial *= 2*n + 1;
    xtothe2nplus1 *= x*x;
    oneplusxsquaredtonplus1 *= oneplusxsquared;
    ydelta = twotothe2n * nfactorial * nfactorial / twonplus1factorial;
    ydelta *= xtothe2nplus1 / oneplusxsquaredtonplus1;

    y += ydelta;

    if (n >= MIN_LIMIT) {
      if (ydelta == 0.0) {
        printf("iterations: %i", n);
        break;
      }
      yd = y / ydelta;
      if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
        printf("iterations: %i", n);
        break;
      }
    }

  }
  return y;
*/  
  if (x >= -3.0 && x <= 3.0) {
// Castellanos 1988
//
//            y        2 y   2*4 y^2   2*4*6 y^3
// tan-1(x) = _ * (1 + ___ + _______ + _________ + ...)
//            x        3     3*5       3*5*7
//
//             x^2
// where y = _______
//           1 + x^2
//
// Doesn't seem to compute well outside of range -3 <= x <= +3

    float y = 1;
    float yratio = x*x/(1 + x*x);
    float yrationpowern = 1.0;
    float fraction = 1.0;
    float ydelta = 0.0;
    float yd = 1.0;
    
    for (int n=1; n<=LIMIT; n++) {
      yrationpowern *= yratio;
      fraction = fraction * 2.0*n / (2.0*n + 1);
      ydelta = fraction * yrationpowern;
  
      y += ydelta;
  
      if (n >= MIN_LIMIT) {
        if (ydelta == 0.0) {
//          printf("sm_atan iterations: %i", n);
          break;
        }
        yd = y / ydelta;
        if (yd > INVERSE_LIMIT || yd < -INVERSE_LIMIT) {
//          printf("sm_atan iterations: %i", n);
          break;
        }
      }
    }
    y = y * yratio / x;
    return y;
    
  } else {

//  tan-1(x) - sign(x) * cos-1(1/sqrt(x*x + 1))
// 
// Due to lack of precision, this is not good between -1 and 1
//
    if (x == 0.0) return 0.0;
    float y = sm_acos(1.0/(sm_sqrt(x*x + 1)));
    if (x<0) return -y;
    else return y;
    }
}
