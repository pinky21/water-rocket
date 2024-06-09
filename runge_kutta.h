#ifndef _RUNGEKUTTA_H_
#define _RUNGEKUTTA_H_

#include <iostream>
#include <math.h>
#include <vector>

namespace runge_kutta {

static const double SAFETY = 0.9;
static const double EXPSHRINK = -0.25;
static const double EXPGROW = -0.2;

void rk45_cashkarp(const std::vector<double>& y,
            const int n,
            const double x,
            const double h,
            std::vector<double>& yout,
            std::vector<double>& yerr,
            void (*dxdy)(const double, const std::vector<double>&, std::vector<double>&)) {

  static const double a2 = 0.2;
  static const double a3 = 0.3;
  static const double a4 = 0.6; 
  static const double a5 = 1.0; 
  static const double a6 = 0.875;
  static const double b21 = 0.2;
  static const double b31 = 3.0 / 40.0;
  static const double b32 = 9.0 / 40.0;
  static const double b41 = 0.3;
  static const double b42 = -0.9;
  static const double b43 = 1.2;
  static const double b51 = -11.0 / 54.0;
  static const double b52 = 2.5;
  static const double b53 = -70.0 / 27.0;
  static const double b54 = 35.0 / 27.0;
  static const double b61 = 1631.0 / 55296.0;
  static const double b62 = 175.0 / 512.0;
  static const double b63 = 575.0 / 13824.0;
  static const double b64 = 44275.0 / 110592.0;
  static const double b65 = 253.0 / 4096.0;
  static const double c1 = 37.0 / 378.0;
  static const double c3 = 250.0 / 621.0;
  static const double c4 = 125.0 / 594.0;
  static const double c6 = 512.0 / 1771.0;
  static const double dc1 = c1 - 2825.0 / 27648.0;
  static const double dc3 = c3 - 18575.0 / 48384.0;
  static const double dc4 = c4 - 13525.0 / 55296.0;
  static const double dc5 = -277.00 / 14336.0;
  static const double dc6 = c6 - 0.25;

  std::vector<double> ak2(n);
  std::vector<double> ak3(n);
  std::vector<double> ak4(n);
  std::vector<double> ak5(n);
  std::vector<double> ak6(n);
  std::vector<double> dydx(n);
  std::vector<double> ytemp(n);

  (*dxdy)(x, y, dydx);
  for (int i = 0; i < n; i++) {
    ytemp[i] = y[i] + b21 * h * dydx[i];
  }
  (*dxdy)(x + a2 * h, ytemp, ak2);
  for (int i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
  }
  (*dxdy)(x + a3 * h, ytemp, ak3);
  for (int i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
  }
  (*dxdy)(x + a4 * h, ytemp, ak4);
  for (int i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
  }
  (*dxdy)(x + a5 * h, ytemp, ak5);
  for (int i = 0; i < n; i++) {
    ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
  }
  (*dxdy)(x + a6 * h, ytemp, ak6);
  for (int i = 0; i < n; i++) {
    yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
  }
  for (int i = 0; i < n; i++) {
    yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
  }
}

void rk_driver(std::vector<double>& y,
          int n,
          double *x,
          double h_try,
          double eps,
          const std::vector<double>& yscal,
          double *h_next,
          void (*dydx)(double, const std::vector<double>&, std::vector<double>&)) {

  double errmax, h, htemp, xnew;
  std::vector<double> yerr(n);
  std::vector<double> yresult(n);
  h = h_try; // Set stepsize to the initial trial value.
  
  while(true) {
    rk45_cashkarp(y, n, *x, h, yresult, yerr, dydx);
    errmax = 0.0;
    for (int i = 0; i < n; i++) {
      errmax = fmax(errmax, fabs(yerr[i] / yscal[i]));
    }
    errmax /= eps; // Scale relative to required tolerance.
    if (errmax <= 1.0) {
      break; // Step succeeded. Compute size of next step.
    }
    htemp = SAFETY * h * pow(errmax, EXPSHRINK);
    // Truncation error too large, reduce stepsize.
    // But not more than a factor of 10.
    if (h >= 0.0) {
      h = fmax(htemp, 0.1 * h);
    }
    else {
      h = fmin(htemp, 0.1 * h);
    }
    xnew = (*x) + h;
    if (xnew == *x) {
      std::cerr << "stepsize underflow in rk_driver\n";
      exit(1);
    }
  }

  *h_next = SAFETY * h * fmin(pow(errmax, EXPGROW), 5.0);
  *x += h;
  for (int i = 0; i < n; i++) {
    y[i] = yresult[i];
  }
}

}

#endif // _RUNGEKUTTA_H_