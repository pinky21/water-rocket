#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>

#include "runge_kutta.h"

namespace {

  static const double HWATERMIN = 1.e-15;  // minimum water height in bottle.

  // RK parameters
  static const double TINY = 1.e-25;
  static const int nvar = 3;
  static const double eps = 1.e-12;
  static const double t_max = 200.;

  // constants
  static const double rho_water = 1000.;  // kg/m^3
  static const double g = 9.81;  // m/s^2
  static const double p_ambient = 100000.;  // Pa
  static const double rho_ambient = 1.22;  // kg/m^3
  static const double T_ambient = 285.547439;  // K, chosen so that R_s * T * rho_A0 = p_atm
  static const double R_s = 287.052874;  // J/(kg K) specific Gas constant for air
  static const double poly = 1.4;  // dimensionless polytropic index of air

  struct Rocket {
  
    // rocket parameters
    const double m_hull = 0.05;  // kg
    const double c_D = 0.82;  // dimensionless drag coeff
    const double A_expel = 0.000314;  // m^2
    const double A_rocket = 0.005;  // m^2
    const double L_rocket = 0.3;  // m
    double p_air_init = 700000.;  // Pa
    double h_water_init = 0.1;  // m (must be less than L_rocket)

  };

  Rocket rocket;
  double h_rocket_max = 0.;
}


double calc_m_rocket(double h_water) {
  if (h_water <= HWATERMIN) {
    return rocket.m_hull;
  }
  return rocket.m_hull + rho_water * rocket.A_rocket * h_water;

}


double calc_p_air(double h_water) {
  return rocket.p_air_init * pow((rocket.L_rocket - rocket.h_water_init) / (rocket.L_rocket - h_water), poly);
}


double calc_v_expel_water(double h_water) {

  // check if the pressure in the bottle is lower or equal atm pressure
  // 1) if so, the exausting velocity is zero
  double p_air = calc_p_air(h_water);
  double delta_p = p_air - p_ambient;
  if (delta_p <= 0.) {
    return 0.0;
  }

  // 2) if not, the finite exausting velocity can be calculated
  return sqrt(2. / rho_water / (1. - (rocket.A_expel / rocket.A_rocket) * (rocket.A_expel / rocket.A_rocket)) * delta_p);
}


double calc_acceleration_drag(double m_rocket, double v_rocket) {
  return 0.5 * rho_ambient / m_rocket * rocket.c_D * rocket.A_rocket * v_rocket * fabs(v_rocket);
}


void deriv_water_thrust(double x, const std::vector<double>& y, std::vector<double>& dydx) {

  double m_rocket = calc_m_rocket(y[2]);

  // acceleration due to drag
  double a_drag = calc_acceleration_drag(m_rocket, y[0]);

  // acceleration due to thrust
  double v_expel = calc_v_expel_water(y[2]);
  double a_thrust = rocket.A_expel * rho_water * v_expel * v_expel / m_rocket;

  // derivatives
  dydx[0] = a_thrust - g - a_drag;
  dydx[1] = y[0];
  dydx[2] = -rocket.A_expel / rocket.A_rocket * v_expel;

}


double calc_v_expel_air(double rho_air) {
  
  // check if the pressure in the bottle is lower or equal atm pressure
  // 1) if so, the exausting velocity is zero
  double p_air = rho_air * R_s * T_ambient;
  double delta_p = p_air - p_ambient;
  if (delta_p <= 0.) {
    return 0.0;
  }
  
  // 2) if not, the finite exausting velocity can be calculated
  return sqrt(2. / rho_air * delta_p);
}


void deriv_air_thrust(double x, const std::vector<double>& y, std::vector<double>& dydx) {
  
  double m_rocket = rocket.m_hull + y[2] * rocket.A_rocket * rocket.L_rocket;
  
  // acceleration due to drag
  double a_drag = calc_acceleration_drag(m_rocket, y[0]);

  // acceleration due to thrust
  double v_expel = calc_v_expel_air(y[2]);
  double a_thrust = rocket.A_expel * y[2] * v_expel * v_expel / m_rocket;

  // derivatives
  dydx[0] = a_thrust - g - a_drag;
  dydx[1] = y[0];
  dydx[2] = -y[2] * v_expel * rocket.A_expel / rocket.A_rocket / rocket.L_rocket;

}


void deriv_ballistic(double x, const std::vector<double>& y, std::vector<double>& dydx) {

  double m_rocket = calc_m_rocket(y[2]);
  
  // acceleration due to drag
  double a_drag = 0.5 * rho_ambient / m_rocket * rocket.c_D * rocket.A_rocket * y[0] * fabs(y[0]);

  // derivatives
  dydx[0] = -g - a_drag;
  dydx[1] = y[0];
  dydx[2] = 0.;

}


void write_traj(bool header,
                const double t, 
                const double v_rocket, 
                const double h_rocket, 
                const double h_water) {
  
  if (header) {
    std::ofstream out("trajectory.csv",std::ios::out);
    out << "t;v_rocket;h_rocket;h_water" << std::endl;
    out.close();
  }
  else {
    std::ofstream out("trajectory.csv",std::ios::app);
    out << t << ";" << v_rocket << ";" << h_rocket << ";" << h_water << std::endl;
    out.close();
  }
}


void rk_step(double* x, double* h, double x_max, std::vector<double>& y,
              std::vector<double>& dydx, std::vector<double>& yscal,
              void (*derivs)(double, const std::vector<double>&, std::vector<double>&)) {

  double hnext;

  (*derivs)(*x, y, dydx);
  for (int i = 0; i < nvar; i++) {
    yscal[i] = fabs(y[i]) + fabs(dydx[i] * (*h)) + TINY;
  }

  if ((*x) + (*h) > x_max) {
    *h = x_max - (*x);  // Adjust step size if we are overshooting x_max
  }
  runge_kutta::rk_driver(y, nvar, x, *h, eps, yscal, &hnext, derivs);
  *h = hnext;

}


void calculate_trajectory(Rocket rocket, bool traj=true) {

  // the trajectory is divided into three stages: 1) water thrust, 2) air thrust, 3) ballistic motion

  // step size and running variable
  double h = 2.e-6;
  double t = 0.;

  // y[0]: v_rocket 
  // y[1]: h_rocket
  // y[2]: h_water or rho_air for water/air thrust
  // dydt[0]: dv_rocket/dt
  // dydt[1]: dh_rocket/dt
  // dydt[2]: dh_water/dt or drho_air/dt for water/air thrust
  std::vector<double> y{0.0, 0.0, rocket.h_water_init};
  std::vector<double> dydt(nvar);
  std::vector<double> yscal(nvar);

  if (traj) {
    write_traj(true, t, y[0], y[1], y[2]);
    write_traj(false, t, y[0], y[1], y[2]);
  }

  double p_air = rocket.p_air_init;
  double delta_p = p_air - p_ambient;

  // water thrust if any
  if (rocket.h_water_init > HWATERMIN) {
    if (traj) {
      std::cout << "Water thrust starting at t=" << t << "s." << std::endl;
    }
    while (t < t_max) {

      p_air = calc_p_air(y[2]);
      delta_p = p_air - p_ambient;

      if (y[2] <= HWATERMIN || delta_p <= 0.) {
        break;
      }

      rk_step(&t, &h, t_max, y, dydt, yscal, deriv_water_thrust);
      
      if (traj) {
        write_traj(false, t, y[0], y[1], y[2]);
      }
    }
  }

  // air thrust if any
  double water_level = y[2];  // store water level for ballistic motion later
  if (y[2] <= HWATERMIN && delta_p > 0.) {
    if (traj) {
      std::cout << "Air thrust starting at t=" << t << "s." << std::endl;
    }
    y[2] = p_air / (T_ambient * R_s);  // use third variable as rho_air now
    while (t < t_max) {
      p_air = y[2] * T_ambient * R_s;
      delta_p = p_air - p_ambient;

      if (delta_p <= 0.) {
        break;
      }

      rk_step(&t, &h, t_max, y, dydt, yscal, deriv_air_thrust);

      if (traj) {
        write_traj(false, t, y[0], y[1], y[2]);
      }
    }
  }

  // Ballistic motion
  if (traj) {
    std::cout << "Ballisitc motion starting at t=" << t << "s." << std::endl;
  }
  y[2] = water_level;  // use third variable again as water_level
  while (t < t_max) {

    if (y[1] <= 0.) {
      break;
    }

    rk_step(&t, &h, t_max, y, dydt, yscal, deriv_ballistic);

    if (traj) {
      write_traj(false, t, y[0], y[1], y[2]);
    }
    
    // assign peak value
    if (y[1] > h_rocket_max) {
      h_rocket_max = y[1];
    }
  }

}


void write_MC(bool header,
              double p_air_ratio_init,
              double h_water_ratio_init,
              double h_rocket_max) {
  
  if (header) {
    std::ofstream out("max_heights.csv",std::ios::out);
    out << "p_air_ratio_init;h_water_ratio_init;h_rocket_max" << std::endl;
    out.close();
  }
  else {
    std::ofstream out("max_heights.csv",std::ios::app);
    out << p_air_ratio_init << ";" << h_water_ratio_init << ";" << h_rocket_max << std::endl;
    out.close();
  }

}


int main() {

  // calculate a trajectory for the standard parameters specified above
  calculate_trajectory(rocket, true);


  // apply very simple Monte Carlo approach to see influence of pressure and water filling level
  write_MC(true, -1., -1., -1.);
  int n_MC = 100000;
  for (int k = 0; k < n_MC; k++) {
    h_rocket_max = 0.;
    double rn_1 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    double rn_2 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    rocket.p_air_init = p_ambient + 100000. * 7. * rn_1;
    rocket.h_water_init = rn_2 * rocket.L_rocket;
    calculate_trajectory(rocket, false);
    write_MC(false, (rocket.p_air_init - p_ambient) / p_ambient, rocket.h_water_init / rocket.L_rocket, h_rocket_max);
  }

  return 0;
}
