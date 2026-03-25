#include <iostream>
#include <cmath>
using namespace std;

// constant
const double g  = 9.8;
const double rho = 1.9;
const double Cd = 0.4;
const double r = 0.0366;
const double m = 0.145;

const double A = M_PI * r * r;
const double k = (0.5 * Cd * rho * A) / m;

//state: [x, vx, y, vy]
struct State {
    double x, vx, y, vy;
};

// ODE system
State derivatives(const State& s) {
    double v = sqrt(s.vx * s.vx + s.vy * s.vy);

    State ds;
    ds.x = s.vx;
    ds.vx = -k * v * s.vy;
    ds.y = s.vy;
    ds.vy = -g -k * v * s.vy;

    return ds;
}

// RK4 step
State rk4_step(State s, double dt) {
    State k1 = derivatives(s);

    State s2 = {s.x + 0.5*dt*k1.x, s.vx + 0.5*dt*k1.vx,
                s.y + 0.5*dt*k1.y, s.vy + 0.5*dt*k1.vy};
    State k2 = derivatives(s2);

    State s3 = {s.x + 0.5*dt*k2.x, s.vx + 0.5*dt*k2.vx,
                s.y + 0.5*dt*k2.y, s.vy + 0.5*dt*k2.vy};
    State k3 = derivatives(s3);

    State s4 = {s.x + dt*k3.x, s.vx + dt*k3.vx,
                s.y + dt*k3.y, s.vy + dt*k3.vy};
    State k4 = derivatives(s4);

    State result;
    result.x  = s.x  + dt*(k1.x  + 2*k2.x  + 2*k3.x  + k4.x)/6;
    result.vx = s.vx + dt*(k1.vx + 2*k2.vx + 2*k3.vx + k4.vx)/6;
    result.y  = s.y  + dt*(k1.y  + 2*k2.y  + 2*k3.y  + k4.y)/6;
    result.vy = s.vy + dt*(k1.vy + 2*k2.vy + 2*k3.vy + k4.vy)/6;

    return result;
}

// Simulate until time T
State simulate(double v, double theta, double T) {
    State s = {0, v*cos(theta), 0, v*sin(theta)};
    double dt = 0.001;

    for (double t = 0; t < T; t += dt) {
        s = rk4_step(s, dt);
    }
    return s;
}

int main() {
    cout << "k = " << k << endl;

    // --- Catcher ---
    double theta_h = 60 * M_PI / 180.0;
    double T_h = 77.0 / 30.0;

    double v_h = 30.0;

    // Simple search instead of fsolve
    for (int i = 0; i < 100; i++) {
        State s = simulate(v_h, theta_h, T_h);
        double error = s.y;

        v_h -= error * 0.1; // adjust step
    }

    State final_h = simulate(v_h, theta_h, T_h);
    double D = final_h.x;

    cout << "Catcher v_h = " << v_h << ", D = " << D << endl;

    // --- Thrower ---
    double T_w = 21.0 / 30.0;
    double v_w = 20.0;
    double theta_w = 0.5;

    for (int i = 0; i < 200; i++) {
        State s = simulate(v_w, theta_w, T_w);

        double err_x = s.x - D;
        double err_y = s.y;

        v_w -= err_x * 0.01;
        theta_w -= err_y * 0.01;
    }

    cout << "Thrower v_w = " << v_w 
         << ", theta_w (deg) = " << theta_w * 180.0 / M_PI << endl;

    return 0;
}