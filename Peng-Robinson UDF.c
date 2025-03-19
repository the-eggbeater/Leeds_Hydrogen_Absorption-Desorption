#include "udf.h"

#define R_GAS 8314.4621    // J/kmol·K
#define TC_H2 33.15        // K 
#define PC_H2 13.0e5       // Pa 
#define OMEGA_H2 -0.215    
#define HENRY_H2 7.12e6    // Pa·m³/mol
#define M_PI 3.14159265358979323846

// Cubic equation solver (C89-compliant)
int solve_cubic(double a3, double a2, double a1, double a0, double roots[3]) {
    double Q, R, Q3, R2, theta, A, B;
    int i;

    Q = (3 * a1 - a2 * a2) / 9.0;
    R = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54.0;
    Q3 = Q * Q * Q;
    R2 = R * R;

    if (R2 < Q3) {
        theta = acos(R / sqrt(Q3));
        roots[0] = -2 * sqrt(Q) * cos(theta / 3.0) - a2 / 3.0;
        roots[1] = -2 * sqrt(Q) * cos((theta + 2 * M_PI) / 3.0) - a2 / 3.0;
        roots[2] = -2 * sqrt(Q) * cos((theta - 2 * M_PI) / 3.0) - a2 / 3.0;
        return 3;
    }
    else {
        double sign_R = (R >= 0) ? 1.0 : -1.0;
        A = -sign_R * pow(fabs(R) + sqrt(R2 - Q3), 1.0 / 3.0);
        B = (A != 0) ? Q / A : 0;
        roots[0] = (A + B) - a2 / 3.0;
        return 1;
    }
}

real compute_phi(real P, real T, real y_H2) {
    real kappa, alpha, a, b, A, B, Z, sqrt2, term, ln_phi;
    double coeff[4];
    double roots[3];
    int num_roots, i;

    kappa = 0.37464 + 1.54226 * OMEGA_H2 - 0.26992 * OMEGA_H2 * OMEGA_H2;
    alpha = pow(1 + kappa * (1 - sqrt(T / TC_H2)), 2);
    a = 0.45724 * (R_GAS * R_GAS * TC_H2 * TC_H2) / PC_H2 * alpha;
    b = 0.07780 * R_GAS * TC_H2 / PC_H2;

    A = (a * P) / (R_GAS * R_GAS * T * T);
    B = (b * P) / (R_GAS * T);

    coeff[0] = -(A * B - B * B - B * B * B);
    coeff[1] = A - 3 * B * B - 2 * B;
    coeff[2] = -(1 - B);
    coeff[3] = 1.0;

    num_roots = solve_cubic(coeff[3], coeff[2], coeff[1], coeff[0], roots);

    Z = roots[0];
    for (i = 1; i < num_roots; i++) {
        if (roots[i] > Z) Z = roots[i];
    }

    sqrt2 = sqrt(2.0);
    term = log((Z + (1 + sqrt2) * B) / (Z + (1 - sqrt2) * B));
    ln_phi = (Z - 1) - log(Z - B) - (A / (2 * sqrt2 * B)) * term;
    return exp(ln_phi);
}

DEFINE_LINEARIZED_MASS_TRANSFER(h2_mass_transfer, cell, thread,
    from_index, from_species_index,
    to_index, to_species_index,
    lin_from, lin_to)
{
    real m_transfer = 0.0, K = 1e-5, area_vol = 1e3;
    Thread* gas_th = THREAD_SUB_THREAD(thread, to_index);
    Thread* liq_th = THREAD_SUB_THREAD(thread, from_index);
    real T = C_T(cell, gas_th);
    real P_gas = C_P(cell, gas_th);
    real y_H2 = C_YI(cell, gas_th, to_species_index);
    real phi = compute_phi(P_gas, T, y_H2);
    real ff_gas = phi * y_H2 * P_gas;
    real mw_H2 = 2.016e-3, mw_H2O = 18.015e-3;
    real mass_frac_H2 = C_YI(cell, liq_th, from_species_index);
    real x_H2 = (mass_frac_H2 / mw_H2) / (mass_frac_H2 / mw_H2 + (1 - mass_frac_H2) / mw_H2O);
    real f_liq = HENRY_H2 * x_H2;

    m_transfer = K * area_vol * (ff_gas - f_liq);
    *lin_from = -K * area_vol * HENRY_H2 * (1.0 / mw_H2) / pow((mass_frac_H2 / mw_H2 + (1 - mass_frac_H2) / mw_H2O), 2);
    *lin_to = K * area_vol * phi * P_gas;

    return m_transfer;
}




