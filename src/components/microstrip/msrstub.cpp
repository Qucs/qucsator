/*
 * msrstub.cpp - microstrip radial stub class implementation
 *
 * Copyright (C) 2009 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * $Id$
 *
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "component.h"
#include "substrate.h"
#include "msrstub.h"

using namespace qucs;

msrstub::msrstub () : circuit (1) {
  type = CIR_MSRSTUB;
}

// Returns the microstrip radial stub reactance.
nr_double_t msrstub::calcReactance (nr_double_t r1, nr_double_t r2,
				    nr_double_t phi, nr_double_t er,
				    nr_double_t h, nr_double_t frequency) {

  nr_double_t l0 = C0 / frequency;
  nr_double_t W = (r1 + (r2 - r1) / 2) * (phi);
  nr_double_t ereff = (er + 1.0) / 2 + (er - 1.0) /
    (2.0 * qucs::sqrt (1 + 10.0 * h / W));
  nr_double_t k = 2.0 * pi * qucs::sqrt (ereff) / l0;
  nr_double_t a = k * r1;
  nr_double_t b = k * r2;
  nr_double_t Z_0 = Z0 / qucs::sqrt (ereff) * qucs::sqrt (sqr (j0 (a)) + sqr (y0 (a))) /
    qucs::sqrt (sqr (j1 (a)) + sqr (y1 (a)));
  nr_double_t theta_1 = qucs::atan (y0 (a) / j0 (a));
  //  nr_double_t theta_2 = atan (y0 (b) / j0 (b));
  nr_double_t phi_1 = qucs::atan (-j1 (a) / y1 (a));
  nr_double_t phi_2 = qucs::atan (-j1 (b) / y1 (b));

  nr_double_t X1 = h * Z_0 / (2.0 * pi * r1) * (2.0 * pi / phi) *
    qucs::cos (theta_1 - phi_2) / qucs::sin (phi_1 - phi_2);

  return X1;
}

void msrstub::calcSP (nr_double_t frequency) {
  setS (NODE_1, NODE_1, ztor (calcZ (frequency)));
}

nr_complex_t msrstub::calcZ (nr_double_t frequency) {

  //radial stub - get parameters
  substrate * subst = getSubstrate ();
  nr_double_t r_i = getPropertyDouble ("ri"); //inner radius
  nr_double_t r_o = getPropertyDouble ("ro"); //outer radius
  nr_double_t phi = deg2rad(getPropertyDouble ("alpha")); //stub radius (called phi in formulas, although left alpha in UI)
  nr_double_t w = getPropertyDouble ("Wf"); //feed line width
  nr_double_t h = subst->getPropertyDouble ("h"); //substrate thickness
  nr_double_t t = subst->getPropertyDouble ("t"); //metal thickness
  nr_double_t eps_r = subst->getPropertyDouble ("er"); //substrate relative permittivity
  nr_double_t resist = subst->getPropertyDouble ("rho"); //metal resistivity
  nr_double_t loss_tan = subst->getPropertyDouble ("tand"); //loss tangent of substrate
  const char * EffDimens = getPropertyString ("EffDimens");
  const char * Model = getPropertyString ("Model");

  ////////////////////////////////////////////
  //corrected dimensions of the stub to account for fringing fields:
  nr_double_t p = r_i * qucs::cos(phi/2); //the following intermediate values can be put into conditions to speed up calculations of siple models, like no corrections at all
  nr_double_t u = w / h;
  nr_double_t del_u1 = t / pi * qucs::log(1 + 4*euler / (t * qucs::pow(qucs::coth(qucs::sqrt(6.517*u)), 2.0)));
  nr_double_t w_e = w + h * del_u1 / 2 * (1 + 1 / (qucs::cosh(qucs::sqrt(eps_r - 1))));
  nr_double_t r_s = (2 * p + w_e - w) / (2 * qucs::cos(phi/2));

  nr_double_t u_rad = (phi * r_s) / h;
  nr_double_t del_u1_rad = t / pi * qucs::log(1 + 4*euler / (t * qucs::pow(qucs::coth(qucs::sqrt(6.517*u_rad)), 2.0)));
  nr_double_t w_e_phi_r_s = (phi * r_s) + h * del_u1_rad / 2 * (1 + 1 / (qucs::cosh(qucs::sqrt(eps_r - 1))));
  nr_double_t del_rs = (w_e_phi_r_s - phi * r_s) / 2;

  nr_double_t w_g = (2 * p + w_e - w) * qucs::tan(phi/2);
  nr_double_t w_ge = w_g + 2 * del_rs / qucs::cos(phi/2);

  nr_double_t r_ie;
  nr_double_t r_oe;

  if (!strcmp (EffDimens, "Chew_Kong")) {
    r_ie = r_s + del_rs / qucs::tan(phi/2) + del_rs * qucs::tan(phi/2);
    r_oe = r_o * qucs::sqrt(1 + 2 * h / (pi * eps_r * r_o) * (qucs::log(pi * r_o / (2 * h)) + 1.7726));
  }
  else if (!strcmp (EffDimens, "Giannini")) {
    r_ie = w_e / (2 * qucs::sin(phi/2));
    r_oe = r_o * qucs::sqrt(1 + 2 * h / (pi * r_o) * (qucs::log(pi * r_o / (2 * h)) + 1.7726)) + (w_e - w) / (2 * qucs::sin(phi/2));
  }
  else
  {
    r_ie = r_i;
    r_oe = r_o;
  }


  if (!strcmp (Model, "March")) {
    //complex propagation constant (with losses):
    nr_double_t depth = qucs::sqrt(resist / (pi * frequency * 1 * MU0)); //calculate skin depth (assume mu relative 1)
    if(depth > t)
      depth = t;
    nr_double_t Rsurf = resist / depth; //surface resistivity
    nr_double_t k_beta = 2 * pi * frequency * qucs::sqrt(eps_r) / C0; //real part of complex wavenumber
    nr_double_t k_alpha = Rsurf * qucs::sqrt(eps_r) / (Z0) + k_beta / 2 * loss_tan; //imaginary part of complex wavenumber (losses)
    nr_complex_t k = nr_complex_t(k_beta, k_alpha);

    //impedance calculation 
    nr_double_t Z0_rie = Z0 * h / (r_ie * phi * qucs::sqrt(eps_r));
    nr_complex_t cot_krie_kroe = (yn(0, k * r_ie) * jn(1, k * r_oe) - jn(0, k * r_ie) * yn(1, k * r_oe)) /
    (jn(1, k * r_ie) * yn(1, k * r_oe) - yn(1, k * r_ie) * jn(1, k * r_oe));
    nr_complex_t Zin = nr_complex_t(0, -1) * Z0_rie * cot_krie_kroe;
    return Zin;
  }

  else if (!strcmp (Model, "Giannini")) {

    // Solving equation for k0n:
    // J'0(k0n*r_oe) * Y'0(k0n*r_ie) - J'0(k0n*r_ie) * Y'0(k0n*r_oe) = 0
    // note that J'0(x)=-J1(x)
    int n_max = 20; //how long is the series to calculate - adjust for acceptable accuracy and computation time
    nr_double_t k0n_steps = 2 * pi / (4 * r_o) * 0.002; // estimate step in finding the roots numerically (estimate wavenember of first resonance of lambda/4 stub, step is 0.2%)
    nr_double_t k0n_root[n_max], k0n_x, k0n_x_inc;
    for(int i_root = 0, zeros_cnt = 0; zeros_cnt < n_max; ++i_root)
    {
      k0n_x = i_root * k0n_steps;
      k0n_x_inc = (i_root+1) * k0n_steps;
      if ((jn(1, k0n_x*r_oe) * yn(1, k0n_x*r_ie) - jn(1, k0n_x*r_ie) * yn(1, k0n_x*r_oe)) *
          (jn(1, k0n_x_inc*r_oe) * yn(1, k0n_x_inc*r_ie) - jn(1, k0n_x_inc*r_ie) * yn(1, k0n_x_inc*r_oe)) < 0)
      {
        k0n_root[zeros_cnt] = k0n_x;
        ++zeros_cnt;
      }
    }


    // free space wavenumber
    nr_double_t k = 2.0 * pi * frequency / C0;

    // kg - wavenumber of feed line
    nr_double_t a_feed = 1.0 + 1.0/49.0 * qucs::log((qucs::pow(u, 4.0) + qucs::pow(u/52.0, 2.0)) / (qucs::pow(u, 4.0) + 0.432)) + 1.0/18.7 * qucs::log(1.0 + qucs::pow(u/18.1, 3.0));
    nr_double_t b_feed = 0.564 * qucs::pow((eps_r - 0.9) / (eps_r + 3), 0.053);
    nr_double_t eps_e1_feed_eff = (eps_r + 1.0) / 2.0 + (eps_r - 1.0) / 2.0 * qucs::pow(1.0 + 10.0/u, -a_feed * b_feed);
    nr_double_t eps_e_feed_eff = eps_e1_feed_eff; // TODO, eq 4.6 - skipped the correction (copper thickness?), to clarify, maybe check https://pure.tue.nl/ws/files/46936050/640261-1.pdf
    nr_double_t kg = 2.0 * pi * frequency * qucs::sqrt(eps_e_feed_eff) / C0; // feed line wavenumber (should be complex?)

    // characteristic impedance of the radial stub
    nr_double_t Z0rs = Z0 * h / (r_ie * phi * qucs::sqrt(eps_r));

    // static mode Qt and P coefficients
    nr_double_t Qt0 = 1.0 / loss_tan;
    nr_double_t P0 = qucs::sqrt((2.0 * w_ge) / (phi * (qucs::pow(r_oe, 2.0) - qucs::pow(r_ie, 2.0)))); // is the alpha the same as phi (angle of stub)?


    // impedance and phase velocity in annular and radial directions for capacitances calculations
    // annular:
    nr_double_t w_ann = r_o - r_i;
    nr_double_t del_u1_ann = t / pi * qucs::log(1.0 + 4.0 * euler / (t * qucs::pow((qucs::coth(qucs::sqrt(6.517 * w_ann/h))), 2.0))); //not sure if correction for effective width is needed here, probably doesn't matter much
    nr_double_t we_ann = w_ann + h * del_u1_ann / 2.0 * (1.0 + 1.0 / (qucs::cosh(qucs::sqrt(eps_r - 1))));
    nr_double_t ue_ann = we_ann / h;
    nr_double_t u_ann = w_ann / h; // TODO need separate w/h and weffective/h ?
    nr_double_t f_ann = 6.0 + (2.0 * pi - 6.0) * qucs::exp(-qucs::pow(30.666/ue_ann, 0.7528));
    nr_double_t a_ann = 1.0 + 1.0/49.0 * qucs::log((qucs::pow(u_ann, 4.0) + qucs::pow((u_ann/52), 2.0)) / (qucs::pow(u_ann, 4.0) + 0.432)) + 1/18.7 * qucs::log(1.0 + qucs::pow(u_ann/18.1, 3.0));
    nr_double_t b_ann = 0.564 * qucs::pow((eps_r - 0.9) / (eps_r + 3), 0.053);
    nr_double_t Z01_we_ann = 60.0 * qucs::log(f_ann / ue_ann + qucs::sqrt(1.0 + 4.0 / qucs::pow(ue_ann, 2.0)));
    nr_double_t eps_e1_ann = (eps_r + 1.0) / 2.0 + (eps_r - 1.0) / 2.0 * qucs::pow(1.0 + 10.0/ue_ann, -a_ann*b_ann);
    nr_double_t Z0_ann = Z01_we_ann / qucs::sqrt(eps_e1_ann);
    nr_double_t eps_e_ann = eps_e1_ann; // TODO, eq 4.6 - skipped the correction (copper thickness?), to clarify, maybe check https://pure.tue.nl/ws/files/46936050/640261-1.pdf
    nr_double_t Zphi = Z0_ann;
    nr_double_t vphi = C0 / qucs::sqrt(eps_e_ann);
    // as above but free for free space (eps_r=1)
    nr_double_t we_ann_eps0 = w_ann + h * del_u1_ann;
    nr_double_t ue_ann_eps0 = we_ann_eps0 / h;
    nr_double_t f_ann_eps0 = 6.0 + (2.0 * pi - 6.0) * qucs::exp(-qucs::pow(30.666/ue_ann_eps0, 0.7528));
    nr_double_t a_ann_eps0 = 1.0 + 1.0/49.0 * qucs::log((qucs::pow(ue_ann_eps0, 4.0) + qucs::pow((ue_ann_eps0/52), 2.0)) / (qucs::pow(ue_ann_eps0, 4.0) + 0.432)) + 1/18.7 * qucs::log(1.0 + qucs::pow(ue_ann_eps0/18.1, 3.0));
    nr_double_t b_ann_eps0 = 0.564 * qucs::pow((1 - 0.9) / (1 + 3), 0.053);
    nr_double_t Z01_we_ann_eps0 = 60.0 * qucs::log(f_ann_eps0 / ue_ann_eps0 + qucs::sqrt(1.0 + 4.0 / qucs::pow(ue_ann_eps0, 2.0)));
    nr_double_t eps_e1_ann_eps0 = (1.0 + 1.0) / 2.0 + (1.0 - 1.0) / 2.0 * qucs::pow(1.0 + 10.0/ue_ann_eps0, -a_ann_eps0*b_ann_eps0);
    nr_double_t Z0_ann_eps0 = Z01_we_ann_eps0 / qucs::sqrt(eps_e1_ann_eps0);
    nr_double_t eps_e_ann_eps0 = eps_e1_ann_eps0; // TODO, eq 4.6 - skipped the correction (copper thickness?), to clarify, maybe check https://pure.tue.nl/ws/files/46936050/640261-1.pdf
    nr_double_t Zphi_eps0 = Z0_ann_eps0;
    nr_double_t vphi_eps0 = C0 / qucs::sqrt(eps_e_ann_eps0);
    // radial (calculate average width):
    int int_steps = 100;
    nr_double_t Zr_vr_int_sum = 0;
    nr_double_t step_rad = phi * (r_i + (r_o - r_i) * 2.0 / (int_steps-1.0)) - phi * (r_i + (r_o - r_i) * 1.0 / (int_steps-1.0));
    nr_double_t w_rad, del_u1_rad, we_rad, ue_rad, u_rad, f_rad, a_rad, b_rad, Z01_we_rad, eps_e1_rad, Z0_rad, eps_e_rad, Zrad, vrad;
    for(int i_rad = 0; i_rad < int_steps; ++i_rad)
    {
      w_rad = phi * (r_i + (r_o - r_i) * i_rad / (int_steps-1.0)); // width in radial direction
      del_u1_rad = t / pi * qucs::log(1.0 + 4.0 * euler / (t * qucs::pow((qucs::coth(qucs::sqrt(6.517 * w_rad/h))), 2.0))); //not sure if correction for effective width is needed here, probably doesn't matter much
      we_rad = w_rad + h * del_u1_rad / 2.0 * (1.0 + 1.0 / (qucs::cosh(qucs::sqrt(eps_r - 1))));
      ue_rad = we_rad / h;
      u_rad = w_rad / h; // TODO need separate w/h and weffective/h ?
      f_rad = 6.0 + (2.0 * pi - 6.0) * qucs::exp(-qucs::pow(30.666/ue_rad, 0.7528));
      a_rad = 1.0 + 1.0/49 * qucs::log((qucs::pow(u_rad, 4.0) + qucs::pow((u_rad/52.0), 2.0)) / (qucs::pow(u_rad, 4.0) + 0.432)) + 1/18.7 * qucs::log(1.0 + qucs::pow(u_rad/18.1, 3.0));
      b_rad = 0.564 * qucs::pow((eps_r - 0.9) / (eps_r + 3), 0.053);
      Z01_we_rad = 60.0 * qucs::log(f_rad / ue_rad + qucs::sqrt(1.0 + 4.0 / qucs::pow(ue_rad, 2.0)));
      eps_e1_rad = (eps_r + 1.0) / 2.0 + (eps_r - 1.0) / 2.0 * qucs::pow(1.0 + 10.0/ue_rad, -a_rad*b_rad);
      Z0_rad = Z01_we_rad / qucs::sqrt(eps_e1_rad);
      eps_e_rad = eps_e1_rad; // TODO, eq 4.6 - skipped the correction (copper thickness?), to clarify, maybe check https://pure.tue.nl/ws/files/46936050/640261-1.pdf
      Zrad = Z0_rad;
      vrad = C0 / qucs::sqrt(eps_e_rad);
      Zr_vr_int_sum += 1.0 / (Zrad * vrad) * step_rad;
    }
    nr_double_t ZrVr = 1.0 / (1.0 / (phi * (r_o - r_i)) * Zr_vr_int_sum);
    // as above but free for free air (eps_r=1)
    nr_double_t Zr_vr_int_sum_eps0 = 0;
    nr_double_t w_rad_eps0, del_u1_rad_eps0, we_rad_eps0, ue_rad_eps0, u_rad_eps0, f_rad_eps0, a_rad_eps0, b_rad_eps0, Z01_we_rad_eps0, eps_e1_rad_eps0, Z0_rad_eps0, eps_e_rad_eps0, Zrad_eps0, vrad_eps0;
    for(int i_rad = 0; i_rad < int_steps; ++i_rad)
    {
      w_rad_eps0 = phi * (r_i + (r_o - r_i) * i_rad / (int_steps-1.0)); // width in radial direction
      del_u1_rad_eps0 = t / pi * qucs::log(1.0 + 4.0 * euler / (t * qucs::pow((qucs::coth(qucs::sqrt(6.517 * w_rad_eps0/h))), 2.0))); //not sure if correction for effective width is needed here, probably doesn't matter much
      we_rad_eps0 = w_rad_eps0 + h * del_u1_rad_eps0;
      ue_rad_eps0 = we_rad_eps0 / h;
      u_rad_eps0 = w_rad_eps0 / h; // TODO need separate w/h and weffective/h ?
      f_rad_eps0 = 6.0 + (2.0 * pi - 6.0) * qucs::exp(-qucs::pow(30.666/u_rad_eps0, 0.7528));
      a_rad_eps0 = 1.0 + 1.0/49 * qucs::log((qucs::pow(u_rad_eps0, 4.0) + qucs::pow((u_rad_eps0/52), 2.0)) / (qucs::pow(u_rad_eps0, 4.0) + 0.432)) + 1/18.7 * qucs::log(1.0 + qucs::pow(u_rad_eps0/18.1, 3.0));
      b_rad_eps0 = 0.564 * qucs::pow((1 - 0.9) / (1 + 3), 0.053);
      Z01_we_rad_eps0 = 60.0 * qucs::log(f_rad_eps0 / ue_rad_eps0 + qucs::sqrt(1.0 + 4.0 / qucs::pow(ue_rad_eps0, 2.0)));
      eps_e1_rad_eps0 = (1.0 + 1.0) / 2.0 + (1.0 - 1.0) / 2.0 * qucs::pow(1.0 + 10.0/ue_rad_eps0, -a_rad_eps0*b_rad_eps0);
      Z0_rad_eps0 = Z01_we_rad_eps0 / qucs::sqrt(eps_e1_rad_eps0);
      eps_e_rad_eps0 = eps_e1_rad_eps0; // TODO, eq 4.6 - skipped the correction (copper thickness?), to clarify, maybe check https://pure.tue.nl/ws/files/46936050/640261-1.pdf
      Zrad_eps0 = Z0_rad_eps0;
      vrad_eps0 = C0 / qucs::sqrt(eps_e_rad_eps0);
      Zr_vr_int_sum_eps0 += 1.0 / (Zrad_eps0 * vrad_eps0) * step_rad;
    }
    nr_double_t ZrVr_eps0 = 1.0 / (1.0 / (phi * (r_o - r_i)) * Zr_vr_int_sum_eps0);

    // calculations of capacitances for dynamic permittivity calculations
    nr_double_t Cd0_eps_r = E0 * eps_r * phi * (qucs::pow(r_o, 2.0) - qucs::pow(r_i, 2.0)) / (2.0 * h);
    nr_double_t Cdf_eps_r = (r_o + r_i) * phi / 2.0 * (1.0 / (Zphi * vphi) - E0 * eps_r * (r_o - r_i) / h);
    nr_double_t Cdf_eps0 = (r_o + r_i) * phi / 2.0 * (1.0 / (Zphi_eps0 * vphi_eps0) - E0 * 1 * (r_o - r_i) / h);
    nr_double_t Cds_eps_r = (r_o - r_i) / 2.0 * (1.0 / (ZrVr) - E0 * eps_r * (r_o + r_i) * phi / (2.0 * h));
    nr_double_t Cds_eps0 = (r_o - r_i) / 2.0 * (1.0 / (ZrVr_eps0) - E0 * 1 * (r_o + r_i) * phi / (2.0 * h));
    
    // static mode dynamic permittivity
    nr_double_t eps_d0 = (Cd0_eps_r + Cdf_eps_r + 2 * Cds_eps_r) / (Cd0_eps_r / eps_r + Cdf_eps0 + 2 * Cds_eps0);

    // higher modes
    nr_complex_t sum_part_zin = nr_complex_t(0, 0);
    nr_double_t Kn, A0n, B0n, P0n, Cd0n_eps_r, eps_d0n, omega0n, Rs, K, Itheta, dItheta, Qc0n, Qr0n, Qt0n;
    for(int i_loop = 0; i_loop < n_max; ++i_loop)
    {  
      // coefficients ...
      Kn = -jn(1, k0n_root[i_loop]*r_oe) / yn(1, k0n_root[i_loop]*r_oe);
      A0n = qucs::sqrt(2/phi) / qucs::sqrt(qucs::pow(r_oe, 2.0) * qucs::pow(jn(0, k0n_root[i_loop]*r_oe) + Kn*yn(0, k0n_root[i_loop]*r_oe), 2.0) - qucs::pow(r_ie, 2.0) * qucs::pow(jn(0, k0n_root[i_loop]*r_ie) + Kn*yn(0, k0n_root[i_loop]*r_ie), 2.0));
      B0n = Kn * A0n;
      P0n = qucs::sqrt(w_ge) * (A0n * jn(0, k0n_root[i_loop]*r_ie) + B0n * yn(0, k0n_root[i_loop]*r_ie));
  
      // calculation of capacitance Cd0n and dynamic permittivity (higher modes)
      Cd0n_eps_r = E0 * eps_r * phi / (2.0 * h) * (qucs::pow(r_o, 2.0) - qucs::pow(r_i, 2.0) * qucs::pow((jn(0, k0n_root[i_loop] * r_oe) + Kn * yn(0, k0n_root[i_loop] * r_oe)) / (jn(0, k0n_root[i_loop] * r_ie) + Kn * yn(0, k0n_root[i_loop] * r_ie)), 2.0));
      eps_d0n = (Cd0n_eps_r + Cdf_eps_r + Cds_eps_r) / (Cd0n_eps_r / eps_r + Cdf_eps0 + Cds_eps0);
  
      // quality factors calculations
      omega0n = k0n_root[i_loop] * C0 / qucs::sqrt(eps_d0n); //TODO - correct?
      Rs = qucs::sqrt(omega0n * MU0 / (2.0 * 1.0 / resist)); //TODO - what kind of resistivity? what with skin effect?
      K = 1.0 - 
        (qucs::pow(r_ie, 2.0) * qucs::pow(jn(0, k0n_root[i_loop]*r_ie) * yn(1, k0n_root[i_loop]*r_ie) - jn(1, k0n_root[i_loop]*r_ie) * yn(0, k0n_root[i_loop]*r_ie), 2.0)) /
        (qucs::pow(r_oe, 2.0) * qucs::pow(jn(0, k0n_root[i_loop]*r_oe) * yn(1, k0n_root[i_loop]*r_oe) - jn(1, k0n_root[i_loop]*r_ie) * yn(0, k0n_root[i_loop]*r_oe), 2.0));
      Itheta = 0;
      for(int i_theta = 1; i_theta <= 180; ++i_theta) // integrate numerically (each step 1deg)
      {
        dItheta = 1.0/180.0 * qucs::pow(jn(1, K * r_oe * qucs::sin(deg2rad(i_theta))), 2.0) * qucs::sin(deg2rad(i_theta));
        Itheta += dItheta;
      }
      Qc0n = omega0n * MU0 * h / (2.0 * Rs);
      Qr0n = 1.0 / (omega0n * h * Itheta / (eps_d0 * K * C0));
      Qt0n = 1.0 / (1.0 / Qt0 + 1.0 / Qc0n + 1.0 / Qr0n);
      
      // part of stub impedance (higher modes)
      sum_part_zin += (kg * qucs::pow(P0n, 2.0)) / (qucs::pow(k0n_root[i_loop], 2.0) * (1.0 + nr_complex_t(0,1) / Qt0n) - qucs::pow(k, 2.0) * eps_d0n);
    }

    //normalized impedance and impedance calculations
    nr_complex_t z_in = (1.0 / Qt0 - nr_complex_t(0,1) * kg * qucs::pow(P0, 2.0)) / (qucs::pow(k, 2.0) * eps_d0) + nr_complex_t(0,1) * sum_part_zin;
    nr_complex_t Zin = Z0rs * z_in;
    return Zin;
  }

  else {
    return nr_complex_t (0, calcReactance (r_ie, r_oe, phi, eps_r, h, frequency)); //old model used in qucs
  }
}

void msrstub::initDC (void) {
  allocMatrixMNA ();
  setY (NODE_1, NODE_1, 0);
}

void msrstub::calcAC (nr_double_t frequency) {
  setY (NODE_1, NODE_1, 1.0 / calcZ (frequency));
}

// properties
PROP_REQ [] = {
  { "ri", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "ro", PROP_REAL, { 10e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "Wf", PROP_REAL, { 1e-3, PROP_NO_STR }, PROP_POS_RANGE },
  { "alpha", PROP_REAL, { 90, PROP_NO_STR }, PROP_RNGII (0, 180) },
  { "Subst", PROP_STR, { PROP_NO_VAL, "Subst1" }, PROP_NO_RANGE },
  { "EffDimens", PROP_STR, { PROP_NO_VAL, "OldQucsNoCorrection" }, PROP_RNG_MRSCORR },
  { "Model", PROP_STR, { PROP_NO_VAL, "OldQucsModel" }, PROP_RNG_MRSMOD },
  PROP_NO_PROP };
PROP_OPT [] = {
  PROP_NO_PROP };
struct define_t msrstub::cirdef =
  { "MRSTUB", 1, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR, PROP_DEF };


// Literature
// Design of a rectenna for wireless low-power transmission, Akkermans, J.A.G.
// A New Simple and Accurate Formula for Microstrip Radial Stub, Roberto Sorrentino, Luca Roselli
// A Crooked U-Slot Dual-Band Antenna With RadialStub Feeding, Hyo Rim Bae, Soon One So, Choon Sik Cho, Member, IEEE, Jae W. Lee, Member, IEEE, and Jaeheung Kim, Member, IEEE
// Performance Characterization of Radial Stub Microstrip Bow-Tie Antenna, B.T.P.Madhav, 2S.S.Mohan Reddy, 3Neha Sharma, 3J. Ravindranath Chowdary 3Bala Rama Pavithra, 3K.N.V.S. Kishore, 3G. Sriram, 3B. Sachin Kumar 1Associate Professor, Department of ECE, K L University, Guntur DT, AP, India 2Associate Professor, Department of ECE, SRKR Engineering College, Bhimavaram, AP, India 3Project Students, Department of ECE, K L University, Guntur DT, AP, India, Performance Characterization of Radial Stub Microstrip Bow-Tie Antenna. Available from: https://www.researchgate.net/publication/286889155_Performance_Characterization_of_Radial_Stub_Microstrip_Bow-Tie_Antenna [accessed Oct 19 2021].
