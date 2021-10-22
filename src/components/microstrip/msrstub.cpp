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
  nr_double_t w = getPropertyDouble ("Wf"); //feed line width needed for calculation of effective dimensions
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
