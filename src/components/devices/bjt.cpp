/*
 * bjt.cpp - bipolar junction transistor class implementation
 *
 * Copyright (C) 2004 Stefan Jahn <stefan@lkcc.org>
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
 * the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.  
 *
 * $Id: bjt.cpp,v 1.21 2004-10-25 21:01:32 ela Exp $
 *
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "complex.h"
#include "matrix.h"
#include "object.h"
#include "node.h"
#include "circuit.h"
#include "net.h"
#include "component_id.h"
#include "constants.h"
#include "device.h"
#include "bjt.h"

#define NEWSGP 1

#define NODE_B 1 /* base node       */
#define NODE_C 2 /* collector node  */
#define NODE_E 3 /* emitter node    */
#define NODE_S 4 /* substrate node  */

bjt::bjt () : circuit (4) {
  cbcx = rb = re = rc = NULL;
  type = CIR_BJT;
}

void bjt::calcSP (nr_double_t frequency) {
  setMatrixS (ytos (calcMatrixY (frequency)));
}

matrix bjt::calcMatrixY (nr_double_t frequency) {

  // fetch computed operating points
  nr_double_t Cbe  = getOperatingPoint ("Cbe");
  nr_double_t gbe  = getOperatingPoint ("gpi");
  nr_double_t Cbci = getOperatingPoint ("Cbci");
  nr_double_t gbc  = getOperatingPoint ("gmu");
  nr_double_t Ccs  = getOperatingPoint ("Ccs");
  nr_double_t gmfr = getOperatingPoint ("gmf");
  nr_double_t gmr  = getOperatingPoint ("gmr");
  nr_double_t Ptf  = getPropertyDouble ("Ptf");
  nr_double_t Tf   = getPropertyDouble ("Tf");

  // compute admittance matrix entries
  complex Ybe = rect (gbe, 2.0 * M_PI * frequency * Cbe);
  complex Ybc = rect (gbc, 2.0 * M_PI * frequency * Cbci);
  complex Ycs = rect (0.0, 2.0 * M_PI * frequency * Ccs);

  // compute influence of excess pase
  nr_double_t phase = rad (Ptf) * Tf * 2 * M_PI * frequency;
  complex gmf = polar (gmfr, -phase);

  // build admittance matrix and convert it to S-parameter matrix
  matrix y (4);
  y.set (NODE_B, NODE_B, Ybc + Ybe);
  y.set (NODE_B, NODE_C, -Ybc);
  y.set (NODE_B, NODE_E, -Ybe);
  y.set (NODE_B, NODE_S, 0);
  y.set (NODE_C, NODE_B, -Ybc + gmf - gmr);
  y.set (NODE_C, NODE_C, Ybc + gmr + Ycs);
  y.set (NODE_C, NODE_E, -gmf);
  y.set (NODE_C, NODE_S, -Ycs);
  y.set (NODE_E, NODE_B, -Ybe - gmf + gmr);
  y.set (NODE_E, NODE_C, -gmr);
  y.set (NODE_E, NODE_E, Ybe + gmf);
  y.set (NODE_E, NODE_S, 0);
  y.set (NODE_S, NODE_B, 0);
  y.set (NODE_S, NODE_C, -Ycs);
  y.set (NODE_S, NODE_E, 0);
  y.set (NODE_S, NODE_S, Ycs);

  return y;
}

void bjt::calcNoise (nr_double_t frequency) {

  // fetch computed operating points
  nr_double_t Ibe = getOperatingPoint ("Ibe");
  nr_double_t Ice = getOperatingPoint ("Ice");

  // get model properties
  nr_double_t Kf  = getPropertyDouble ("Kf");
  nr_double_t Af  = getPropertyDouble ("Af");
  nr_double_t Ffe = getPropertyDouble ("Ffe");
  nr_double_t Kb  = getPropertyDouble ("Kb");
  nr_double_t Ab  = getPropertyDouble ("Ab");
  nr_double_t Fb  = getPropertyDouble ("Fb");

  nr_double_t ib = 2 * Ibe * QoverkB / T0 +            // shot noise
    (Kf * pow (Ibe, Af) / pow (frequency, Ffe) +       // flicker noise
     Kb * pow (Ibe, Ab) / (1 + sqr (frequency / Fb)))  // burst noise
    / kB / T0;
  nr_double_t ic = 2 * Ice * QoverkB / T0;             // shot noise

  /* build noise current correlation matrix and convert it to
     noise-wave correlation matrix */
  matrix y = matrix (4);
  y.set (NODE_B, NODE_B, ib);
  y.set (NODE_B, NODE_E, -ib);
  y.set (NODE_C, NODE_C, ic);
  y.set (NODE_C, NODE_E, -ic);
  y.set (NODE_E, NODE_B, -ib);
  y.set (NODE_E, NODE_C, -ic);
  y.set (NODE_E, NODE_E, ic + ib);
  setMatrixN (cytocs (y * z0, getMatrixS ()));
}

void bjt::initDC (void) {

  // apply polarity of BJT
  char * type = getPropertyString ("Type");
  pol = !strcmp (type, "pnp") ? -1 : 1;

  // get simulation temperature
  nr_double_t T = getPropertyDouble ("Temp");

  // initialize starting values
  UbePrev = real (getV (NODE_B) - getV (NODE_E)) * pol;
  UbcPrev = real (getV (NODE_B) - getV (NODE_C)) * pol;

  // disable additional base-collector capacitance
  if (deviceEnabled (cbcx)) {
    disableCapacitance (this, cbcx, getNet ());
  }

  // possibly insert series resistance at emitter
  nr_double_t Re = getPropertyDouble ("Re");
  if (Re != 0.0) {
    // create additional circuit if necessary and reassign nodes
    re = splitResistance (this, re, getNet (), "Re", "emitter", NODE_E);
    re->setProperty ("R", Re);
    re->setProperty ("Temp", T);
  }
  // no series resistance at emitter
  else {
    disableResistance (this, re, getNet (), NODE_E);
  }

  // possibly insert series resistance at collector
  nr_double_t Rc = getPropertyDouble ("Rc");
  if (Rc != 0.0) {
    // create additional circuit if necessary and reassign nodes
    rc = splitResistance (this, rc, getNet (), "Rc", "collector", NODE_C);
    rc->setProperty ("R", Rc);
    rc->setProperty ("Temp", T);
  }
  // no series resistance at collector
  else {
    disableResistance (this, rc, getNet (), NODE_C);
  }

  // possibly insert base series resistance
  nr_double_t Rb  = getPropertyDouble ("Rb");
  nr_double_t Rbm = getPropertyDouble ("Rbm");
  if (Rbm <= 0.0) Rbm = Rb; // Rbm defaults to Rb if zero
  if (Rb < Rbm)   Rbm = Rb; // Rbm must be less or equal Rb
  setProperty ("Rbm", Rbm);
  if (Rbm != 0.0) {
    // create additional circuit and reassign nodes
    rb = splitResistance (this, rb, getNet (), "Rbb", "base", NODE_B);
    rb->setProperty ("R", Rb);
    rb->setProperty ("Temp", T);
  }
  // no series resistance at base
  else {
    disableResistance (this, rb, getNet (), NODE_B);
    Rbb = 0.0;                 // set this operating point
    setProperty ("Xcjc", 1.0); // other than 1 is senseless here
  }
}

void bjt::calcDC (void) {

  // fetch device model parameters
  nr_double_t Is   = getPropertyDouble ("Is");
  nr_double_t Nf   = getPropertyDouble ("Nf");
  nr_double_t Nr   = getPropertyDouble ("Nr");
  nr_double_t Vaf  = getPropertyDouble ("Vaf");
  nr_double_t Var  = getPropertyDouble ("Var");
  nr_double_t Ikf  = getPropertyDouble ("Ikf");
  nr_double_t Ikr  = getPropertyDouble ("Ikr");
  nr_double_t Bf   = getPropertyDouble ("Bf");
  nr_double_t Br   = getPropertyDouble ("Br");
  nr_double_t Ise  = getPropertyDouble ("Ise");
  nr_double_t Isc  = getPropertyDouble ("Isc");
  nr_double_t Ne   = getPropertyDouble ("Ne");
  nr_double_t Nc   = getPropertyDouble ("Nc");
  nr_double_t Rb   = getPropertyDouble ("Rb");
  nr_double_t Rbm  = getPropertyDouble ("Rbm");
  nr_double_t Irb  = getPropertyDouble ("Irb");
  nr_double_t T    = getPropertyDouble ("Temp");

  nr_double_t Ut, Ube, Ubc, Q1, Q2;
  nr_double_t Iben, Ibcn, Ibei, Ibci, Ibc, gbe, gbc, gtiny;
  nr_double_t Uce, IeqB, IeqC, IeqE, IeqS, UbeCrit, UbcCrit;
  nr_double_t gm, go;

  // interpret zero as infinity for these model parameters
  Ikf = Ikf > 0 ? 1.0 / Ikf : 0;
  Ikr = Ikr > 0 ? 1.0 / Ikr : 0;
  Vaf = Vaf > 0 ? 1.0 / Vaf : 0;
  Var = Var > 0 ? 1.0 / Var : 0;

  T = kelvin (T);
  Ut = T * kBoverQ;
  Ube = real (getV (NODE_B) - getV (NODE_E)) * pol;
  Ubc = real (getV (NODE_B) - getV (NODE_C)) * pol;

  // critical voltage necessary for bad start values
  UbeCrit = pnCriticalVoltage (Is, Nf * Ut);
  UbcCrit = pnCriticalVoltage (Is, Nr * Ut);
  UbePrev = Ube = pnVoltage (Ube, UbePrev, Ut * Nf, UbeCrit);
  UbcPrev = Ubc = pnVoltage (Ubc, UbcPrev, Ut * Nr, UbcCrit);

  Uce = Ube - Ubc;

  // base-emitter diodes
  gtiny = Ube < - 10 * Ut * Nf ? (Is + Ise) : 0;
  If = pnCurrent (Ube, Is, Ut * Nf);
  Ibei = If / Bf;
  gif = pnConductance (Ube, Is, Ut * Nf);
  gbei = gif / Bf;
  Iben = pnCurrent (Ube, Ise, Ut * Ne);
  gben = pnConductance (Ube, Ise, Ut * Ne);
  Ibe = Ibei + Iben + gtiny * Ube;
  gbe = gbei + gben + gtiny;

  // base-collector diodes
  gtiny = Ubc < - 10 * Ut * Nr ? (Is + Isc) : 0;
  Ir = pnCurrent (Ubc, Is, Ut * Nr);
  Ibci = Ir / Br;
  gir = pnConductance (Ubc, Is, Ut * Nr);
  gbci = gir / Br;
  Ibcn = pnCurrent (Ubc, Isc, Ut * Nc);
  gbcn = pnConductance (Ubc, Isc, Ut * Nc);
  Ibc = Ibci + Ibcn + gtiny * Ubc;
  gbc = gbci + gbcn + gtiny;

  // compute base charge quantities
  Q1 = 1 / (1 - Ubc * Vaf - Ube * Var);
  Q2 = If * Ikf + Ir * Ikr;
  nr_double_t SArg = 1.0 + 4.0 * Q2;
  nr_double_t Sqrt = SArg > 0 ? sqrt (SArg) : 1;
  Qb = Q1 * (1 + Sqrt) / 2;
  dQbdUbe = Q1 * (Qb * Var + gif * Ikf / Sqrt);
  dQbdUbc = Q1 * (Qb * Vaf + gir * Ikr / Sqrt);

  // compute transfer current
  It = (If - Ir) / Qb;

  // compute forward and backward transconductance
  gitf = (gif - If * dQbdUbe / Qb) / Qb;
  gitr = (gir - Ir * dQbdUbc / Qb) / Qb;

  // compute old SPICE values
  go = (gir + It * dQbdUbc) / Qb;
  gm = (gif - It * dQbdUbe) / Qb - go;
  setOperatingPoint ("gm", gm);
  setOperatingPoint ("go", go);

  // calculate current-dependent base resistance
  if (Rbm != 0.0) {
    if (Irb != 0.0) {
      nr_double_t a, b, z;
      a = (Ibci + Ibcn + Ibei + Iben) / Irb;
      a = MAX (a, 1e-12); // enforce positive values
      z = (sqrt (1 + 144 / sqr (M_PI) * a) - 1) / 24 * sqr (M_PI) / sqrt (a);
      b = tan (z);
      Rbb = Rbm + 3 * (Rb - Rbm) * (b - z) / z / sqr (b);
    }
    else {
      Rbb = Rbm + (Rb - Rbm) / Qb;
    }
    rb->setProperty ("R", Rbb);
    rb->calcDC ();
  }

  // compute autonomic current sources
  IeqB = Ibe - gbe * Ube;
  IeqC = Ibc - gbc * Ubc;
#if NEWSGP
  IeqE = It - gitf * Ube + gitr * Ubc;
#else
  IeqE = It - gm * Ube - go * Uce;
#endif
  IeqS = 0;
  setI (NODE_B, (-IeqB - IeqC) * pol);
  setI (NODE_C, (+IeqC - IeqE - IeqS) * pol);
  setI (NODE_E, (+IeqB + IeqE) * pol);
  setI (NODE_S, (+IeqS) * pol);

  // apply admittance matrix elements
#if NEWSGP
  setY (NODE_B, NODE_B, gbc + gbe);
  setY (NODE_B, NODE_C, -gbc);
  setY (NODE_B, NODE_E, -gbe);
  setY (NODE_B, NODE_S, 0);
  setY (NODE_C, NODE_B, -gbc + gitf - gitr);
  setY (NODE_C, NODE_C, gbc + gitr);
  setY (NODE_C, NODE_E, -gitf);
  setY (NODE_C, NODE_S, 0);
  setY (NODE_E, NODE_B, -gbe - gitf + gitr);
  setY (NODE_E, NODE_C, -gitr);
  setY (NODE_E, NODE_E, gbe + gitf);
  setY (NODE_E, NODE_S, 0);
  setY (NODE_S, NODE_B, 0);
  setY (NODE_S, NODE_C, 0);
  setY (NODE_S, NODE_E, 0);
  setY (NODE_S, NODE_S, 0);
#else
  setY (NODE_B, NODE_B, gbc + gbe);
  setY (NODE_B, NODE_C, -gbc);
  setY (NODE_B, NODE_E, -gbe);
  setY (NODE_B, NODE_S, 0);
  setY (NODE_C, NODE_B, -gbc + gm);
  setY (NODE_C, NODE_C, go + gbc);
  setY (NODE_C, NODE_E, -go - gm);
  setY (NODE_C, NODE_S, 0);
  setY (NODE_E, NODE_B, -gbe - gm);
  setY (NODE_E, NODE_C, -go);
  setY (NODE_E, NODE_E, gbe + go + gm);
  setY (NODE_E, NODE_S, 0);
  setY (NODE_S, NODE_B, 0);
  setY (NODE_S, NODE_C, 0);
  setY (NODE_S, NODE_E, 0);
  setY (NODE_S, NODE_S, 0);
#endif
}

void bjt::calcOperatingPoints (void) {

  // fetch device model parameters
  nr_double_t Cje0 = getPropertyDouble ("Cje");
  nr_double_t Vje  = getPropertyDouble ("Vje");
  nr_double_t Mje  = getPropertyDouble ("Mje");
  nr_double_t Cjc0 = getPropertyDouble ("Cjc");
  nr_double_t Vjc  = getPropertyDouble ("Vjc");
  nr_double_t Mjc  = getPropertyDouble ("Mjc");
  nr_double_t Xcjc = getPropertyDouble ("Xcjc");
  nr_double_t Cjs0 = getPropertyDouble ("Cjs");
  nr_double_t Vjs  = getPropertyDouble ("Vjs");
  nr_double_t Mjs  = getPropertyDouble ("Mjs");
  nr_double_t Fc   = getPropertyDouble ("Fc");
  nr_double_t Vtf  = getPropertyDouble ("Vtf");
  nr_double_t Tf   = getPropertyDouble ("Tf");
  nr_double_t Xtf  = getPropertyDouble ("Xtf");
  nr_double_t Itf  = getPropertyDouble ("Itf");
  nr_double_t Tr   = getPropertyDouble ("Tr");

  nr_double_t Cbe, Ube, Ubc, Cbci, Cbcx, Ucs, Cbc, Ccs;

  // interpret zero as infinity for that model parameter
  Vtf = Vtf > 0 ? 1.0 / Vtf : 0;

  Ube = real (getV (NODE_B) - getV (NODE_E)) * pol;
  Ubc = real (getV (NODE_B) - getV (NODE_C)) * pol;
  Ucs = real (getV (NODE_S) - getV (NODE_C)) * pol;

  // depletion capacitance of base-emitter diode
  Cbe = pnCapacitance (Ube, Cje0, Vje, Mje, Fc);
  Qbe = pnCharge (Ube, Cje0, Vje, Mje, Fc);

  // diffusion capacitance of base-emitter diode
  nr_double_t e, Tff, dTffdUbe;
  e = 2 * exp (Ubc * Vtf);
  Tff = Tf * (1 + Xtf * sqr (If / (If + Itf)) * e);
  dTffdUbe = Tf * Xtf * 2 * gif * If * Itf / cubic (If + Itf) * e;
  Cbe += (If * dTffdUbe + Tff * (gif - If / Qb * dQbdUbe)) / Qb;
  Qbe += If * Tff / Qb;

  // depletion and diffusion capacitance of base-collector diode
  nr_double_t Qbc;
  Cbc = pnCapacitance (Ubc, Cjc0, Vjc, Mjc, Fc);
  Cbci = Cbc * Xcjc + Tr * gir;
  Qbc = pnCharge (Ubc, Cjc0, Vjc, Mjc, Fc);
  Qbci = Xcjc * Qbc + Tr * Ir;
  Cbcx = Cbc * (1 - Xcjc);
  Qbcx = Qbc * (1 - Xcjc);

  // depletion capacitance of collector-substrate diode
  Ccs = pnCapacitance (Ucs, Cjs0, Vjs, Mjs);
  Qcs = pnCharge (Ucs, Cjs0, Vjs, Mjs);

  // finally save the operating points
  setOperatingPoint ("Cbe", Cbe);
  setOperatingPoint ("Cbci", Cbci);
  setOperatingPoint ("Cbcx", Cbcx);
  setOperatingPoint ("Ccs", Ccs);
  setOperatingPoint ("gmf", gitf);
  setOperatingPoint ("gmr", gitr);
  setOperatingPoint ("gmu", gbci + gbcn);
  setOperatingPoint ("gpi", gbei + gben);
  setOperatingPoint ("Vbe", Ube);
  setOperatingPoint ("Vbc", Ubc);
  setOperatingPoint ("Vce", Ube - Ubc);
  setOperatingPoint ("Vcs", Ucs);
  setOperatingPoint ("Rbb", Rbb);
  setOperatingPoint ("Ibe", Ibe);
  setOperatingPoint ("Ice", It);
}

void bjt::initSP (void) {
  processCbcx ();
}

void bjt::processCbcx (void) {
  nr_double_t Xcjc = getPropertyDouble ("Xcjc");
  nr_double_t Rbm  = getPropertyDouble ("Rbm");
  nr_double_t Cjc0 = getPropertyDouble ("Cjc");

  /* if necessary then insert external capacitance between internal
     collector node and external base node */
  if (Rbm != 0.0 && Cjc0 != 0.0 && Xcjc != 1.0) {
    if (!deviceEnabled (cbcx)) {
      cbcx = splitCapacitance (this, cbcx, getNet (), "Cbcx", rb->getNode (1),
			       getNode (NODE_C));
    }
    cbcx->setProperty ("C", getOperatingPoint ("Cbcx"));
  }
  else {
    disableCapacitance (this, cbcx, getNet ());
  }
}

void bjt::initAC (void) {
  processCbcx ();
  clearI ();
}

void bjt::calcAC (nr_double_t frequency) {
  setMatrixY (calcMatrixY (frequency));
}

#define qbeState 0 // base-emitter charge state
#define cbeState 1 // base-emitter current state
#define qbcState 2 // base-collector charge state
#define cbcState 3 // base-collector current state
#define qcsState 4 // collector-substrate charge state
#define ccsState 5 // collector-substrate current state

#define qbxState 0 // external base-collector charge state
#define cbxState 1 // external base-collector current state

void bjt::initTR (void) {
  setStates (6);
  initDC ();

  // handle external base-collector capacitance appropriately
  processCbcx ();
  if (deviceEnabled (cbcx)) cbcx->setProperty ("Controlled", getName ());
}

void bjt::calcTR (nr_double_t t) {
  calcDC ();
  calcOperatingPoints ();

  nr_double_t Ube = getOperatingPoint ("Vbe");
  nr_double_t Ubc = getOperatingPoint ("Vbc");
  nr_double_t Ucs = getOperatingPoint ("Vcs");
  nr_double_t Cbe = getOperatingPoint ("Cbe");
  nr_double_t Ccs = getOperatingPoint ("Ccs");
  nr_double_t Cbci = getOperatingPoint ("Cbci");
  nr_double_t Cbcx = getOperatingPoint ("Cbcx");

  // handle Rbb and Cbcx appropriately
  if (Rbb != 0.0) {
    rb->setProperty ("R", Rbb);
    rb->calcTR (t);
    if (deviceEnabled (cbcx)) {
      cbcx->clearI ();
      cbcx->clearY ();
      cbcx->transientCapacitance (qbxState, 1, 2, Cbcx, Ubc, Qbcx);
    }
  }

  // TODO: excess phase
  transientCapacitance (qbeState, NODE_B, NODE_E, Cbe, Ube, Qbe);
  transientCapacitance (qbcState, NODE_B, NODE_C, Cbci, Ubc, Qbci);
  transientCapacitance (qcsState, NODE_S, NODE_C, Ccs, Ucs, Qcs);
}
