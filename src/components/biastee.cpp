/*
 * biastee.cpp - bias T class implementation
 *
 * Copyright (C) 2003, 2004, 2005, 2008 Stefan Jahn <stefan@lkcc.org>
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
 * $Id: biastee.cpp,v 1.11 2008/10/05 17:52:11 ela Exp $
 *
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "component.h"
#include "biastee.h"

biastee::biastee () : circuit (3) {
  type = CIR_BIASTEE;
}

void biastee::initSP (void) {
  allocMatrixS ();
  setS (NODE_1, NODE_1, 0.0);
  setS (NODE_1, NODE_2, 1.0);
  setS (NODE_1, NODE_3, 0.0);
  setS (NODE_2, NODE_1, 1.0);
  setS (NODE_2, NODE_2, 0.0);
  setS (NODE_2, NODE_3, 0.0);
  setS (NODE_3, NODE_1, 0.0);
  setS (NODE_3, NODE_2, 0.0);
  setS (NODE_3, NODE_3, 1.0);
}

void biastee::initDC (void) {
  setISource (false);
  setVoltageSources (1);
  allocMatrixMNA ();
  clearB ();
  clearC ();
  voltageSource (VSRC_1, NODE_2, NODE_3);
}

void biastee::initAC (void) {
  setISource (false);
  setVoltageSources (1);
  allocMatrixMNA ();
  clearB ();
  clearC ();
  voltageSource (VSRC_1, NODE_2, NODE_1);
}

#define fState 0 // flux state
#define vState 1 // voltage state
#define qState 2 // charge state
#define cState 3 // current state

void biastee::initTR (void) {
  initDC ();
  setStates (4);
  setISource (true);
}

void biastee::calcTR (nr_double_t) {
  nr_double_t l = getPropertyDouble ("L");
  nr_double_t c = getPropertyDouble ("C");
  nr_double_t g, r, v;
  nr_double_t i = real (getJ (VSRC_1));

  setState (fState, i * l);
  integrate (fState, l, r, v);
  setD (VSRC_1, VSRC_1, -r);
  setE (VSRC_1, v);

  v = real (getV (NODE_1) - getV (NODE_2));
  setState (qState, c * v);
  integrate (qState, c, g, i);
  setY (NODE_1, NODE_1, +g); setY (NODE_2, NODE_2, +g);
  setY (NODE_1, NODE_2, -g); setY (NODE_2, NODE_1, -g);
  setI (NODE_1 , -i);
  setI (NODE_2 , +i);
}

// properties
struct define_t biastee::cirdef =
  { "BiasT", 3, PROP_COMPONENT, PROP_NO_SUBSTRATE, PROP_LINEAR,
    { PROP_NO_PROP },
    { { "C", PROP_REAL, { 1e-6, PROP_NO_STR }, PROP_POS_RANGE },
      { "L", PROP_REAL, { 1e-6, PROP_NO_STR }, PROP_POS_RANGE },
      PROP_NO_PROP }
  };
