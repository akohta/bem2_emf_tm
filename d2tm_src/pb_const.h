/*
 * d2te_const.h
 *
 *  Created on: Aug 2, 2018
 *      Author: ohta
 */

#ifndef D2TE_CONST_H_
#define D2TE_CONST_H_

// epsilon
#define MEPS   1.0e-15 // machine epsilon
#define IEPS   1.0e-8  // epsilon for Gauss-lendre, Gauss-kronrod quadrature
#define DEPS   1.0e-12 // Double Exponential epsilon
#define CDFV   1.0e2   // cut off parameter for derivative boundary integral equation

// Gauss-Legendre setting
#define GLH    32    // for high precision coefficient

// finite difference settings for boundary value
#define FDH   1.0e-3 // finite difference H
#define FDO   2      // finite difference order, available 2,4,6
#define SFDM  0      // 0:derivative boundary integral equation, 1: finite difference method

// Gmsh data
#define MSHVER   2.2 // version number, 2.2 later
#define MSHASCI  0   // file type, 0:ascii
#define MSHPREC  8   // data size(precision), 8:double
#define ELEMTYPE 8   // element type 8:3-node second order line
#define OPENDID  99  // id for open region 

// 3 point Gauss-Legender and Kronrod extension node and weight
#define A_GL   0.774596669241483377036
#define B_GK   0.434243749346802558002
#define G_GK   0.960491268708020283424
#define W0_GL  0.888888888888888888889
#define WA_GL  0.555555555555555555556
#define W0_GK  0.450916538658474142345
#define WA_GK  0.268488089868333440729
#define WB_GK  0.401397414775962222905
#define WG_GK  0.104656226026467265194
#define i_A_GL 1.290994448735805628393 // 1.0/A_GL

#endif /* D2TE_CONST_H_ */
