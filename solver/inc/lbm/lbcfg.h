/*! \file lbcfg.h
 *
 * Lattice Boltzmann Model configuration.
 *
 * \author Yu Zhou
 */

/*!
 * \brief SCSP reference-density equilibrium and velocity when zero.
 * ColorFluid uses local total density for its velocity reconstruction.
 * The historical macro spelling and default are preserved.
 */
#define ISCOMPRESSABLE 0
