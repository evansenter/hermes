#ifndef VIENNA_DATA_STRUCTURES_H
#define VIENNA_DATA_STRUCTURES_H

#include "energy_par.h"

typedef struct {
  float energy;
  char* structure;
} SOLUTION;

#define   VRNA_GQUAD_MAX_STACK_SIZE     7
#define   VRNA_GQUAD_MIN_STACK_SIZE     2
#define   VRNA_GQUAD_MAX_LINKER_LENGTH  15
#define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
#define   VRNA_GQUAD_MIN_BOX_SIZE       ((4*VRNA_GQUAD_MIN_STACK_SIZE)+(3*VRNA_GQUAD_MIN_LINKER_LENGTH))
#define   VRNA_GQUAD_MAX_BOX_SIZE       ((4*VRNA_GQUAD_MAX_STACK_SIZE)+(3*VRNA_GQUAD_MAX_LINKER_LENGTH))

/**
 *  \brief The data structure that contains the complete model details used throughout the calculations
 *
 */
typedef struct {
  int     dangles;      /**<  \brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
                              \note   Some function do not implement all dangle model but only a subset of
                                      (0,1,2,3). Read the documentaion of the particular recurrences or
                                      energy evaluation function for information about the provided dangle
                                      model.
                        */
  int     special_hp;   /**<  \brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int     noLP;         /**<  \brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int     noGU;         /**<  \brief  Do not allow GU pairs */
  int     noGUclosure;  /**<  \brief  Do not allow loops to be closed by GU pair */
  int     logML;        /**<  \brief  Use logarithmic scaling for multi loops */
  int     circ;         /**<  \brief  Assume molecule to be circular */
  int     gquad;        /**<  \brief  Include G-quadruplexes in structure prediction */
} model_detailsT;

/**
 *  \brief The datastructure that contains temperature scaled energy parameters.
 */
typedef struct {
  int id;
  int stack[NBPAIRS + 1][NBPAIRS + 1];
  int hairpin[31];
  int bulge[MAXLOOP + 1];
  int internal_loop[MAXLOOP + 1];
  int mismatchExt[NBPAIRS + 1][5][5];
  int mismatchI[NBPAIRS + 1][5][5];
  int mismatch1nI[NBPAIRS + 1][5][5];
  int mismatch23I[NBPAIRS + 1][5][5];
  int mismatchH[NBPAIRS + 1][5][5];
  int mismatchM[NBPAIRS + 1][5][5];
  int dangle5[NBPAIRS + 1][5];
  int dangle3[NBPAIRS + 1][5];
  int int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  int int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  int int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  int ninio[5];
  double  lxc;
  int     MLbase;
  int     MLintern[NBPAIRS + 1];
  int     MLclosing;
  int     TerminalAU;
  int     DuplexInit;
  int     Tetraloop_E[200];
  char    Tetraloops[1401];
  int     Triloop_E[40];
  char    Triloops[241];
  int     Hexaloop_E[40];
  char    Hexaloops[1801];
  int     TripleC;
  int     MultipleCA;
  int     MultipleCB;
  int     gquad [VRNA_GQUAD_MAX_STACK_SIZE + 1]
  [3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  
  double  temperature;            /**<  \brief  Temperature used for loop contribution scaling */
  
  model_detailsT model_details;   /**<  \brief  Model details to be used in the recursions */
  
} paramT;

/**
 *  \brief  The datastructure that contains temperature scaled Boltzmann weights of the energy parameters.
 */
typedef struct {
  int     id;
  double  expstack[NBPAIRS + 1][NBPAIRS + 1];
  double  exphairpin[31];
  double  expbulge[MAXLOOP + 1];
  double  expinternal[MAXLOOP + 1];
  double  expmismatchExt[NBPAIRS + 1][5][5];
  double  expmismatchI[NBPAIRS + 1][5][5];
  double  expmismatch23I[NBPAIRS + 1][5][5];
  double  expmismatch1nI[NBPAIRS + 1][5][5];
  double  expmismatchH[NBPAIRS + 1][5][5];
  double  expmismatchM[NBPAIRS + 1][5][5];
  double  expdangle5[NBPAIRS + 1][5];
  double  expdangle3[NBPAIRS + 1][5];
  double  expint11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  double  expint21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  double  expint22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  double  expninio[5][MAXLOOP + 1];
  double  lxc;
  double  expMLbase;
  double  expMLintern[NBPAIRS + 1];
  double  expMLclosing;
  double  expTermAU;
  double  expDuplexInit;
  double  exptetra[40];
  double  exptri[40];
  double  exphex[40];
  char    Tetraloops[1401];
  double  expTriloop[40];
  char    Triloops[241];
  char    Hexaloops[1801];
  double  expTripleC;
  double  expMultipleCA;
  double  expMultipleCB;
  double  expgquad[VRNA_GQUAD_MAX_STACK_SIZE + 1]
  [3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  
  double  kT;
  double  pf_scale;     /**<  \brief    Scaling factor to avoid over-/underflows */
  
  double  temperature;  /**<  \brief    Temperature used for loop contribution scaling */
  double  alpha;        /**<  \brief    Scaling factor for the thermodynamic temperature
                              \details  This allows for temperature scaling in Boltzmann
                                        factors independently from the energy contributions.
                                        The resulting Boltzmann factors are then computed by
                                        \f$ e^{-E/(\alpha \cdot K \cdot T)} \f$
                        */

  model_detailsT model_details; /**<  \brief  Model details to be used in the recursions */
  
}  pf_paramT;

#endif
