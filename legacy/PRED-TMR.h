/*
********************************************************************************
* PRED-TMR.h                                                   Claude Pasquier *
* Created 10/06/98                                        University of Athens *
*                                                                              *
* Predict transmembrane domains in a protein sequence                          *
*                                                                              *
*                                                                              *
*                                                                              *
********************************************************************************
*/

#include "seqfile.h"
#include "sequence.h"
#include "common.h"

bool main(int argc, char *argv[]);
bool ProcessFile(char* fileName);
bool ComputeIndicators();
bool MarkTM();
void DisplayResult();
bool ParseArg(int argc, char *argv[]);
bool ParseArg();
bool rating(float& Qa, float& Ca, float& SM);
void Generate_CoPreTHi_Page();
