#ifndef PRO_H
#define PRO_H

#include "ins.h"

// ins process
int ins_process(string filename, string filename2);
// read results
int read_resultfile(string filename, InsState *ins_state);
// compare result
void compare_result(InsState *test, InsState *real, InsState *res, double *rms[9]);
// compute residual
void compute_residual(InsState *res, int endepoch, string outfile);
#endif