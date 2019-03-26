#ifndef PRO_H
#define PRO_H

#include "ins.h"

// ins process
int ins_process(char* filename, char* filename2, int test_flag);
// read results
int read_resultfile(char* filename, InsState *ins_state);
// compare result
void compare_result(InsState *test, InsState *real, InsState *res, double rms[9]);
// compute residual
void compute_residual(InsState *res, char* outfile, double acc);
// compute result
void compute_result(InsState *test);

#endif