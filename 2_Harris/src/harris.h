#include "ap_int.h"
#define CPU

typedef ap_uint<8> DTYPE;

#define G_LEN 16
#define GB_LEN 12
#define R_LEN 16

void harris(DTYPE* imgSrc, DTYPE* imgDst, double alpha);
void harris_cpu(DTYPE* imgSrc, DTYPE* imgDst, double alpha);
