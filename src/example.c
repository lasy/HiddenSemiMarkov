#include <R.h>
#include <Rinternals.h>

void add_(double* x, double* y, double* out) {
  out[0] = x[0] + y[0];
}
