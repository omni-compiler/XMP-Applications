#include <stdio.h>
#include <math.h>
#include "perf.h"

#define FUNC1 0
#define FUNC2 1

double func1(int n);
double func2(int n);

int main(int argc, char **argv)
{
  int n = 1000000;
  double a, b;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

  perf_start(FUNC1);
  a =  func1(n);
  perf_stop(FUNC1);

  perf_start(FUNC2);
  b =  func2(n);
  perf_stop(FUNC2);

  printf("a+b = %f\n", a + b);

#ifdef USE_MPI
  perf_print_time_mpi(FUNC1, "func1: ");
  perf_print_time_mpi_full(FUNC2, "func2: ");
#else
  perf_print_time(FUNC1, "func1: ");
  perf_print_time(FUNC2, "func2: ");
#endif

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}


double func1(int n)
{
  int i;
  double x;
  x = 0.0;
  for (i = 0; i < n; i++) {
    x += exp(-pow(sqrt(fabs(sin(i*3.14))), 1.234));
  }
  return x;
}


double func2(int n)
{
  return func1(n);
}

