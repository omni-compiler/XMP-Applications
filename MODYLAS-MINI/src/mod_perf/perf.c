#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#ifdef USE_MPI
#include "mpi.h"
#endif
#include "perf.h"


typedef struct {
  int id;

  /* timing */
  double start;
  double time;
  long count;
  int started;

  /* operation counters */
  double fp_ops;      /* floating point */
  double ld_ops;      /* load */
  double st_ops;      /* sotre */
  double ld_min_ops;  /* load (effective) */
  double st_min_ops;  /* sotre (effective) */
} Section;


static Section sections[MAX_SECTIONS];


static double get_current_time()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
}


static void error_exit(const char *fmt, ...)
{
  const int error_code = 1;

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD, error_code);
#else
  exit(error_code);
#endif 
}


static double get_value(double val, PERF_ValueType type)
{
  switch (type) {
    case PERF_ROOT:
      return val;
#ifdef USE_MPI
    case PERF_AVE:
    {
      int nprocs;
      double sum = 0.0;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Reduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return sum / nprocs;
    }
    case PERF_MIN:
    {
      double min = 0.0;
      MPI_Reduce(&val, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      return min;
    }
    case PERF_MAX:
    {
      double max = 0.0;
      MPI_Reduce(&val, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      return max;
    }
    case PERF_SD:
    {
      int nprocs;
      double sum2 = 0;
      double val2= val*val;
      double ave = get_value(val, PERF_AVE);
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Reduce(&val2, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      return sqrt(sum2/nprocs - ave*ave);
    }
#endif
    default:
      error_exit("***error: get_value: unknown Value Type<%d>\n", type);
  }
  /* NOTREACHED */
}


void perf_start(int id)
{
  sections[id].started = 1;
  sections[id].start = get_current_time();
}


void perf_stop(int id)
{
  if (!sections[id].started) error_exit("***error: perf_stop: section<%d> not started\n", id);
  sections[id].time += get_current_time() - sections[id].start;
  sections[id].count++;
  sections[id].started = 0;
}


void perf_add_fp_ops(int id, double ops)
{
  sections[id].fp_ops += ops;
}


void perf_add_ld_ops(int id, double ops)
{
  sections[id].ld_ops += ops;
}


void perf_add_st_ops(int id, double ops)
{
  sections[id].st_ops += ops;
}


void perf_add_ld_min_ops(int id, double ops)
{
  sections[id].ld_min_ops += ops;
}


void perf_add_st_min_ops(int id, double ops)
{
  sections[id].st_min_ops += ops;
}


double perf_get_time(int id, PERF_ValueType type)
{
  return get_value(sections[id].time, type);
}


double perf_get_flops(int id, PERF_ValueType type)
{
  double val_local;
  if (sections[id].time > 0.0) {
    val_local = sections[id].fp_ops / sections[id].time 
                * 1.0e-9;  /* in unit of GFLOPS */
  } else {
    val_local = 0.0;
  }
  return get_value(val_local, type);
}


double perf_get_throughput(int id, PERF_ValueType type)
{
  double val_local;
  if (sections[id].time > 0.0) {
    val_local = (sections[id].ld_ops + sections[id].st_ops) / sections[id].time
                * 8 * 1.0e-9;  /* in unit of GB/s */
  } else {
    val_local = 0.0;
  }
  return get_value(val_local, type);
}


double perf_get_effective_throughput(int id, PERF_ValueType type)
{
  double val_local;
  if (sections[id].time > 0.0) {
    val_local = (sections[id].ld_min_ops + sections[id].st_min_ops) / sections[id].time 
                 * 8 * 1.0e-9;  /* in unit of GB/s */
  } else {
    val_local = 0.0;
  }
  return get_value(val_local, type);
}


void perf_print_time(int id, const char *name)
{
  int myrank = 0;
  double time = perf_get_time(id, PERF_ROOT);
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
  if (myrank == 0) printf("%s %f [sec]\n", name, time);
}


#ifdef USE_MPI
void perf_print_time_mpi(int id, const char *name)
{
  int myrank;
  double time_root = perf_get_time(id, PERF_ROOT);
  double time_min  = perf_get_time(id, PERF_MIN);
  double time_max  = perf_get_time(id, PERF_MAX);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    printf("%s (rank0)%f / (min)%f / (max)%f [sec]\n",
           name, time_root, time_min, time_max);
  }
}


void perf_print_time_mpi_full(int id, const char *name)
{
  int myrank;
  double time_root = perf_get_time(id, PERF_ROOT);
  double time_min  = perf_get_time(id, PERF_MIN);
  double time_max  = perf_get_time(id, PERF_MAX);
  double time_ave  = perf_get_time(id, PERF_AVE);
  double time_sd   = perf_get_time(id, PERF_SD);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    printf("%s\n", name);
    printf("    rank0: %f [sec]\n", time_root);
    printf("    min:   %f [sec]\n", time_min);
    printf("    max:   %f [sec]\n", time_max);
    printf("    ave:   %f [sec]\n", time_ave);
    printf("    sd:    %f [sec]\n", time_sd);
  }
}
#endif
