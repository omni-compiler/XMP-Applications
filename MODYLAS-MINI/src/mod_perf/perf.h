#ifndef PERF_H
#define PERF_H

/**
 * value type constans for getter functions.
 */
typedef enum {
  PERF_ROOT,   /**< rank0 */
  PERF_AVE,    /**< average */
  PERF_MIN,    /**< minimum */
  PERF_MAX,    /**< maximum */
  PERF_SD,     /**< standard deviation */
} PERF_ValueType;


#ifdef __cplusplus
extern "C" {
#endif

void perf_start(int id);

void perf_stop(int id);

void perf_add_fp_ops(int id, double ops);

void perf_add_ld_ops(int id, double ops);

void perf_add_st_ops(int id, double ops);

void perf_add_ld_min_ops(int id, double ops);

void perf_add_st_min_ops(int id, double ops);

double perf_get_time(int id, PERF_ValueType type);

double perf_get_flops(int id, PERF_ValueType type);

double perf_get_throughput(int id, PERF_ValueType type);

double perf_get_effective_throughput(int id, PERF_ValueType type);

void perf_print_time(int id, const char *name);

void perf_print_time_mpi(int id, const char *name);

void perf_print_time_mpi_full(int id, const char *name);

#ifdef __cplusplus
}
#endif

#endif /* PERF_H */
