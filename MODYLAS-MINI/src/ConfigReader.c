/*
 * Copyright (C) 2014 RIKEN AICS
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "ConfigReader.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

typedef enum {
  CONFIG_INT,
  CONFIG_DOUBLE,
  CONFIG_STR,
  CONFIG_BOOL,
  CONFIG_INT3,
  CONFIG_DOUBLE3,
} ConfigType;

typedef struct {
  const char *key;
  ConfigType type;
  const char *val;
  union {
    int int_val;
    double double_val;
    const char *str_val;
    int bool_val;
    int int_array[3];
    double double_array[3];
  } data;
} ConfigParam;

static ConfigParam params[] = {
#include "params.h"
};

static int n_params = sizeof(params)/sizeof(params[0]);

static int n_errors = 0;

#define BUF_SIZE 1024
#define MAX_KEY_LEN 512
#define MAX_VAL_LEN 512

char section[MAX_KEY_LEN] = "";

#ifdef USE_MPI
static MPI_Comm comm = MPI_COMM_WORLD;
static int myrank;
#endif

/*
static void print_message(const char *fmt, ...)
{
#ifdef USE_MPI
  if (myrank != 0) return;
#endif
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}
*/

static void abort_prog(int code)
{
#ifdef USE_MPI
  MPI_Abort(comm, code);
#else
  exit(code);
#endif
}

static char *trim_left(char *s)
{
  for (; *s != '\0' && isspace(*s); s++);
  return s;
}

static char *trim_right(char *s)
{
  int i;
  for (i = strlen(s)-1; i >= 0 && isspace(s[i]); i--) {
    s[i] = '\0';
  }
  return s;
}

#define STRCAT(s1, s2)  strncat((s1), (s2), sizeof(s1)-strlen(s1)-1)
#define STRCPY(s1, s2)  (strncpy((s1), (s2), sizeof(s1)-1), s1[sizeof(s1)-1] = '\0')

static int check_section(char *s)
{
  if (sscanf(s, "[%[^]]", section) == 1) {
    STRCAT(section, ".");
//  printf("[%s]\n", section);
    return 1;
  } else {
    return 0;
  }
}

static int scan_line(FILE *f, char **key, char **val)
{
  static char buf[BUF_SIZE];
  static char k[MAX_KEY_LEN];
  static char v[MAX_VAL_LEN];
  char *p, *p0, *p1;

  *key = k;
  *val = v;

  while (fgets(buf, sizeof(buf), f) != NULL) {
    if ((p = strchr(buf, '#')) != NULL) *p = '\0';
    p0 = trim_left(buf);
    if (*p0 == '\0') continue;
    if (check_section(p0)) continue;
    p1 = strchr(buf, '=');
    if (p1 != NULL) {
      *p1 = '\0';
      p1++;
      p0 = trim_right(p0);
      p1 = trim_right(trim_left(p1));
      STRCPY(k, section);
      STRCAT(k, p0);
      STRCPY(v, p1);
//    printf("key = [%s], val = [%s]\n", k, v);
      return 1;
    } else {
//      print_message("*** illegal line: %s", buf);
      printf("*** illegal line: %s", buf);
      n_errors++;
    }
  }

  return 0;
}

static void parse_param(ConfigParam *param)
{
  switch (param->type) {
  case CONFIG_INT:
    if (sscanf(param->val, "%d", &param->data.int_val) != 1) {
//      print_message("*** cannot convert to int: key = %s, val = %s\n",
      printf("*** cannot convert to int: key = %s, val = %s\n",
                    param->key, param->val);
      n_errors++;
    }
    break;
  case CONFIG_DOUBLE:
    if (sscanf(param->val, "%lf", &param->data.double_val) != 1) {
//      print_message("*** cannot convert to double: key = %s, val = %s\n",
      printf("*** cannot convert to double: key = %s, val = %s\n",
                    param->key, param->val);
      n_errors++;
    }
    break;
  case CONFIG_STR:
    param->data.str_val = param->val;
    break;
  case CONFIG_BOOL:
    if (strcasecmp(param->val, "yes") == 0 ||
        strcasecmp(param->val, "on") == 0 ||
        strcasecmp(param->val, "true") == 0) {
      param->data.bool_val = 1;
    }
    else if (strcasecmp(param->val, "no") == 0 ||
             strcasecmp(param->val, "off") == 0 ||
             strcasecmp(param->val, "false") == 0) {
      param->data.bool_val = 0;
    }
    else {
//      print_message("*** cannot convert to bool: key = %s, val = %s\n",
      printf("*** cannot convert to bool: key = %s, val = %s\n",
                    param->key, param->val);
      n_errors++;
    }
    break;
  case CONFIG_INT3:
    if (sscanf(param->val, "%d %d %d",
               &param->data.int_array[0],
               &param->data.int_array[1],
               &param->data.int_array[2]) != 3) {
//      print_message("*** cannot convert to int[3]: key = %s, val = %s\n",
      printf("*** cannot convert to int[3]: key = %s, val = %s\n",
                    param->key, param->val);
      n_errors++;
    }
    break;
  case CONFIG_DOUBLE3:
    if (sscanf(param->val, "%lf %lf %lf",
               &param->data.double_array[0],
               &param->data.double_array[1],
               &param->data.double_array[2]) != 3) {
//      print_message("*** cannot convert to double[3]: key = %s, val = %s\n",
      printf("*** cannot convert to double[3]: key = %s, val = %s\n",
                    param->key, param->val);
      n_errors++;
    }
    break;
  }
}

static ConfigParam *search_param(const char *key)
{
  int i;
  for (i = 0; i < n_params; i++) {
    if (strcmp(key, params[i].key) == 0) return &params[i];
  }
  return NULL;
}

static ConfigParam *search_param_or_abort(const char *key)
{
  ConfigParam *param = search_param(key);
  if (!param) {
//    print_message("*** invalid key: %s\n", key);
    printf("*** invalid key: %s\n", key);
    abort_prog(1);
  }
  return param;
}

static void store_key_val(const char *key, const char *val)
{
  ConfigParam *param = search_param(key);
  if (param) {
    param->val = strdup(val);
  } else {
//    print_message("*** unknown key: %s\n", key);
    printf("*** unknown key: %s\n", key);
    n_errors++;
  }
}

#ifdef USE_MPI
static char *dump_vals(int *size)
{
  int i;
  int n = 0;
  char *vals, *d;
  const char *s;
  for (i = 0; i < n_params; i++) {
    assert(params[i].val);
    n += strlen(params[i].val) + 1;
  }
  vals = malloc(n);
//  printf("*** [%03d] vals= %s, address= %p, n= %d ***\n", myrank, vals, vals, n);
//  assert( n==90 );
  d = vals;
  for (i = 0; i < n_params; i++) {
    for (s = params[i].val; *s != '\0'; s++, d++) {
      *d = *s;
    }
    *d++ = '\0';
  }
  assert(d == &vals[n]);
  *size = n;

  return vals;
}

static void store_vals(const char* vals)
{
  int i;
  const char *s = vals;
  for (i = 0; i < n_params; i++) {
    params[i].val = strdup(s);
    s += strlen(s) + 1;
  }
}
#endif

static void parse_params()
{
  int i;
  for (i = 0; i < n_params; i++) {
    // printf("key = [%s], val = [%s]\n", params[i].key, params[i].val);
    if (params[i].val != NULL) {
      parse_param(&params[i]);
    } else {
//      print_message("*** not specified key: %s\n", params[i].key);
      printf("*** not specified key: %s\n", params[i].key);
      n_errors++;
    }
  }
}

void ConfigReader_parse(const char *file)
{
  FILE *f;
  char *key, *val;
#ifdef USE_MPI
  char *vals;
  int size;
  MPI_Comm_rank(comm, &myrank);
  if (myrank == 0) {
#endif
    if ((f = fopen(file, "r")) == NULL) {
      perror(file);
      exit(1);
    }
    while (scan_line(f, &key, &val)) {
      store_key_val(key, val);
    }
    fclose(f);

    parse_params();

    if (n_errors > 0) {
//      print_message("*** %d error(s) detected in config file %s\n",
      printf("*** %d error(s) detected in config file %s\n",
                    n_errors, file);
      abort_prog(1);
    }

#ifdef USE_MPI
    /* myrank == 0 */
    vals = dump_vals(&size);
//    printf("[%03d] rank = 0, (1)MPI_Bcast: vals= %s, addr= %p\n", myrank, vals, vals);
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    MPI_Bcast(vals, size, MPI_CHAR, 0, comm);
//    printf("[%03d] rank = 0, (2)MPI_Bcast: vals= %s, addr= %p\n", myrank, vals, vals);
    free(vals);
  } else {
    /* myrank != 0 */
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    vals = malloc(size);
//    printf("[%03d] rank /= 0, (1)MPI_Bcast: vals= %s, addr= %p\n", myrank, vals, vals);
    MPI_Bcast(vals, size, MPI_CHAR, 0, comm);
//    printf("[%03d] rank /= 0, (2)MPI_Bcast: vals= %s, addr= %p\n", myrank, vals, vals);
    store_vals(vals);
    free(vals);
    parse_params();
  }
#endif

}
void parse_(const char *file)
{
  ConfigReader_parse(file);
}

int ConfigReader_get_int(const char *key)
{
  const ConfigParam *param = search_param_or_abort(key);
  if (param->type != CONFIG_INT) {
//    print_message("*** key %s has no int type\n", key);
    printf("*** key %s has no int type\n", key);
    abort_prog(1);
  }
  return param->data.int_val;
}
int get_int_(const char *key)
{
  return ConfigReader_get_int(key);
}

double ConfigReader_get_double(const char *key)
{
  const ConfigParam *param = search_param_or_abort(key);
  if (param->type != CONFIG_DOUBLE) {
//    print_message("*** key %s has no double type\n", key);
    printf("*** key %s has no double type\n", key);
    abort_prog(1);
  }
  return param->data.double_val;
}
double get_double_(const char *key)
{
  return ConfigReader_get_double(key);
}

const char *ConfigReader_get_str(const char *key)
{
  const ConfigParam *param = search_param_or_abort(key);
  if (param->type != CONFIG_STR) {
//    print_message("*** key %s has no str type\n", key);
    printf("*** key %s has no str type\n", key);
    abort_prog(1);
  }
  return param->data.str_val;
}

//int ConfigReader_get_bool(const char *key)
int ConfigRead_get_bool(const char *key)
{
  const ConfigParam *param = search_param_or_abort(key);
  if (param->type != CONFIG_BOOL) {
//    print_message("*** key %s has no bool type\n", key);
    printf("*** key %s has no bool type\n", key);
    abort_prog(1);
  }
  return param->data.bool_val;
}
int get_bool_(const char *key)
{
  return ConfigRead_get_bool(key);
}

/*---
const int *ConfigReader_get_int3(const char *key)
{
  const ConfigParam *param = search_param_or_abort(key);
  if (param->type != CONFIG_INT3) {
    print_message("*** key %s has no int[3] type\n", key);
    abort_prog(1);
  }
  return param->data.int_array;
}

const double *ConfigReader_get_double3(const char *key)
{
  const ConfigParam *param = search_param_or_abort(key);
  if (param->type != CONFIG_DOUBLE3) {
    print_message("*** key %s has no double[3] type\n", key);
    abort_prog(1);
  }
  return param->data.double_array;
}
---*/
