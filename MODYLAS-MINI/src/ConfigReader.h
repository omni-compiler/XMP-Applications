/*
 * Copyright (C) 2014 RIKEN AICS
 */
#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#ifdef __cplusplus
extern "C" {
#endif

void ConfigReader_parse(const char *file);

int ConfigReader_get_int(const char *key);

double ConfigReader_get_double(const char *key);

const char *ConfigReader_get_str(const char *key);

int ConfigReader_get_bool(const char *key);

const int *ConfigReader_get_int3(const char *key);

const double *ConfigReader_get_double3(const char *key);

#ifdef __cplusplus
}
#endif

#endif /* CONFIG_READER_H */
