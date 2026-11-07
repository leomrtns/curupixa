/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_STATUS_H
#define CURUPIXA_STATUS_H

#ifdef __cplusplus
extern "C" {
#endif

/* New APIs return a status and leave ownership with the caller. No process-global
 * last-error variable, callbacks, logging, MPI calls or process termination. */
typedef enum crpx_status {
  CRPX_OK = 0,
  CRPX_ERROR_INVALID_ARGUMENT,
  CRPX_ERROR_ALLOCATION,
  CRPX_ERROR_OVERFLOW,
  CRPX_ERROR_OUT_OF_RANGE,
  CRPX_ERROR_PARSE,
  CRPX_ERROR_TOPOLOGY,
  CRPX_ERROR_NOT_FOUND,
  CRPX_ERROR_IO,
  CRPX_ERROR_UNSUPPORTED,
  CRPX_ERROR_CONVERGENCE
} crpx_status;

const char *crpx_status_string(crpx_status status);

#ifdef __cplusplus
}
#endif
#endif
