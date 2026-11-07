/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "status.h"

const char *crpx_status_string(crpx_status status)
{
  switch (status) {
    case CRPX_OK: return "success";
    case CRPX_ERROR_INVALID_ARGUMENT: return "invalid argument";
    case CRPX_ERROR_ALLOCATION: return "memory allocation failed";
    case CRPX_ERROR_OVERFLOW: return "size or arithmetic overflow";
    case CRPX_ERROR_OUT_OF_RANGE: return "index out of range";
    case CRPX_ERROR_PARSE: return "invalid input syntax";
    case CRPX_ERROR_TOPOLOGY: return "invalid tree topology";
    case CRPX_ERROR_NOT_FOUND: return "item not found";
    case CRPX_ERROR_IO: return "input/output failure";
    case CRPX_ERROR_UNSUPPORTED: return "unsupported operation";
    case CRPX_ERROR_CONVERGENCE: return "numerical iteration did not converge";
  }
  return "unknown status";
}
