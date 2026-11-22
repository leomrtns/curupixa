/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file quasi_sequence.h
    \brief Context-free, seekable low-discrepancy sequences with deterministic reset.
    Bratley and Fox (1988), doi:10.1145/42288.214372. */
#ifndef CURUPIXA_QUASI_SEQUENCE_H
#define CURUPIXA_QUASI_SEQUENCE_H
#include "status.h"
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { CRPX_QUASI_HALTON, CRPX_QUASI_REVERSE_HALTON, CRPX_QUASI_SOBOL } crpx_quasi_kind;
typedef struct crpx_quasi_sequence crpx_quasi_sequence;
/* Halton supports 1..2048 distinct prime bases; Sobol supports 1..8 dimensions.
   No random context or random initialisation is used. Unsupported dimensions fail
   rather than reuse bases. The first next() returns point 1 (point 0 is all zero).
   Values are doubles in [0,1). Index is limited to 2^53-1; exhaustion returns OVERFLOW. */
crpx_status crpx_quasi_sequence_new (crpx_quasi_kind kind, size_t dimensions, crpx_quasi_sequence **out);
void crpx_quasi_sequence_free (crpx_quasi_sequence *sequence);
crpx_status crpx_quasi_sequence_seek (crpx_quasi_sequence *sequence, uint64_t next_index);
crpx_status crpx_quasi_sequence_next (crpx_quasi_sequence *sequence, double *point, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
