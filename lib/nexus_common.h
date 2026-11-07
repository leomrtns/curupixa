/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_NEXUS_COMMON_H
#define CURUPIXA_NEXUS_COMMON_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Remove nested [] comments while respecting single-quoted labels/doubled
 * quotes. Preserve line breaks and replace comments with spaces. Caller frees
 * result. Unmatched brackets/quotes are parse errors. */
crpx_status crpx_nexus_strip_comments(const char *text, char **result, size_t *error_offset);
#ifdef __cplusplus
}
#endif
#endif
