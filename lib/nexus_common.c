/* Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "nexus_common.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

crpx_status crpx_nexus_strip_comments(const char *text, char **result, size_t *error_offset)
{
  size_t depth = 0, length;
  int quoted = 0;
  char *copy;
  if (!result) return CRPX_ERROR_INVALID_ARGUMENT;
  *result = NULL;
  if (!text) return CRPX_ERROR_INVALID_ARGUMENT;
  length = strlen(text);
  if (length == SIZE_MAX) return CRPX_ERROR_OVERFLOW;
  copy = malloc(length + 1);
  if (!copy) return CRPX_ERROR_ALLOCATION;
  for (size_t i = 0; i < length; ++i) {
    char c = text[i];
    if (depth) {
      if (c == '[') ++depth;
      else if (c == ']') --depth;
      copy[i] = c == '\n' || c == '\r' ? c : ' ';
      continue;
    }
    copy[i] = c;
    if (c == '\'') {
      if (quoted && text[i+1] == '\'') { copy[i+1] = text[i+1]; ++i; }
      else quoted = !quoted;
    } else if (!quoted && c == '[') { depth = 1; copy[i] = ' '; }
    else if (!quoted && c == ']') { if (error_offset) *error_offset = i; free(copy); return CRPX_ERROR_PARSE; }
  }
  if (depth || quoted) { if (error_offset) *error_offset = length; free(copy); return CRPX_ERROR_PARSE; }
  copy[length] = '\0';
  *result = copy;
  return CRPX_OK;
}
