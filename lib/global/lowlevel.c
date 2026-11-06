/* 
 * This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file lowlevel.c 
 *  \brief Lowest level basic functions, that should be available to all other modules. 
 */

#include "lowlevel.h"
#include "runtime_internal.h"

void *
crpx_malloc_with_errmsg (const char *file, int line, crpx_global_t context, size_t size)
{
  if (!size) return NULL;
  void *value = malloc(size);
  if (!value) crpx_logger_message(CRPX_LOGLEVEL_ERROR, file, line, context,
                                 "Failed to allocate %zu bytes", size);
  return value;
}

void *
crpx_calloc_with_errmsg (const char *file, int line, crpx_global_t context, size_t count, size_t size)
{
  if (!count || !size) return NULL;
  if (count > SIZE_MAX / size) { errno = ENOMEM; return NULL; }
  void *value = calloc(count, size);
  if (!value) crpx_logger_message(CRPX_LOGLEVEL_ERROR, file, line, context,
                                 "Failed to allocate %zu members of %zu bytes", count, size);
  return value;
}

void *
crpx_realloc_with_errmsg (const char *file, int line, crpx_global_t context, void *pointer, size_t size)
{
  if (!size) { free(pointer); return NULL; }
  void *value = realloc(pointer, size);
  if (!value) crpx_logger_message(CRPX_LOGLEVEL_ERROR, file, line, context,
                                 "Failed to reallocate %zu bytes", size);
  return value;
}

void *
crpx_reallocarray_with_errmsg (const char *file, int line, crpx_global_t context, void *pointer,
                             size_t count, size_t size)
{
  if (size && count > SIZE_MAX / size) {
    errno = ENOMEM;
    crpx_logger_message(CRPX_LOGLEVEL_ERROR, file, line, context, "Allocation size overflow");
    return NULL;
  }
  return crpx_realloc_with_errmsg(file, line, context, pointer, count * size);
}

void
crpx_free_with_errmsg (const char *file, int line, crpx_global_t context, void *pointer)
{
  (void)file; (void)line; (void)context;
  free(pointer);
}

void
crpx_fprintf_colour (FILE *stream, int regular, int colour, const char *message, const char *normaltext, ...)
{
  va_list arguments;
  if (!stream || !message || !normaltext) return;
  if (regular < 0 || regular > 6) regular = 0;
  if (colour < 0 || colour > 7) colour = 1;
  fprintf(stream, "\033[%d;%dm%s\033[0m", regular == 1 ? 1 : 0, 30 + colour, message);
  va_start(arguments, normaltext);
  vfprintf(stream, normaltext, arguments);
  va_end(arguments);
}

static const char *const level_names[] = { "FATAL", "ERROR", "WARNING", "INFO", "VERBOSE", "DEBUG" };

void
crpx_logger_message (uint8_t level, const char *file, int line, crpx_global_t context, const char *format, ...)
{
  if (!context || !format || level > CRPX_LOGLEVEL_DEBUG) return;
  crpx_runtime_lock(&context->runtime_state->logger_lock);
  if (level == CRPX_LOGLEVEL_FATAL) context->error = 2;
  else if (level == CRPX_LOGLEVEL_ERROR && !context->error) context->error = 1;
  FILE *outputs[2] = { stderr, context->logfile };
  uint8_t levels[2] = { context->loglevel_stderr, context->loglevel_file };
  for (size_t i = 0; i < 2; i++) {
    if (!outputs[i] || level > levels[i]) continue;
    va_list arguments;
    fprintf(outputs[i], "[%s] ", level_names[level]);
    va_start(arguments, format);
    vfprintf(outputs[i], format, arguments);
    va_end(arguments);
    if (level <= CRPX_LOGLEVEL_ERROR || level == CRPX_LOGLEVEL_DEBUG)
      fprintf(outputs[i], " [file %s line %d]", file ? file : "?", line);
    fputc('\n', outputs[i]);
    fflush(outputs[i]);
  }
  crpx_runtime_unlock(&context->runtime_state->logger_lock);
}

void
crpx_logger_set_level (crpx_global_t context, uint8_t level)
{
  if (!context) return;
  if (level > CRPX_LOGLEVEL_DEBUG) level = CRPX_LOGLEVEL_DEBUG;
  crpx_runtime_lock(&context->runtime_state->logger_lock);
  context->loglevel_stderr = level;
  crpx_runtime_unlock(&context->runtime_state->logger_lock);
}

void
crpx_logger_set_file (crpx_global_t context, const char *filename, const char *level_string)
{
  if (!context || !filename) return;
  FILE *replacement = fopen(filename, "a");
  if (!replacement) {
    crpx_logger_error(context, "Could not open log file %s", filename);
    return;
  }
  crpx_runtime_lock(&context->runtime_state->logger_lock);
  FILE *previous = context->logfile;
  context->logfile = replacement;
  context->loglevel_file = crpx_get_logger_level_number(level_string, NULL);
  if (previous) fclose(previous);
  crpx_runtime_unlock(&context->runtime_state->logger_lock);
}

uint8_t
crpx_get_logger_level_number (const char *level_string, char *level_output)
{
  uint8_t level = CRPX_LOGLEVEL_ERROR;
  if (level_string) switch (level_string[0]) {
    case 'f': case 'F': level = CRPX_LOGLEVEL_FATAL; break;
    case 'w': case 'W': level = CRPX_LOGLEVEL_WARN; break;
    case 'i': case 'I': level = CRPX_LOGLEVEL_INFO; break;
    case 'v': case 'V': level = CRPX_LOGLEVEL_VERBOSE; break;
    case 'd': case 'D': level = CRPX_LOGLEVEL_DEBUG; break;
    default: break;
  }
  if (level_output) strcpy(level_output, level_names[level]);
  return level;
}
