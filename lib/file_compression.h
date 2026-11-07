/* Reimplementation of biomcmc-lib file_compression (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CURUPIXA_FILE_COMPRESSION_H
#define CURUPIXA_FILE_COMPRESSION_H
#include <stddef.h>
#include "status.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct crpx_file crpx_file;
typedef enum { CRPX_FILE_AUTO, CRPX_FILE_RAW, CRPX_FILE_GZIP, CRPX_FILE_BZIP2, CRPX_FILE_XZ } crpx_file_format;
/* mode is "r" or "w". Read AUTO detects magic bytes, write AUTO uses suffix.
 * Unavailable compression returns UNSUPPORTED before creating output. Each handle
 * belongs to one worker. Compressed readers process concatenated members and reject
 * truncated streams and trailing garbage when read to EOF. XZ permits standard stream
 * padding. Explicit compression formats require valid headers. No shell commands,
 * process exit or global error state. */
crpx_status crpx_file_open(const char *path, const char *mode, crpx_file_format format, crpx_file **out);
/* Close always destroys the handle, also on I/O failure. */
crpx_status crpx_file_close(crpx_file *file);
crpx_status crpx_file_read(crpx_file *file, void *buffer, size_t capacity, size_t *read_size);
crpx_status crpx_file_write(crpx_file *file, const void *buffer, size_t length);
/* Grow caller-owned *line (initially NULL, *capacity=0); includes newline and a
 * terminating NUL. NOT_FOUND denotes clean EOF, distinct from IO failure. The
 * caller frees *line. A read/allocation error may consume part of the input. */
crpx_status crpx_file_getline(crpx_file *file, char **line, size_t *capacity, size_t *length);
crpx_file_format crpx_file_get_format(const crpx_file *file);
int crpx_file_format_available(crpx_file_format format);

#ifdef __cplusplus
}
#endif
#endif
