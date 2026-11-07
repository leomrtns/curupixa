/* Reimplementation of biomcmc-lib file_compression (Leonardo de Oliveira Martins).
 * Copyright (C) 2019-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "file_compression.h"
#include "config.h"
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif
#ifdef HAVE_BZIP2
#include <bzlib.h>
#endif
#ifdef HAVE_LZMA
#include <lzma.h>
#endif

struct crpx_file {
  crpx_file_format format;
  bool writing, input_eof, read_finished, member_end, codec_initialised;
  crpx_status failure;
  FILE *raw;
#ifdef HAVE_ZLIB
  z_stream gzip;
#endif
#ifdef HAVE_BZIP2
  bz_stream bzip2;
#endif
#ifdef HAVE_LZMA
  lzma_stream xz;
#endif
  unsigned char input[8192], lines[8192];
  size_t input_position, input_size, line_position, line_size;
};

int crpx_file_format_available(crpx_file_format format)
{
  switch (format) {
    case CRPX_FILE_AUTO: case CRPX_FILE_RAW: return 1;
#ifdef HAVE_ZLIB
    case CRPX_FILE_GZIP: return 1;
#endif
#ifdef HAVE_BZIP2
    case CRPX_FILE_BZIP2: return 1;
#endif
#ifdef HAVE_LZMA
    case CRPX_FILE_XZ: return 1;
#endif
    default: return 0;
  }
}

static int has_suffix(const char *text, const char *suffix)
{
  size_t a = strlen(text), b = strlen(suffix);
  return a >= b && !strcmp(text + a - b, suffix);
}

static crpx_status initialise_codec(crpx_file *file)
{
  switch (file->format) {
    case CRPX_FILE_RAW: return CRPX_OK;
#ifdef HAVE_ZLIB
    case CRPX_FILE_GZIP: {
      int code = file->writing ? deflateInit2(&file->gzip, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 31, 8, Z_DEFAULT_STRATEGY) :
                                 inflateInit2(&file->gzip, 31);
      if (code != Z_OK) return code == Z_MEM_ERROR ? CRPX_ERROR_ALLOCATION : CRPX_ERROR_IO;
      break;
    }
#endif
#ifdef HAVE_BZIP2
    case CRPX_FILE_BZIP2: {
      int code = file->writing ? BZ2_bzCompressInit(&file->bzip2, 9, 0, 30) : BZ2_bzDecompressInit(&file->bzip2, 0, 0);
      if (code != BZ_OK) return code == BZ_MEM_ERROR ? CRPX_ERROR_ALLOCATION : CRPX_ERROR_IO;
      break;
    }
#endif
#ifdef HAVE_LZMA
    case CRPX_FILE_XZ: {
      lzma_ret code = file->writing ? lzma_easy_encoder(&file->xz, 6, LZMA_CHECK_CRC64) :
                                     lzma_stream_decoder(&file->xz, UINT64_MAX, LZMA_CONCATENATED);
      if (code != LZMA_OK) return code == LZMA_MEM_ERROR ? CRPX_ERROR_ALLOCATION : CRPX_ERROR_IO;
      break;
    }
#endif
    default: return CRPX_ERROR_UNSUPPORTED;
  }
  file->codec_initialised = true;
  return CRPX_OK;
}

static void destroy_codec(crpx_file *file)
{
  if (!file->codec_initialised) return;
  switch (file->format) {
#ifdef HAVE_ZLIB
    case CRPX_FILE_GZIP: if (file->writing) deflateEnd(&file->gzip); else inflateEnd(&file->gzip); break;
#endif
#ifdef HAVE_BZIP2
    case CRPX_FILE_BZIP2: if (file->writing) BZ2_bzCompressEnd(&file->bzip2); else BZ2_bzDecompressEnd(&file->bzip2); break;
#endif
#ifdef HAVE_LZMA
    case CRPX_FILE_XZ: lzma_end(&file->xz); break;
#endif
    default: break;
  }
  file->codec_initialised = false;
}

crpx_status crpx_file_open(const char *path, const char *mode, crpx_file_format format, crpx_file **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out = NULL;
  if (!path || !mode || (strcmp(mode, "r") && strcmp(mode, "w")) ||
      format < CRPX_FILE_AUTO || format > CRPX_FILE_XZ) return CRPX_ERROR_INVALID_ARGUMENT;
  bool writing = mode[0] == 'w';
  if (format == CRPX_FILE_AUTO && writing)
    format = has_suffix(path, ".gz") ? CRPX_FILE_GZIP : has_suffix(path, ".bz2") ? CRPX_FILE_BZIP2 :
             has_suffix(path, ".xz") ? CRPX_FILE_XZ : CRPX_FILE_RAW;
  if (format != CRPX_FILE_AUTO && !crpx_file_format_available(format)) return CRPX_ERROR_UNSUPPORTED;
  crpx_file *file = calloc(1, sizeof(*file));
  if (!file) return CRPX_ERROR_ALLOCATION;
  file->writing = writing;
  file->raw = fopen(path, writing ? "wb" : "rb");
  if (!file->raw) { free(file); return CRPX_ERROR_IO; }
  if (format == CRPX_FILE_AUTO) {
    file->input_size = fread(file->input, 1, 6, file->raw);
    if (ferror(file->raw)) { fclose(file->raw); free(file); return CRPX_ERROR_IO; }
    const unsigned char *magic = file->input;
    size_t count = file->input_size;
    format = count >= 2 && magic[0] == 0x1f && magic[1] == 0x8b ? CRPX_FILE_GZIP :
             count >= 3 && !memcmp(magic, "BZh", 3) ? CRPX_FILE_BZIP2 :
             count >= 6 && !memcmp(magic, "\xfd" "7zXZ\0", 6) ? CRPX_FILE_XZ : CRPX_FILE_RAW;
  }
  file->format = format;
  if (!crpx_file_format_available(format)) { fclose(file->raw); free(file); return CRPX_ERROR_UNSUPPORTED; }
  crpx_status status = initialise_codec(file);
  if (status != CRPX_OK) { fclose(file->raw); free(file); return status; }
  *out = file;
  return CRPX_OK;
}

/* One codec step. All library-specific pointer/size conventions stay here. */
static crpx_status codec_step(crpx_file *file, const unsigned char *input, size_t input_size,
                             unsigned char *output, size_t output_size, bool finish,
                             size_t *consumed, size_t *produced, bool *ended)
{
  (void)input; (void)output; (void)finish;
  *consumed = *produced = 0; *ended = false;
  if (input_size > UINT_MAX) input_size = UINT_MAX;
  if (output_size > UINT_MAX) output_size = UINT_MAX;
  switch (file->format) {
#ifdef HAVE_ZLIB
    case CRPX_FILE_GZIP: {
      file->gzip.next_in = (Bytef *)input; file->gzip.avail_in = (uInt)input_size;
      file->gzip.next_out = output; file->gzip.avail_out = (uInt)output_size;
      int code = file->writing ? deflate(&file->gzip, finish ? Z_FINISH : Z_NO_FLUSH) : inflate(&file->gzip, Z_NO_FLUSH);
      *consumed = input_size - file->gzip.avail_in; *produced = output_size - file->gzip.avail_out;
      *ended = code == Z_STREAM_END;
      if (code == Z_MEM_ERROR) return CRPX_ERROR_ALLOCATION;
      return code == Z_OK || code == Z_BUF_ERROR || *ended ? CRPX_OK : CRPX_ERROR_IO;
    }
#endif
#ifdef HAVE_BZIP2
    case CRPX_FILE_BZIP2: {
      file->bzip2.next_in = (char *)input; file->bzip2.avail_in = (unsigned int)input_size;
      file->bzip2.next_out = (char *)output; file->bzip2.avail_out = (unsigned int)output_size;
      int code = file->writing ? BZ2_bzCompress(&file->bzip2, finish ? BZ_FINISH : BZ_RUN) : BZ2_bzDecompress(&file->bzip2);
      *consumed = input_size - file->bzip2.avail_in; *produced = output_size - file->bzip2.avail_out;
      *ended = code == BZ_STREAM_END;
      if (code == BZ_MEM_ERROR) return CRPX_ERROR_ALLOCATION;
      return code == BZ_OK || code == BZ_RUN_OK || code == BZ_FINISH_OK || *ended ? CRPX_OK : CRPX_ERROR_IO;
    }
#endif
#ifdef HAVE_LZMA
    case CRPX_FILE_XZ: {
      file->xz.next_in = input; file->xz.avail_in = input_size;
      file->xz.next_out = output; file->xz.avail_out = output_size;
      lzma_ret code = lzma_code(&file->xz, finish ? LZMA_FINISH : LZMA_RUN);
      *consumed = input_size - file->xz.avail_in; *produced = output_size - file->xz.avail_out;
      *ended = code == LZMA_STREAM_END;
      if (code == LZMA_MEM_ERROR) return CRPX_ERROR_ALLOCATION;
      return code == LZMA_OK || *ended ? CRPX_OK : CRPX_ERROR_IO;
    }
#endif
    default: return CRPX_ERROR_UNSUPPORTED;
  }
}

static crpx_status fill_input(crpx_file *file)
{
  if (file->input_position < file->input_size || file->input_eof) return CRPX_OK;
  file->input_position = 0;
  file->input_size = fread(file->input, 1, sizeof(file->input), file->raw);
  if (ferror(file->raw)) return CRPX_ERROR_IO;
  if (!file->input_size) file->input_eof = true;
  return CRPX_OK;
}

static crpx_status read_backend(crpx_file *file, void *buffer, size_t capacity, size_t *count)
{
  *count = 0;
  if (file->failure != CRPX_OK) return file->failure;
  if (!capacity || file->read_finished) return CRPX_OK;
  unsigned char *output = buffer;
  if (file->format == CRPX_FILE_RAW) {
    size_t pending = file->input_size - file->input_position;
    if (pending > capacity) pending = capacity;
    if (pending) memcpy(output, file->input + file->input_position, pending);
    file->input_position += pending; *count = pending;
    if (pending < capacity) *count += fread(output + pending, 1, capacity - pending, file->raw);
    if (ferror(file->raw)) return file->failure = CRPX_ERROR_IO;
    return CRPX_OK;
  }
  while (*count < capacity && !file->read_finished) {
    crpx_status status = fill_input(file);
    if (status != CRPX_OK) return file->failure = status;
    if (file->member_end) {
      if (file->input_position == file->input_size && file->input_eof) { file->read_finished = true; break; }
      destroy_codec(file);
#ifdef HAVE_ZLIB
      if (file->format == CRPX_FILE_GZIP) memset(&file->gzip, 0, sizeof(file->gzip));
#endif
#ifdef HAVE_BZIP2
      if (file->format == CRPX_FILE_BZIP2) memset(&file->bzip2, 0, sizeof(file->bzip2));
#endif
      status = initialise_codec(file);
      if (status != CRPX_OK) return file->failure = status;
      file->member_end = false;
    }
    size_t consumed, produced;
    bool ended;
    status = codec_step(file, file->input + file->input_position, file->input_size - file->input_position,
                        output + *count, capacity - *count, file->input_eof, &consumed, &produced, &ended);
    file->input_position += consumed; *count += produced;
    if (status != CRPX_OK) return file->failure = status;
    if (ended) {
      if (file->format == CRPX_FILE_XZ) file->read_finished = true;
      else file->member_end = true;
    } else if (!consumed && !produced) return file->failure = CRPX_ERROR_IO;
  }
  return CRPX_OK;
}

crpx_status crpx_file_read(crpx_file *file, void *buffer, size_t capacity, size_t *read_size)
{
  if (!file || (!buffer && capacity) || !read_size || file->writing) return CRPX_ERROR_INVALID_ARGUMENT;
  *read_size = 0;
  if (file->failure != CRPX_OK) return file->failure;
  size_t buffered = file->line_size - file->line_position;
  if (buffered) {
    if (buffered > capacity) buffered = capacity;
    if (buffered) memcpy(buffer, file->lines + file->line_position, buffered);
    file->line_position += buffered; *read_size = buffered;
    return CRPX_OK;
  }
  return read_backend(file, buffer, capacity, read_size);
}

static crpx_status encode(crpx_file *file, const unsigned char *input, size_t length, bool finish)
{
  unsigned char output[8192];
  bool ended = false;
  do {
    size_t consumed, produced;
    crpx_status status = codec_step(file, input, length, output, sizeof(output), finish, &consumed, &produced, &ended);
    if (status != CRPX_OK) return file->failure = status;
    if (produced && fwrite(output, 1, produced, file->raw) != produced) return file->failure = CRPX_ERROR_IO;
    if (consumed) { input += consumed; length -= consumed; }
    if (!consumed && !produced && !ended && (length || finish)) return file->failure = CRPX_ERROR_IO;
  } while (length || (finish && !ended));
  return CRPX_OK;
}

crpx_status crpx_file_write(crpx_file *file, const void *buffer, size_t length)
{
  if (!file || (!buffer && length) || !file->writing) return CRPX_ERROR_INVALID_ARGUMENT;
  if (file->failure != CRPX_OK) return file->failure;
  if (!length) return CRPX_OK;
  if (file->format == CRPX_FILE_RAW) {
    if (fwrite(buffer, 1, length, file->raw) != length) return file->failure = CRPX_ERROR_IO;
    return CRPX_OK;
  }
  return encode(file, buffer, length, false);
}

crpx_status crpx_file_close(crpx_file *file)
{
  if (!file) return CRPX_ERROR_INVALID_ARGUMENT;
  crpx_status status = file->failure;
  if (file->writing && file->format != CRPX_FILE_RAW && status == CRPX_OK) status = encode(file, NULL, 0, true);
  destroy_codec(file);
  if (fclose(file->raw) && status == CRPX_OK) status = CRPX_ERROR_IO;
  free(file); return status;
}

crpx_status crpx_file_getline(crpx_file *file, char **line, size_t *capacity, size_t *length)
{
  if (!file || !line || !capacity || !length || file->writing || (!*line && *capacity))
    return CRPX_ERROR_INVALID_ARGUMENT;
  *length = 0;
  for (;;) {
    if (file->line_position == file->line_size) {
      crpx_status status = read_backend(file, file->lines, sizeof(file->lines), &file->line_size);
      file->line_position = 0;
      if (status != CRPX_OK) return status;
      if (!file->line_size) return *length ? CRPX_OK : CRPX_ERROR_NOT_FOUND;
    }
    if (*length >= SIZE_MAX - 1) return CRPX_ERROR_OVERFLOW;
    if (*capacity < *length + 2) {
      size_t grown = *capacity && *capacity <= SIZE_MAX / 2 ? *capacity * 2 : *length + 2;
      if (grown < 256) grown = 256;
      char *replacement = realloc(*line, grown);
      if (!replacement) return CRPX_ERROR_ALLOCATION;
      *line = replacement; *capacity = grown;
    }
    char character = (char)file->lines[file->line_position++];
    (*line)[(*length)++] = character; (*line)[*length] = '\0';
    if (character == '\n') return CRPX_OK;
  }
}

crpx_file_format crpx_file_get_format(const crpx_file *file) { return file ? file->format : CRPX_FILE_AUTO; }
