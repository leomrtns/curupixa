/* Copyright (C) 2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "file_compression.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(expression) do { if (!(expression)) { \
  fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, #expression); return 1; \
} } while (0)

static const char *one = "curupixa-test-compression-one.bin";
static const char *two = "curupixa-test-compression-two.bin";
static const char *altered = "curupixa-test-compression-altered.bin";

static int compressed_bytes(unsigned char **result, size_t *length)
{
  FILE *file = fopen(one, "rb");
  CHECK(file);
  CHECK(fseek(file, 0, SEEK_END) == 0);
  long count = ftell(file);
  CHECK(count >= 0 && fseek(file, 0, SEEK_SET) == 0);
  unsigned char *bytes = malloc((size_t)count + 1);
  CHECK(bytes && fread(bytes, 1, (size_t)count, file) == (size_t)count && fclose(file) == 0);
  *result = bytes; *length = (size_t)count;
  return 0;
}

static int read_data(const char *path, crpx_file_format format, const unsigned char *expected, size_t length)
{
  crpx_file *file;
  CHECK(crpx_file_open(path, "r", format, &file) == CRPX_OK);
  unsigned char buffer[997];
  size_t position = 0, count;
  for (;;) {
    CHECK(crpx_file_read(file, buffer, sizeof(buffer), &count) == CRPX_OK);
    if (!count) break;
    CHECK(count <= length - position && !memcmp(buffer, expected + position, count));
    position += count;
  }
  CHECK(position == length && crpx_file_close(file) == CRPX_OK);
  return 0;
}

static int read_failure(crpx_file_format format)
{
  crpx_file *file;
  CHECK(crpx_file_open(altered, "r", format, &file) == CRPX_OK);
  unsigned char buffer[1000];
  size_t count;
  crpx_status status;
  do { status = crpx_file_read(file, buffer, sizeof(buffer), &count); } while (status == CRPX_OK && count);
  CHECK(status == CRPX_ERROR_IO);
  CHECK(crpx_file_close(file) == CRPX_ERROR_IO);
  return 0;
}

static int test_format(crpx_file_format format)
{
  if (!crpx_file_format_available(format)) return 0;
  const size_t length = 100003;
  unsigned char *data = malloc(length * 2);
  CHECK(data);
  for (size_t i = 0; i < length; i++) data[i] = (unsigned char)((i * 33) ^ (i >> 4));
  memcpy(data + length, data, length);
  crpx_file *file;
  CHECK(crpx_file_open(one, "w", format, &file) == CRPX_OK);
  CHECK(crpx_file_write(file, data, 1) == CRPX_OK);
  CHECK(crpx_file_write(file, data + 1, 8192) == CRPX_OK);
  CHECK(crpx_file_write(file, data + 8193, length - 8193) == CRPX_OK);
  CHECK(crpx_file_close(file) == CRPX_OK);
  CHECK(read_data(one, CRPX_FILE_AUTO, data, length) == 0);
  CHECK(read_data(one, format, data, length) == 0);
  if (format != CRPX_FILE_RAW) {
    unsigned char *bytes;
    size_t byte_count;
    CHECK(compressed_bytes(&bytes, &byte_count) == 0 && byte_count > 3);
    FILE *raw = fopen(two, "wb");
    CHECK(raw && fwrite(bytes, 1, byte_count, raw) == byte_count);
    CHECK(fwrite(bytes, 1, byte_count, raw) == byte_count && fclose(raw) == 0);
    CHECK(read_data(two, CRPX_FILE_AUTO, data, length * 2) == 0);
    if (format == CRPX_FILE_XZ) {
      const unsigned char padding[4] = {0, 0, 0, 0};
      raw = fopen(two, "wb");
      CHECK(raw && fwrite(bytes, 1, byte_count, raw) == byte_count);
      CHECK(fwrite(padding, 1, sizeof(padding), raw) == sizeof(padding));
      CHECK(fwrite(bytes, 1, byte_count, raw) == byte_count);
      CHECK(fwrite(padding, 1, sizeof(padding), raw) == sizeof(padding) && fclose(raw) == 0);
      CHECK(read_data(two, CRPX_FILE_AUTO, data, length * 2) == 0);
    }
    raw = fopen(altered, "wb");
    CHECK(raw && fwrite(bytes, 1, byte_count - 3, raw) == byte_count - 3 && fclose(raw) == 0);
    CHECK(read_failure(format) == 0);
    raw = fopen(altered, "wb");
    CHECK(raw && fwrite(bytes, 1, byte_count, raw) == byte_count);
    CHECK(fwrite("trailing garbage", 1, 16, raw) == 16 && fclose(raw) == 0);
    CHECK(read_failure(format) == 0);
    bytes[byte_count / 2] ^= 0x40;
    raw = fopen(altered, "wb");
    CHECK(raw && fwrite(bytes, 1, byte_count, raw) == byte_count && fclose(raw) == 0);
    CHECK(read_failure(format) == 0);
    free(bytes);
  }
  const char text[] = "first\nsecond line\nlast";
  CHECK(crpx_file_open(one, "w", format, &file) == CRPX_OK);
  CHECK(crpx_file_write(file, text, sizeof(text) - 1) == CRPX_OK && crpx_file_close(file) == CRPX_OK);
  CHECK(crpx_file_open(one, "r", CRPX_FILE_AUTO, &file) == CRPX_OK);
  char *line = NULL; size_t capacity = 0, count;
  CHECK(crpx_file_getline(file, &line, &capacity, &count) == CRPX_OK && count == 6 && !strcmp(line, "first\n"));
  char prefix[3];
  CHECK(crpx_file_read(file, prefix, 3, &count) == CRPX_OK && count == 3 && !memcmp(prefix, "sec", 3));
  CHECK(crpx_file_getline(file, &line, &capacity, &count) == CRPX_OK && !strcmp(line, "ond line\n"));
  CHECK(crpx_file_getline(file, &line, &capacity, &count) == CRPX_OK && !strcmp(line, "last"));
  CHECK(crpx_file_getline(file, &line, &capacity, &count) == CRPX_ERROR_NOT_FOUND);
  CHECK(crpx_file_close(file) == CRPX_OK);
  free(line);
  CHECK(crpx_file_open(one, "w", format, &file) == CRPX_OK && crpx_file_close(file) == CRPX_OK);
  CHECK(read_data(one, CRPX_FILE_AUTO, data, 0) == 0);
  free(data);
  return 0;
}

int main(void)
{
  unsigned char sentinel;
  crpx_file *invalid = (crpx_file *)&sentinel;
  CHECK(crpx_file_open(NULL, "r", CRPX_FILE_AUTO, &invalid) == CRPX_ERROR_INVALID_ARGUMENT && !invalid);
  invalid = (crpx_file *)&sentinel;
  CHECK(crpx_file_open("unused", "invalid", CRPX_FILE_AUTO, &invalid) == CRPX_ERROR_INVALID_ARGUMENT && !invalid);
  for (int format = CRPX_FILE_RAW; format <= CRPX_FILE_XZ; format++) CHECK(test_format((crpx_file_format)format) == 0);
  remove(one); remove(two); remove(altered);
  return 0;
}
