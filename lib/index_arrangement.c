/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: MIT */
/*! \file index_arrangement.c
    \brief Independent index enumeration and composition, preserving the established interface. */
#include "index_arrangement.h"
static bool step_permutation (crpx_index_permutation_t p, bool forward);
crpx_index_permutation_t
crpx_index_permutation_new (crpx_global_t context, size_t count)
{
  if (count > SIZE_MAX / sizeof (size_t)) return NULL;
  crpx_index_permutation_t p = calloc (1, sizeof (*p));
  if (!p) return NULL;
  p->idx = count ? malloc (count * sizeof (size_t)) : NULL;
  if (count && !p->idx) {
    free (p);
    return NULL;
  }
  p->size = count;
  p->cglob = crpx_global_retain (context);
  crpx_index_permutation_reset (p);
  return p;
}
void
del_crpx_index_permutation (crpx_index_permutation_t p)
{
  if (p) {
    free (p->idx);
    crpx_global_finalise (p->cglob);
    free (p);
  }
}
void
crpx_index_permutation_reset (crpx_index_permutation_t p)
{
  if (p)
    for (size_t i = 0; i < p->size; i++)
      p->idx[i] = i;
}
bool
crpx_index_permutation_swap (crpx_index_permutation_t p, size_t i, size_t j)
{
  if (!p || i >= p->size || j >= p->size) return false;
  crpx_index_permutation_swap_no_checks (p, i, j);
  return true;
}
void
crpx_index_permutation_swap_no_checks (crpx_index_permutation_t p, size_t i, size_t j)
{
  size_t value = p->idx[i];
  p->idx[i] = p->idx[j];
  p->idx[j] = value;
}
void
crpx_index_permutation_reverse (crpx_index_permutation_t p)
{
  if (p)
    for (size_t i = 0; i < p->size / 2; i++)
      crpx_index_permutation_swap_no_checks (p, i, p->size - 1 - i);
}
bool
crpx_index_permutation_inverse (crpx_index_permutation_t out, const crpx_index_permutation_t input)
{
  if (!out || !input || out->size != input->size || !crpx_index_permutation_is_valid (input)) return false;
  size_t *values = out->idx;
  if (out == input && out->size) {
    values = malloc (out->size * sizeof (size_t));
    if (!values) return false;
  }
  for (size_t i = 0; i < input->size; i++)
    values[input->idx[i]] = i;
  if (values != out->idx) {
    memcpy (out->idx, values, out->size * sizeof (size_t));
    free (values);
  }
  return true;
}
bool
crpx_index_permutation_combine (crpx_index_permutation_t out, const crpx_index_permutation_t first,
                                const crpx_index_permutation_t second)
{
  if (!out || !first || !second || out->size != first->size || out->size != second->size ||
      !crpx_index_permutation_is_valid (first) || !crpx_index_permutation_is_valid (second))
    return false;
  size_t *values = out->idx;
  if ((out == first || out == second) && out->size) {
    values = malloc (out->size * sizeof (size_t));
    if (!values) return false;
  }
  for (size_t i = 0; i < out->size; i++)
    values[i] = second->idx[first->idx[i]];
  if (values != out->idx) {
    memcpy (out->idx, values, out->size * sizeof (size_t));
    free (values);
  }
  return true;
}
crpx_index_combination_t
crpx_index_combination_new (crpx_global_t context, size_t n, size_t k)
{
  if (k > n || k > SIZE_MAX / sizeof (size_t)) return NULL;
  crpx_index_combination_t c = calloc (1, sizeof (*c));
  if (!c) return NULL;
  c->idx = k ? malloc (k * sizeof (size_t)) : NULL;
  if (k && !c->idx) {
    free (c);
    return NULL;
  }
  c->n = n;
  c->k = k;
  c->cglob = crpx_global_retain (context);
  crpx_index_combination_reset_first (c);
  return c;
}
void
del_crpx_index_combination (crpx_index_combination_t c)
{
  if (c) {
    free (c->idx);
    crpx_global_finalise (c->cglob);
    free (c);
  }
}
void
crpx_index_combination_reset_first (crpx_index_combination_t c)
{
  if (c)
    for (size_t i = 0; i < c->k; i++)
      c->idx[i] = i;
}
void
crpx_index_combination_reset_last (crpx_index_combination_t c)
{
  if (c)
    for (size_t i = 0; i < c->k; i++)
      c->idx[i] = c->n - c->k + i;
}
bool
crpx_index_permutation_is_valid (const crpx_index_permutation_t p)
{
  if (!p || (p->size && !p->idx)) return false;
  if (p->size <= 64) {
    uint64_t seen = 0;
    for (size_t i = 0; i < p->size; i++) {
      if (p->idx[i] >= p->size || (seen & (UINT64_C (1) << p->idx[i]))) return false;
      seen |= UINT64_C (1) << p->idx[i];
    }
    return true;
  }
  /* Large permutations use a membership set instead of comparing every pair. */
  unsigned char *seen = p->size ? calloc (p->size, 1) : NULL;
  if (p->size && !seen) return false;
  bool valid = true;
  for (size_t i = 0; i < p->size; i++) {
    if (p->idx[i] >= p->size || seen[p->idx[i]]) {
      valid = false;
      break;
    }
    seen[p->idx[i]] = 1;
  }
  free (seen);
  return valid;
}
bool
crpx_index_permutation_next (crpx_index_permutation_t p)
{
  return step_permutation (p, true);
}
bool
crpx_index_permutation_prev (crpx_index_permutation_t p)
{
  return step_permutation (p, false);
}
bool
crpx_index_combination_is_valid (crpx_index_combination_t c)
{
  if (!c || c->k > c->n || (c->k && !c->idx)) return false;
  for (size_t i = 0; i < c->k; i++)
    if (c->idx[i] >= c->n || (i && c->idx[i] <= c->idx[i - 1])) return false;
  return true;
}
bool
crpx_index_combination_next (crpx_index_combination_t c)
{
  if (!c || !c->k || c->k > c->n) return false;
  size_t pos = c->k;
  while (pos && c->idx[pos - 1] == c->n - c->k + pos - 1)
    pos--;
  if (!pos) return false;
  c->idx[--pos]++;
  for (size_t i = pos + 1; i < c->k; i++)
    c->idx[i] = c->idx[i - 1] + 1;
  return true;
}
bool
crpx_index_combination_prev (crpx_index_combination_t c)
{
  if (!c || !c->k || c->k > c->n) return false;
  size_t pos = c->k;
  while (pos && c->idx[pos - 1] == (pos == 1 ? 0 : c->idx[pos - 2] + 1))
    pos--;
  if (!pos) return false;
  c->idx[--pos]--;
  for (size_t i = pos + 1; i < c->k; i++)
    c->idx[i] = c->n - c->k + i;
  return true;
}
static bool
step_permutation (crpx_index_permutation_t p, bool forward)
{ /* The suffix is monotone: its rightmost admissible element is the successor. */
  if (!p || p->size < 2) return false;
  size_t pivot = p->size - 1;
  while (pivot && (forward ? p->idx[pivot - 1] >= p->idx[pivot] : p->idx[pivot - 1] <= p->idx[pivot]))
    pivot--;
  if (!pivot) return false;
  size_t successor = p->size - 1;
  while (forward ? p->idx[successor] <= p->idx[pivot - 1] : p->idx[successor] >= p->idx[pivot - 1])
    successor--;
  crpx_index_permutation_swap_no_checks (p, pivot - 1, successor);
  for (size_t last = p->size - 1; pivot < last; pivot++, last--)
    crpx_index_permutation_swap_no_checks (p, pivot, last);
  return true;
}
