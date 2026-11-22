/* Copyright (C) 2026 Leonardo de Oliveira Martins
   SPDX-License-Identifier: GPL-3.0-or-later */
/*! \file tip_correspondence.h
    \brief Match complete sets of unique tip names once for ordinary tree comparisons. */
#ifndef CURUPIXA_TIP_CORRESPONDENCE_H
#define CURUPIXA_TIP_CORRESPONDENCE_H
#include "tree.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
  size_t taxa;
  crpx_node_id *first, *second;
} crpx_tip_correspondence;
/* Exact label correspondence. Both trees must contain the same nonempty unique
   tip names; unmatched tips cause an error, never automatic pruning. The output
   owns arrays and retains no tree pointers. Use explicit maps for gene copies. */
crpx_status crpx_tip_correspondence_new (const crpx_tree *first, const crpx_tree *second,
                                         crpx_tip_correspondence **out);
void crpx_tip_correspondence_free (crpx_tip_correspondence *correspondence);
#ifdef __cplusplus
}
#endif
#endif
