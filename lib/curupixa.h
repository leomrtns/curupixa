/* 
 *This file is part of curupixa, a low-level library for phylogenomic analysis.
 * Copyright (C) 2022-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file curupixa.h 
 *  \brief headers exposed to other programs
 */

#ifndef _curupixa_toplevel_h_ 
#define _curupixa_toplevel_h_

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "global/global_variable.h"
#include "index_arrangement.h"
#include "quasi_random.h"
#include "status.h"
#include "char_vector.h"
#include "bipartition.h"
#include "hashtable.h"
#include "distance_matrix.h"
#include "tree.h"
#include "quickselect_quantile.h"
#include "empirical_frequency.h"
#include "file_compression.h"
#include "nexus_common.h"
#include "alignment.h"
#include "parsimony.h"
#include "reconciliation.h"
#include "ancestry_index.h"
#include "reconciliation_roots.h"
#include "taxon_map.h"
#include "spr_distance.h"
#include "lowlevel_algorithms.h"
#include "splitset_distances.h"
#include "topology_distance.h"
#include "topology_randomise.h"
#include "upgma.h"
#include "newick_space.h"
#include "prob_distribution.h"
#include "svd.h"
#include "distance_generator.h"
#include "clustering_goptics.h"
#include "gff3_format.h"
#include "hll.h"
#include "kmerhash.h"
#include "rolling_hash.h"
#include "suffix_tree.h"
#include "minimiser_simplex.h"
#include "character_data.h"
#include "substitution_model.h"
#include "likelihood.h"
#include "continuous_traits.h"
#include "reconstruction.h"
#include "character_analysis.h"

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* if header not defined */
