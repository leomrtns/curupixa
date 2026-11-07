/* Adapted from genefam-dist/lib/topology_splitset.c and guenomu.
 * Copyright (C) 2009-2026 Leonardo de Oliveira Martins
 * SPDX-License-Identifier: GPL-3.0-or-later */
#include "spr_distance.h"
#include "hashtable.h"
#include "lowlevel_algorithms.h"
#include <stdlib.h>
#include <string.h>

typedef struct { uint64_t *words; size_t ones, nwords; } spr_row;
typedef struct { spr_row *rows; uint64_t *storage; size_t count, capacity, stride; } spr_set;
struct crpx_spr_workspace {
  spr_set a,b,agree,disagree;
  double *costs;
  size_t *assignment, *indices;
  size_t assignment_capacity, index_capacity, bits;
  uint64_t *scratch;
  size_t scratch_words;
};
static size_t popcount(uint64_t x)
{
  /* Portable constant-work population count; no target-specific CPU requirement. */
  x-=(x>>1)&UINT64_C(0x5555555555555555);
  x=(x&UINT64_C(0x3333333333333333))+((x>>2)&UINT64_C(0x3333333333333333));
  x=(x+(x>>4))&UINT64_C(0x0f0f0f0f0f0f0f0f);
  return (size_t)((x*UINT64_C(0x0101010101010101))>>56);
}
static uint64_t mask(size_t bits) { return bits%64 ? (UINT64_C(1)<<(bits%64))-1 : UINT64_MAX; }
static size_t words_for(size_t bits) { return bits/64+(bits%64 != 0); }
static void set_free(spr_set *s) { free(s->rows); free(s->storage); memset(s,0,sizeof(*s)); }
void crpx_spr_workspace_clear(crpx_spr_workspace *w)
{
  if (!w) return;
  set_free(&w->a); set_free(&w->b); set_free(&w->agree); set_free(&w->disagree);
  free(w->costs); free(w->assignment); free(w->indices); free(w->scratch); memset(w,0,sizeof(*w));
}
crpx_status crpx_spr_workspace_new(crpx_spr_workspace **out)
{
  if (!out) return CRPX_ERROR_INVALID_ARGUMENT;
  *out=calloc(1,sizeof(**out)); return *out ? CRPX_OK : CRPX_ERROR_ALLOCATION;
}
void crpx_spr_workspace_free(crpx_spr_workspace *w) { if (w) { crpx_spr_workspace_clear(w); free(w); } }
size_t crpx_spr_workspace_storage_size(const crpx_spr_workspace *w)
{
  if (!w) return 0;
  const spr_set *sets[]={&w->a,&w->b,&w->agree,&w->disagree}; size_t n=sizeof(*w);
  for (size_t i=0;i<4;++i) n+=sets[i]->capacity*(sizeof(spr_row)+sets[i]->stride*sizeof(uint64_t));
  return n+w->assignment_capacity*w->assignment_capacity*sizeof(double)+
    (w->assignment_capacity+w->index_capacity)*sizeof(size_t)+w->scratch_words*sizeof(uint64_t);
}
static crpx_status set_reserve(spr_set *s,size_t rows,size_t words)
{
  if (rows<=s->capacity && words<=s->stride) { s->count=0; return CRPX_OK; }
  size_t capacity=rows>s->capacity?rows:s->capacity, stride=words>s->stride?words:s->stride;
  if (capacity>SIZE_MAX/sizeof(spr_row) || (stride && capacity>SIZE_MAX/stride/sizeof(uint64_t)))
    return CRPX_ERROR_OVERFLOW;
  spr_row *r=capacity?calloc(capacity,sizeof(*r)):NULL;
  uint64_t *storage=capacity&&stride?calloc(capacity*stride,sizeof(*storage)):NULL;
  if (capacity && (!r || !storage)) { free(r); free(storage); return CRPX_ERROR_ALLOCATION; }
  set_free(s); s->rows=r; s->storage=storage; s->capacity=capacity; s->stride=stride;
  for (size_t i=0;i<capacity;++i) r[i].words=storage+i*stride;
  return CRPX_OK;
}
static void normalise(spr_row *r,size_t bits)
{
  r->nwords=words_for(bits); r->ones=0;
  if (!r->nwords) return;
  r->words[r->nwords-1]&=mask(bits);
  for (size_t i=0;i<r->nwords;++i) r->ones+=popcount(r->words[i]);
  /* Match the reference orientation: equal-size sides retain bit zero. */
  if (r->ones>bits-r->ones || (r->ones==bits-r->ones && !(r->words[0]&1))) {
    for (size_t i=0;i<r->nwords;++i) r->words[i]=~r->words[i];
    r->words[r->nwords-1]&=mask(bits); r->ones=bits-r->ones;
  }
}
static bool equal(const spr_row *a,const spr_row *b)
{ return a->ones==b->ones && !memcmp(a->words,b->words,a->nwords*sizeof(uint64_t)); }
static int compare_rows(const void *left,const void *right)
{
  const spr_row *a=left,*b=right;
  if (a->ones!=b->ones) return a->ones>b->ones?1:-1;
  for (size_t i=a->nwords;i>0;--i) if (a->words[i-1]!=b->words[i-1]) return a->words[i-1]>b->words[i-1]?1:-1;
  return 0;
}
static void row_copy(spr_row *to,const spr_row *from)
{ memcpy(to->words,from->words,from->nwords*sizeof(uint64_t)); to->nwords=from->nwords; to->ones=from->ones; }
static void erase(spr_set *s,size_t i)
{ spr_row tmp=s->rows[i]; s->rows[i]=s->rows[--s->count]; s->rows[s->count]=tmp; }
static void simplify(spr_set *s,size_t bits)
{
  for (size_t i=0;i<s->count;) { normalise(s->rows+i,bits); if (s->rows[i].ones<2) erase(s,i); else ++i; }
  if (s->count>1) qsort(s->rows,s->count,sizeof(spr_row),compare_rows);
  /* Do not overwrite row descriptors: each owns a distinct slice of storage. */
  for (size_t i=1;i<s->count;) {
    if (equal(s->rows+i,s->rows+i-1)) {
      spr_row removed=s->rows[i]; memmove(s->rows+i,s->rows+i+1,(s->count-i-1)*sizeof(spr_row));
      s->rows[--s->count]=removed;
    } else ++i;
  }
}
static void agreement(crpx_spr_workspace *w)
{
  for (size_t i=0;i<w->a.count;) {
    size_t j=0; for (;j<w->b.count;++j) if (equal(w->a.rows+i,w->b.rows+j)) break;
    if (j==w->b.count) { ++i; continue; }
    row_copy(w->agree.rows+w->agree.count++,w->a.rows+i); erase(&w->a,i); erase(&w->b,j);
  }
}
static void bit_set(spr_row *r,size_t i,bool value)
{
  uint64_t bit=UINT64_C(1)<<(i%64);
  if (value) r->words[i/64]|=bit; else r->words[i/64]&=~bit;
}
static bool bit_get(const spr_row *r,size_t i) { return (r->words[i/64]>>(i%64))&1; }
static void replace_bit(crpx_spr_workspace *w,size_t to,size_t from)
{
  spr_set *sets[]={&w->a,&w->b,&w->agree};
  for (size_t k=0;k<3;++k) for (size_t i=0;i<sets[k]->count;++i)
    bit_set(sets[k]->rows+i,to,bit_get(sets[k]->rows+i,from));
}
static void resize_sets(crpx_spr_workspace *w)
{ simplify(&w->a,w->bits); simplify(&w->b,w->bits); simplify(&w->agree,w->bits); }
static void compress_agreement(crpx_spr_workspace *w)
{
  for (;;) {
    size_t i=0; for (;i<w->agree.count;++i) if (w->agree.rows[i].ones==2) break;
    if (i==w->agree.count) break;
    size_t first=SIZE_MAX,second=SIZE_MAX;
    for (size_t bit=0;bit<w->bits;++bit) if (bit_get(w->agree.rows+i,bit)) {
      if (first==SIZE_MAX) first=bit; else { second=bit; break; }
    }
    if (second<w->bits-1) replace_bit(w,second,w->bits-1);
    --w->bits; resize_sets(w);
  }
}
static crpx_status assign(crpx_spr_workspace *w,uint64_t *score)
{
  size_t n=w->a.count>w->b.count?w->a.count:w->b.count;
  *score=0; w->disagree.count=0;
  if (!n) return CRPX_OK;
  if (n>SIZE_MAX/n || n*n>SIZE_MAX/sizeof(double) || (long double)n*(w->bits/2)>9007199254740992.0L)
    return CRPX_ERROR_OVERFLOW;
  if (n>w->assignment_capacity) {
    double *costs=malloc(n*n*sizeof(double)); size_t *a=malloc(n*sizeof(size_t));
    if (!costs || !a) { free(costs); free(a); return CRPX_ERROR_ALLOCATION; }
    free(w->costs); free(w->assignment); w->costs=costs; w->assignment=a; w->assignment_capacity=n;
  }
  memset(w->costs,0,n*n*sizeof(double));
  for (size_t i=0;i<w->a.count;++i) for (size_t j=0;j<w->b.count;++j) {
    size_t ones=0;
    for (size_t k=0;k<words_for(w->bits);++k) ones+=popcount(w->a.rows[i].words[k]^w->b.rows[j].words[k]);
    if (w->bits-ones<ones) ones=w->bits-ones;
    w->costs[i*n+j]=(double)ones;
  }
  double total; crpx_status status=crpx_assignment_solve(w->costs,n,n,w->assignment,&total);
  if (status!=CRPX_OK) return status;
  *score=(uint64_t)total;
  for (size_t i=0;i<w->a.count;++i) if (w->assignment[i]<w->b.count) {
    spr_row *r=w->disagree.rows+w->disagree.count++;
    for (size_t k=0;k<words_for(w->bits);++k)
      r->words[k]=w->a.rows[i].words[k]^w->b.rows[w->assignment[i]].words[k];
    normalise(r,w->bits);
  }
  if (w->disagree.count>1) qsort(w->disagree.rows,w->disagree.count,sizeof(spr_row),compare_rows);
  for (size_t i=1;i<w->disagree.count;) {
    if (equal(w->disagree.rows+i,w->disagree.rows+i-1)) {
      spr_row removed=w->disagree.rows[i];
      memmove(w->disagree.rows+i,w->disagree.rows+i+1,(w->disagree.count-i-1)*sizeof(spr_row));
      w->disagree.rows[--w->disagree.count]=removed;
    } else ++i;
  }
  return CRPX_OK;
}
static bool subset(const spr_row *container,const spr_row *part,bool complement,size_t bits)
{
  for (size_t i=0;i<part->nwords;++i) {
    uint64_t word=complement?~container->words[i]:container->words[i];
    if (i+1==part->nwords) word&=mask(bits);
    if ((word&part->words[i])!=part->words[i]) return false;
  }
  return true;
}
static crpx_status load_set(spr_set *to,const crpx_splitset *from,size_t bits)
{
  for (size_t i=0;i<crpx_splitset_count(from);++i) {
    const struct crpx_bipartition_struct *b=crpx_splitset_get(from,i);
    if (crpx_bipartition_size(b)!=bits) return CRPX_ERROR_INVALID_ARGUMENT;
    size_t words; const uint64_t *source=crpx_bipartition_words(b,&words);
    spr_row *r=to->rows+to->count++; memcpy(r->words,source,words*sizeof(uint64_t)); normalise(r,bits);
  }
  simplify(to,bits); return CRPX_OK;
}
crpx_status crpx_splitset_spr_distance(const crpx_splitset *a,const crpx_splitset *b,
                                      crpx_spr_workspace *w,crpx_spr_result *out)
{
  if (!a || !b || !w || !out) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t bits=crpx_splitset_taxa(a), ac=crpx_splitset_count(a),bc=crpx_splitset_count(b);
  if (bits!=crpx_splitset_taxa(b)) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t words=words_for(bits), smaller=ac<bc?ac:bc;
  crpx_status status=set_reserve(&w->a,ac,words);
  if (status==CRPX_OK) status=set_reserve(&w->b,bc,words);
  if (status==CRPX_OK) status=set_reserve(&w->agree,smaller,words);
  if (status==CRPX_OK) status=set_reserve(&w->disagree,smaller,words);
  if (status!=CRPX_OK) return status;
  if (bits>w->index_capacity) {
    if (bits>SIZE_MAX/sizeof(size_t)) return CRPX_ERROR_OVERFLOW;
    size_t *indices=malloc(bits*sizeof(size_t)); if (!indices) return CRPX_ERROR_ALLOCATION;
    free(w->indices); w->indices=indices; w->index_capacity=bits;
  }
  if (words>w->scratch_words) {
    if (words>SIZE_MAX/sizeof(uint64_t)) return CRPX_ERROR_OVERFLOW;
    uint64_t *scratch=malloc(words*sizeof(uint64_t)); if (!scratch) return CRPX_ERROR_ALLOCATION;
    free(w->scratch); w->scratch=scratch; w->scratch_words=words;
  }
  w->bits=bits;
  status=load_set(&w->a,a,bits); if (status==CRPX_OK) status=load_set(&w->b,b,bits);
  if (status!=CRPX_OK) return status;
  agreement(w);
  crpx_spr_result result={0};
  if (w->a.count>SIZE_MAX-w->b.count) return CRPX_ERROR_OVERFLOW;
  result.rf=w->a.count+w->b.count;
  status=assign(w,&result.hdist); if (status!=CRPX_OK) return status;
  /* Rebuild: agreement compression must retain the full trees' common edges. */
  w->a.count=w->b.count=w->agree.count=0;
  status=load_set(&w->a,a,bits); if (status==CRPX_OK) status=load_set(&w->b,b,bits);
  if (status!=CRPX_OK) return status;
  bool first_assignment=true;
  for (;;) {
    agreement(w); compress_agreement(w);
    if (!w->a.count || !w->b.count) break;
    uint64_t cost; status=assign(w,&cost); if (status!=CRPX_OK) return status;
    if (first_assignment) { result.reduced_hdist=cost; first_assignment=false; }
    if (!w->disagree.count) return CRPX_ERROR_TOPOLOGY;
    spr_row prune={w->scratch,0,words_for(w->bits)}; row_copy(&prune,w->disagree.rows);
    if (prune.ones>1) {
      bool found=false;
      for (size_t i=0;i<w->disagree.count && !found;++i) for (size_t j=0;j<w->agree.count;++j)
        if (equal(w->disagree.rows+i,w->agree.rows+j)) { row_copy(&prune,w->disagree.rows+i); found=true; break; }
      for (size_t i=0;i<w->a.count;++i) if (!subset(w->a.rows+i,&prune,false,w->bits) &&
          !subset(w->a.rows+i,&prune,true,w->bits)) { ++result.spr_extra; break; }
    }
    if (!prune.ones || prune.ones>=w->bits) return CRPX_ERROR_TOPOLOGY;
    size_t n=w->bits-prune.ones, used=0;
    for (size_t bit=0;bit<w->bits;++bit) if (bit_get(&prune,bit)) w->indices[used++]=bit;
    /* Reference tail-swap restriction preserves its deterministic tie behaviour. */
    size_t low=0, high=prune.ones;
    for (size_t tail=w->bits;tail>n;) {
      --tail;
      if (w->indices[low]>=n) break;
      if (high && tail==w->indices[high-1]) --high;
      else replace_bit(w,w->indices[low++],tail);
    }
    w->bits=n; ++result.spr; resize_sets(w);
  }
  *out=result; return CRPX_OK;
}
/* Tree adapters share explicit taxon mapping; no analysis state is attached to nodes. */
crpx_status crpx_tree_spr_distance(const crpx_tree *a,const crpx_tree *b,
                                  crpx_spr_workspace *w,crpx_spr_result *result)
{
  if (!a || !b || !w || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t taxa=crpx_tree_leaf_count(a), n=0;
  if (!taxa || taxa!=crpx_tree_leaf_count(b)) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa>SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_node_id *left=malloc(taxa*sizeof(*left)),*right=malloc(taxa*sizeof(*right));
  crpx_hashtable_t labels=NULL; crpx_splitset *x=NULL,*y=NULL;
  if (!left || !right) { free(left); free(right); return CRPX_ERROR_ALLOCATION; }
  crpx_status status=crpx_hashtable_new(taxa,0,&labels);
  for (size_t i=0;i<crpx_tree_node_count(b) && status==CRPX_OK;++i)
    if (crpx_tree_first_child(b,(crpx_node_id)i)==CRPX_NODE_NONE) {
      const char *name=crpx_tree_label(b,(crpx_node_id)i); bool replaced;
      if (!name || !*name) { status=CRPX_ERROR_INVALID_ARGUMENT; break; }
      status=crpx_hashtable_put_string(labels,name,i,&replaced);
      if (status==CRPX_OK && replaced) { status=CRPX_ERROR_INVALID_ARGUMENT; break; }
    }
  crpx_hashtable_t seen=NULL;
  if (status==CRPX_OK) status=crpx_hashtable_new(taxa,0,&seen);
  for (size_t i=0;i<crpx_tree_node_count(a) && status==CRPX_OK;++i)
    if (crpx_tree_first_child(a,(crpx_node_id)i)==CRPX_NODE_NONE) {
      const char *name=crpx_tree_label(a,(crpx_node_id)i); uint64_t mapped; bool replaced;
      if (!name || !*name) { status=CRPX_ERROR_INVALID_ARGUMENT; break; }
      status=crpx_hashtable_get_string(labels,name,&mapped);
      if (status==CRPX_ERROR_NOT_FOUND) status=CRPX_ERROR_INVALID_ARGUMENT;
      if (status!=CRPX_OK) break;
      status=crpx_hashtable_put_string(seen,name,i,&replaced);
      if (status!=CRPX_OK || replaced) { if (status==CRPX_OK) status=CRPX_ERROR_INVALID_ARGUMENT; break; }
      left[n]=(crpx_node_id)i; right[n++]=(crpx_node_id)mapped;
    }
  if (status==CRPX_OK) status=crpx_splitset_from_tree(a,left,taxa,&x);
  if (status==CRPX_OK) status=crpx_splitset_from_tree(b,right,taxa,&y);
  if (status==CRPX_OK) status=crpx_splitset_spr_distance(x,y,w,result);
  crpx_splitset_free(x); crpx_splitset_free(y); crpx_hashtable_free(labels); crpx_hashtable_free(seen);
  free(left); free(right); return status;
}
crpx_status crpx_gene_species_spr_distance(const crpx_tree *gene,const crpx_tree *species,
                                          const crpx_node_id *map,crpx_spr_workspace *w,crpx_spr_result *result)
{
  if (!gene || !species || !map || !w || !result) return CRPX_ERROR_INVALID_ARGUMENT;
  size_t taxa=crpx_tree_leaf_count(gene),n=0;
  if (!taxa) return CRPX_ERROR_INVALID_ARGUMENT;
  if (taxa>SIZE_MAX/sizeof(crpx_node_id)) return CRPX_ERROR_OVERFLOW;
  crpx_node_id *a=malloc(taxa*sizeof(*a)),*b=malloc(taxa*sizeof(*b));
  if (!a || !b) { free(a); free(b); return CRPX_ERROR_ALLOCATION; }
  for (size_t i=0;i<crpx_tree_node_count(gene);++i)
    if (crpx_tree_first_child(gene,(crpx_node_id)i)==CRPX_NODE_NONE) { a[n]=(crpx_node_id)i; b[n++]=map[i]; }
  crpx_splitset *x=NULL,*y=NULL; crpx_status status=crpx_splitset_from_tree(gene,a,taxa,&x);
  if (status==CRPX_OK) status=crpx_splitset_from_tree(species,b,taxa,&y);
  if (status==CRPX_OK) status=crpx_splitset_spr_distance(x,y,w,result);
  crpx_splitset_free(x); crpx_splitset_free(y); free(a); free(b); return status;
}
