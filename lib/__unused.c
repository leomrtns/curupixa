uint64_t
crpx_brent_seed64 (void *vstate)
{// https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c
  uint64_t *v = (uint64_t *) vstate;
 (*v) ^= (*v) << 10; (*v) ^= (*v) >> 15; (*v) ^= (*v) << 4;  (*v) ^= (*v) >> 13;
  return *v;
}

uint64_t
crpx_std61_seed64 (void *vstate) 
{// https://github.com/quadram-institute-bioscience/biomcmc-lib/blob/master/lib/random_number_gen.c
  uint64_t *v = (uint64_t *) vstate;
  (*v) = ((*v) >> 31) + ((*v) << 30) - ((*v) >> 42) - ((*v) << 19);
  return *v;
}

