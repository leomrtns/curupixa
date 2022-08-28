void 
crpx_generate_random_seed_256bits (crpx_global_t cglob, uint64_t seed[4]) /* UNUSED / pilot  */
{
  FILE *fp;
  uint64_t x1[4]; // unsigned long int
  unsigned long long int x = 0; // courtesy of 32bits machines, they're not the same
  unsigned int i, idx, success;
  fp = fopen ("/dev/urandom", "r");
  if (fp == NULL) crpx_logger_info (cglob, "Could not open /dev/urandom to include it in the random seed\n");
  else {
    if (fread (seed, sizeof (uint64_t), 4, fp) != 1) crpx_logger_info (cglob, "Could not read from /dev/urandom to include it in the random seed\n");
    fclose (fp);
  }
  crpx_get_time_128bits (x1);
  seed[0] ^= x1[0];
  seed[1] ^= x1[1];
  seed[2] ^= x1[0];
  seed[3] ^= x1[1];
#ifdef HAVE_RDRND
  for (i=64, idx=0; idx < 4;i--) { // DRNG suggests 10 cycles 
    success = __builtin_ia32_rdrand64_step (&x);
    x1[idx] = x;
    idx += (success > 0); // avoid mispredicted branches (not performance-critical here though) 
  }
#endif
  x = (uint64_t) getpid() + getppid(); // process ID and parent process ID (can be small numbers)
  for (i=0; i<4;i++) seed[i] ^= (x1[i] + ((uint64_t)(i + 1) * x));
  return;
}

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

inline uint64_t
crpx_hashint_zixmix64 (uint64_t h) // identical to fastmix64; in zix/digest it's used in word loop (as is fasthash() in curupixa)
{ // biomcmc and https://github.com/drobilla/zix/blob/main/src/digest.c
  h ^= h >> 23U; h *= 0x2127599BF4325C37ULL; h ^= h >> 47U;
  return h;
}

