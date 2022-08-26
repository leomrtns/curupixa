#include <curupixa.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

/* Testing available PRNGs with dieharder, command line below:
 * `./tests/dieharder_rng 7 | dieharder -g 200 -a`
 * notice that according to dieharder's man page, the _generator_ number 200 is stdin (thus "-g 200", not "-f 200")
 *
 * WARNING: this program will run for a _very_ long time 
 */
#define BUFFERSIZE 512

int main(int argc, char **argv)
{
  uint64_t i, ntries;
  uint64_t seed[16];
  uint64_t (*rng)(void*);
  uint32_t (*rng2)(void*);
  uint8_t algo;
  int j;
  char algoname[16] = {'\0'};

  if (argc == 1) return TEST_SKIPPED;
  crpx_global_t cglob = crpx_global_init (0,"debug");
  crpx_generate_bytesized_random_seeds_from_cpu (cglob, seed, 128);
  sscanf (argv[1], " %hhu ", &algo);

  if (algo < 32) {  
    switch (algo) { 
      case 0: rng = &crpx_rng_romu_seed256;       strcpy (algoname, "romu256"); break; // OK // fast+-
      case 1: rng = &crpx_rng_romu_seed192;       strcpy (algoname, "romu192"); break; // OK // fast
      case 2: rng = &crpx_rng_romu_seed128;       strcpy (algoname, "romu128"); break; // excellent // normal
      case 3: rng = &crpx_xoroshiro_pv6_seed128;  strcpy (algoname, "xoropv6"); break; // OK // fast  
      case 4: rng = &crpx_xoroshiro_pv8_seed128;  strcpy (algoname, "xoropv8"); break; // Ok // fast
      case 5: rng = &crpx_xoroshiro_pp_seed128;   strcpy (algoname, "xoropp1"); break; // OK // normal
      case 6: rng = &crpx_xoroshiro_pp_seed256;   strcpy (algoname, "xoropp2"); break; // +- // normal
      case 7: rng = &crpx_xoroshiro_star_seed256; strcpy (algoname, "xorost2"); break; // OK // normal
      case 8: rng = &crpx_xorshift_star_seed64;   strcpy (algoname, "shif64 "); break; // excellent // slow 
      case 9: rng = &crpx_xorshift_p_seed128;     strcpy (algoname, "shif128"); break; // OK // fast+-

      case 10: rng = &crpx_rng_wyhash_state64;  strcpy (algoname, "wyhash "); break; // excellent // fast+-
      case 11: rng = &crpx_rng_splitmix_seed64; strcpy (algoname, "splitmi"); break; // excellent // normal
      case 12: rng = &crpx_rng_rrmixer_seed64;  strcpy (algoname, "rrmixer"); break; // OK // normal
      case 13: rng = &crpx_rng_moremur_seed64;  strcpy (algoname, "moremur"); break; // excellent // normal
      case 14: rng = &crpx_rng_lehmer_seed128;  strcpy (algoname, "lehmer "); break; // excellent // fast+-
      case 15: rng = &crpx_rng_wyrand_seed64;   strcpy (algoname, "wyrand "); break; // +- // fastest
      case 16: rng = &crpx_rng_pcg_seed256;     strcpy (algoname, "pcg    "); break; // OK (depends on seed) // slow
      case 17: rng = &crpx_rng_jenkins13_seed256; strcpy (algoname, "64jen13"); break; // excellent // normal
      case 18: rng = &crpx_rng_jenkins19_seed256; strcpy (algoname, "64jen19"); break; // excellent // normal
      default: rng = &crpx_rng_wyrand_seed64;   strcpy (algoname, "wyrand "); break; 
    }
  } else {
    switch (algo) { // speed is comparable to 64bits (500mi~440mi numbers/second), which means half the throughput
      case 32: rng2 = &crpx_rng_abyssinian_seed128; strcpy (algoname, "abyssin"); break; // OK // 3rd fastest
      case 33: rng2 = &crps_rng_widynski_seed192;   strcpy (algoname, "widynsk"); break; // excellent // 4th fastest
      case 34: rng2 = &crpx_rng_jenkins8_seed128;   strcpy (algoname, "32jenk8"); break; // OK // 1st fastest (practrand found a couple "unusual")
      case 35: rng2 = &crpx_rng_jenkins13_seed128;  strcpy (algoname, "32jen13"); break; // excellent // 2nd fastest
      default: rng2 = &crpx_rng_jenkins8_seed128;   strcpy (algoname, "32jenk8"); break;
    }
  }

  for (i=0; i < 16; i++) {
    seed[i] |= 1ULL; // make it odd
    fprintf (stderr, "%17lx ", seed[i]); if (!((i+1)%4)) fprintf (stderr, "\n");
  }
  fprintf (stderr, "%s :  %lf seconds to set seed vector\n", algoname, crpx_update_elapsed_time_128bits (cglob->elapsed_time));

  if (argc == 2) { // dieharder test
    uint64_t x[BUFFERSIZE];
    uint32_t *x2 = (uint32_t*) x;
    for (i=0; i < UINT64_MAX; i++) {
      if (algo < 32) for (j=0; j < BUFFERSIZE; j++) x[j]  = rng(seed);
      else       for (j=0; j < 2 * BUFFERSIZE; j++) x2[j] = rng2(seed);
      fwrite (&x, sizeof (uint64_t), BUFFERSIZE, stdout);
    }
  } else { // argv != 2 : timing test
    double t = 0.0;
    sscanf (argv[2], " %lu ", &ntries);
    for (j = 0; j < 10; j++) {
      if (algo < 32) for (i=0; i < ntries; i++) rng(seed); // throws away actual random numbers
      else           for (i=0; i < ntries; i++) rng2(seed); // 32 bits
      t = crpx_update_elapsed_time_128bits (cglob->elapsed_time);
      fprintf (stderr, "%s x %lu :  %lf seconds = %.1lf million numbers/second\n", algoname, ntries, t, ((double)ntries)/(t*1.0e6));
    }
  }
  crpx_global_finalise (cglob);
  return TEST_SKIPPED;
}

