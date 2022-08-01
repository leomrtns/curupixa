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
int main(int argc, char **argv)
{
  uint64_t i;
  uint64_t x, seed[16];
  uint64_t (*rng)(void*);
  uint8_t algo;

  if (argc == 1) return TEST_SKIPPED;
  crpx_global_t cglob = crpx_global_init (0,0,"debug");
  crpx_generate_bytesized_random_seeds (cglob, seed, 128);

  sscanf (argv[1], " %hhd ", &algo);
	
  switch (algo) { 
    case 0: rng = &crpx_rng_romu_seed256; break;
    case 1: rng = &crpx_rng_romu_seed192; break;
    case 2: rng = &crpx_rng_romu_seed128; break;
    case 3: rng = &crpx_xoro128plus_seed128; break;
    case 4: rng = &crpx_xs64star_seed64; break;
    case 5: rng = &crpx_rng_wyhash64_state64; break;
    case 6: rng = &crpx_rng_splitmix64_seed64; break;
    case 7: rng = &crpx_rng_rrmixer_seed64; break;
    case 8: rng = &crpx_rng_moremur_seed64; break;
    default: rng = &crpx_xs128plus_seed128; break;
  }
  for (i=0; i < UINT64_MAX; i++) {
    x = rng(seed);
    fwrite (&x, sizeof (uint64_t), 1, stdout);
  };
  crpx_global_finalise (cglob);
  return TEST_SKIPPED;
}

