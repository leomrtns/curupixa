#include <curupixa.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

/* Testing available PRNGs with dieharder and practrand, command lines below:
 * `./tests/dieharder_rng 7 | dieharder -g 200 -a` 
 * `./tests/dieharder_hash 7 | RNG_test stdin32;
 * notice that according to dieharder's man page, the _generator_ number 200 is stdin (thus "-g 200", not "-f 200")
 *
 * WARNING: this program will run for a _very_ long time 
 */
#define BUFFERSIZE 512

int main(int argc, char **argv)
{
  uint64_t i, ntries;
  uint8_t algo;
  int j;

  if (argc == 1) return TEST_SKIPPED;
  crpx_global_t cglob = crpx_global_init (0, "debug");

  sscanf (argv[1], " %hhu ", &algo);
  if (algo) crpx_set_random_generator (cglob, (uint8_t) (algo & 255), 0);

  for (i=0; i < (cglob->rng_size * cglob->nthreads); i++) {
    fprintf (stderr, "%17lx ", cglob->rng_seed_vector[i]); if (!((i+1)%4)) fprintf (stderr, "\n");
  }
  fprintf (stderr, "%lf seconds to set seed vector\n", crpx_update_elapsed_time_128bits (cglob->elapsed_time));

  if (argc == 2) { // dieharder test
    uint64_t x[BUFFERSIZE];
    for (i=0; i < UINT64_MAX; i++) {
      for (j=0; j < BUFFERSIZE; j++) x[j] = crpx_random_64bits (cglob);
      fwrite (&x, sizeof (uint64_t), BUFFERSIZE, stdout);
    }
  } else { // argv != 2 : timing test
    double t = 0.0;
    sscanf (argv[2], " %lu ", &ntries);
    for (j = 0; j < 10; j++) {
      for (i=0; i < ntries; i++) crpx_random_64bits (cglob); // throws away actual random numbers
      t = crpx_update_elapsed_time_128bits (cglob->elapsed_time);
      printf ("%32s x %lu :  %lf seconds = %.1lf million numbers/second\n", cglob->rng_name, ntries, t, ((double)ntries)/(t*1.0e6));
    }
  }
  crpx_global_finalise (cglob);
  return TEST_SKIPPED;
}

