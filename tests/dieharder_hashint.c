#include <curupixa.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

/* Testing simple hash functions with dieharder and PractRand, command lines below:
 * `./tests/dieharder_hash 7 | dieharder -g 200 -a`
 * `./tests/dieharder_hash 7 | RNG_test stdin32;
 * notice that according to dieharder's man page, the _generator_ number 200 is stdin (thus "-g 200", not "-f 200")
 * WARNING: this program will run for a _very_ long time 
 */

#define BUFFERSIZE 512

int main(int argc, char **argv)
{
  uint8_t algo;
  uint64_t i, k, ntries;
  char algoname[16] = {'\0'};
  int j;
  uint64_t (*h)(uint64_t);
  uint32_t (*h2)(uint32_t);

  if (argc == 1) return TEST_SKIPPED;
  crpx_global_t cglob = crpx_global_init (0,0,"debug");

  sscanf (argv[1], " %hhd ", &algo);

  if (algo < 8) {
    switch (algo) { //fastmix64, murmurmix64 fail
      case 0: h = &crpx_hashint_splitmix64;   strcpy (algoname, "splitmix64"); break; // fails a few dieharder // 2^32 practrand
      case 1: h = &crpx_hashint_degski64;     strcpy (algoname, "degski64  "); break; // fails some dieharder // fail practrand
      case 2: h = &crpx_hashint_nasam64;      strcpy (algoname, "nasam64   "); break; // great dieharder // pass practrand >2^36
      case 3: h = &crpx_hashint_pelican64;    strcpy (algoname, "pelican64 "); break; // great dieharder // pass practrand >2^36  
      case 4: h = &crpx_hashint_rrmixer64;    strcpy (algoname, "rrmixer64 "); break; // ok    dieharder // pass practrand >2^35
      case 5: h = &crpx_hashint_moremur64;    strcpy (algoname, "moremur64 "); break; // great dieharder // 2^32 practrand
      case 6: h = &crpx_hashint_staffordmix64;strcpy (algoname, "stafford64"); break; // fails a few dieharder // 2^32 practrand
      case 7: h = &crpx_hashint_entropy;      strcpy (algoname, "entropy   "); break; // fails a few dieharder // 2^34 practrand
    }
  } else { // 32 bits
    switch (algo) {
      case 8:  h2 = &crpx_hashint_jenkins;     strcpy (algoname, "jenkins   "); break; // fails dieharder // fail practrand
      case 9:  h2 = &crpx_hashint_jenkins_v2;  strcpy (algoname, "jenkins_v2"); break; // fails dieharder // fail practrand
      case 10: h2 = &crpx_hashint_avalanche;   strcpy (algoname, "avalanche "); break; // fails dieharder // fail practrand
      case 11: h2 = &crpx_hashint_murmurmix;   strcpy (algoname, "murmurmix "); break; // fails a few dieharder // fail practrand
      case 12: h2 = &crpx_hashint_wellons3ple; strcpy (algoname, "wellons3pl"); break; // great dieharder // 2^32 practrand
      case 13: h2 = &crpx_hashint_degski;      strcpy (algoname, "degski    "); break; // fails some dieharder // fail practrand
      default: h2 = &crpx_hashint_wellons;     strcpy (algoname, "wellons   "); break; // weak dieharder // 2^32 practrand
    }
  } 
  fprintf (stderr, "%s :  %lf seconds to initialise (and start timer)\n", algoname, crpx_update_elapsed_time_128bits (cglob->elapsed_time));

  if (argc == 2) { // dieharder/practrand test
    uint64_t x[BUFFERSIZE];
    uint32_t *x2 = (uint32_t*) x;
    for (i=0; i < UINT64_MAX; i++) {
      if (algo < 8) for (j=0; j < BUFFERSIZE; j++) x[j]  = h(k++);
      else for (j=0; j < 2 * BUFFERSIZE; j++) x2[j] = h2((uint32_t) k++);
      fwrite (&x, sizeof (uint64_t), BUFFERSIZE, stdout);
      /*DEBUG*/ // for (j=0; j < BUFFERSIZE; j++) {printf ("%16lx ", x[j]); if ((j+1) % 8 == 0) printf ("\n");}
    }
  } else { // timing
    double t = 0.0;
    sscanf (argv[2], " %lu ", &ntries);
    k = 1234;
    for (j = 0; j < 10; j++) {
      if (algo < 8) for (i=0; i < ntries; i++) h(k++); // throws away actual values; k grows indefinitely
      else for (i=0; i < ntries; i++) h2((uint32_t)(k++));
      t = crpx_update_elapsed_time_128bits (cglob->elapsed_time);
      fprintf (stderr, "%s %2d x %lu :  %lf seconds = %.1lf million numbers/second\n", algoname, j, ntries, t, ((double)ntries)/(t*1.0e6)); 
    }
  }
  crpx_global_finalise (cglob);
  return TEST_SKIPPED;
}
