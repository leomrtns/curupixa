#include <curupixa.h> 

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

/* Testing simple hash functions with dieharder, command line below:
 * `./tests/dieharder_hash 7 | dieharder -g 200 -a`
 * notice that according to dieharder's man page, the _generator_ number 200 is stdin (thus "-g 200", not "-f 200")
 *
 * WARNING: this program will run for a _very_ long time 
 */
int main(int argc, char **argv)
{
  uint64_t i, x64, seed = 0x3581cf2a5687e23ULL;
  uint32_t x32;
  uint8_t algo;

  if (argc == 1) return TEST_SKIPPED;
  crpx_global_t cglob = crpx_global_init (0,0,"debug");

  sscanf (argv[1], " %hhd ", &algo);

  if (algo ==0) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_fastmix64 (i); // NG
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else if (algo == 1) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_murmurmix64 (i); // OK 
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else if (algo == 2) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_fnv_hash64 ((void*)(&i), sizeof (uint64_t)); // NG
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else if (algo == 3) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 = crpx_hash_fletcher32 ((void*)(&i), sizeof (uint64_t)); // NG
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 4) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 = crpx_hash_jenkins ((void*)(&i), sizeof (uint64_t)); // OK 
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 5) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 = crpx_hash_jenkins_mailund_seed32 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // NG
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 6) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 = crpx_hash_mailund_seed32 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // OK
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 7) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 = crpx_hash_rotating_seed32 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // NG
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 8) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_fasthash64_seed64 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // OK +- // passed practrand 32bits
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else if (algo == 9) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 = crpx_fnv_hash32 ((void*)(&i), sizeof (uint64_t)); // NG
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 10) {
    for (i=0; i < UINT64_MAX; i++) {
      x32 =crpx_hsieh_hash32_seed32 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // NG 
      fwrite (&x32, sizeof (uint32_t), 1, stdout);
    }
  } else if (algo == 11) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_metrohash64_v1_seed64 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // OK
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else if (algo == 12) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_rrmixer64 (i); // OK // passed practrand 32bits
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else if (algo == 13) {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_moremur64 (i); // OK // suspicious/unusual practrand 32bits
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  } else {
    for (i=0; i < UINT64_MAX; i++) {
      x64 = crpx_metrohash64_v2_seed64 ((void*)(&i), sizeof (uint64_t), (void*)(&seed)); // OK
      fwrite (&x64, sizeof (uint64_t), 1, stdout);
    }
  }

  crpx_global_finalise (cglob);
  return TEST_SKIPPED;
}

