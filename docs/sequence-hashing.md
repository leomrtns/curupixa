# Rolling sequence hashes
The polynomial byte iterator in `rolling_hash.h` remains available, with unchanged hash values (26,740 windows were compared before and after the shared seed-helper change). The current implementation uses a seeded symbol table and a polynomial recurrence modulo $2^{64}$. Its first window takes linear time and each subsequent window constant time. The MMseqs2/linclust reference in the bibliography describes earlier work; the current source contains no MMseqs2 implementation. Its seed-table setup now calls the existing integer mixer instead of maintaining another copy of it.
## Split-rotation alternative
`sequence_hash.h` adds `crpx_sequence_hash_*`, using the 33/31-bit split rotation of [ntHash2 (Kazemi et al., 2022)](https://doi.org/10.1093/bioinformatics/btac564). Rotating the two subwords separately increases the rotation period from 64 to 1023. The initial window costs $O(k)$ and later valid windows $O(1)$, with fixed-size tables rather than a window-sized buffer. These remain finite 64-bit hashes with possible collisions; applications requiring exact identities should compare the underlying words when hashes match.
Choose `CRPX_HASH_DNA`, `CRPX_HASH_PROTEIN` or `CRPX_HASH_BYTES` at construction. DNA accepts A/C/G/T/U in either case, treats U as T and breaks a window at any other byte. It updates forward and reverse-complement hashes together; their sum is invariant to strand orientation. The returned orientation flag records which strand hash is smaller, with ties assigned to forward. DNA windows can exceed 32 bases, whereas `kmerhash.h` retains its exact two-bit encoding for $k\leq32$.
Protein mode maps ASCII letters to 26 case-insensitive codes. Other bytes, including gaps and stops, break the window. Byte mode accepts all 256 values literally, including zero. `crpx_sequence_hash_set_map` can define a reduced alphabet: equal codes share a symbol-table entry, and code 256 breaks the window. This allows protein grouping without embedding a particular biological scoring scheme. DNA maps retain the four complementary codes `0..3`; the complement is `3-code`.
Like the existing rolling and k-mer iterators, the new iterator borrows the input. Keep it unchanged until reset or destruction and give each worker its own iterator. `next` returns `CRPX_ERROR_NOT_FOUND` at exhaustion and leaves the supplied result unchanged. All positions refer to the original byte sequence, including positions after an invalid character.
```c
#include <curupixa/sequence_hash.h>

crpx_sequence_hash *iterator = NULL;
crpx_status status = crpx_sequence_hash_new (51, 0, CRPX_HASH_DNA, &iterator);
if (status == CRPX_OK) status = crpx_sequence_hash_reset (iterator, sequence, sequence_length);
crpx_sequence_hash_value word;
while (status == CRPX_OK && (status = crpx_sequence_hash_next (iterator, &word)) == CRPX_OK) {
  /* word.position and word.hash identify this contiguous window. */
}
crpx_sequence_hash_free (iterator);
/* NOT_FOUND is the ordinary end of iteration; handle other statuses in your application. */
```
With DNA seed zero, the new iterator uses ntHash 2.4.0's original base constants and canonical sums. Other seeds and alphabets use Curupixa's reproducible SplitMix-derived symbol tables; those outputs have their own format. Protein mode uses the generic split-rotation recurrence and a supplied alphabet map; it does not claim aaHash's particular hash values or reduced-alphabet tables. [aaHash (Wong et al., 2023)](https://doi.org/10.1093/bioadv/vbad162) is a related protein-specific method.
## Provenance and checks
The one-step rotations and DNA constants follow the MIT-licensed [ntHash source](https://github.com/bcgsc/ntHash), commit `c26bd4572a19de81e30d55042dbd33c1fd21d4b6`. Its notice is retained in `LICENSES/ntHash-MIT.txt` and installed with the library. The arbitrary-distance rotation, iterator, maps and protein/byte interfaces were written for Curupixa. No MMseqs2 source or GPL rolling-hash implementation was imported.
Tests compare incremental updates with freshly initialised windows for all three alphabets, reverse-complement invariance, invalid-character gaps and reduced alphabets. Window lengths include 1, 2, 31–33, 64–65 and 1023–1024. A separate local comparison compiled the original ntHash implementation and checked 75,591 DNA hash values for tested lengths 3..2050. The ntHash API crashed at length one in that comparison; lengths one and two were checked from the recurrence and Curupixa's direct-window tests instead. This limitation of the external check does not restrict the Curupixa constructor.
