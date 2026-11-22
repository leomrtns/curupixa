# Integer mixers and small-state generators
Integer mixers permute the bits of a fixed-width word. Multiplication by an odd number and XOR with a shifted copy are individually invertible modulo $2^w$. Composing these operations can disperse changes in the input across the output. A mixer used after incrementing a counter gives a small-state pseudo-random generator; this construction does not establish cryptographic security.
Curupixa retains SplitMix/Stafford Mix13, Murmur finalisers, Wellons/hash-prospector permutations and Pelican among its integer transformations. Their public-domain, Unlicense or Apache grants are recorded in `LICENSES/` and [the RNG/hash guide](random-hashing.md). The default worker generator remains xoshiro256++; no default generator changed in the MIT transition.
## Variants removed pending permission
The following implementations have been removed because their original permission could not be established from the cited publications. No third-party GPL derivation was established for these snippets; absence of an explicit grant was sufficient to exclude them from the permissive release.
| Variant | Algorithmic structure and original description |
| --- | --- |
| RRMixer, published as rrmxmx | Combine the input with two rotated copies, followed by two multiplication stages and intervening right-XOR shifts. [Evensen, July 2018](https://mostlymangling.blogspot.com/2018/07/). |
| NASAM | Combine two rotated input copies, then alternate multiplication with a combination of two right-XOR shifts. Curupixa formerly added a fixed XOR offset before that pipeline. [Evensen, January 2020](https://mostlymangling.blogspot.com/2020/01/nasam-not-another-strange-acronym-mixer.html). |
| Moremur | Alternate two odd multiplications with right-XOR shifts, with a final right-XOR shift. [Evensen, December 2019](https://mostlymangling.blogspot.com/2019/12/stronger-better-morer-moremur-better.html). |
| Curupixa's `entropy` variant | A historical offset followed by a Stafford multiply/XOR-shift variant distinct from retained SplitMix/Mix13. [Stafford, September 2011](https://zimbry.blogspot.com/2011/09/better-bit-mixing-improving-on.html). |
The following pseudocode describes the class of transformations, without reproducing any removed implementation or its exact parameter schedule:
```text
mix_word(word, stages, word_width):
  for stage in stages:
    if stage is ROTATED_XOR:
      word = word XOR rotate_right(word, stage.rotation_a)
                  XOR rotate_right(word, stage.rotation_b)
    if stage is ODD_MULTIPLICATION:
      word = (word * stage.odd_multiplier) modulo 2^word_width
    if stage is SHIFTED_XOR:
      previous = word
      word = previous XOR (previous shifted right by stage.shift_a)
      if stage has a second shift:
        word = word XOR (previous shifted right by stage.shift_b)
  return word

counter_generator(state, increment, permitted_mixer):
  state = (state + increment) modulo 2^word_width
  return permitted_mixer(state)
```
A future implementation could use the original parameter schedules if an adequate original grant is established. That step requires recording the author, applicable permission and exact source version, followed by output and state-consumption tests. The pseudocode above is a description of mathematical operations, not a claim that translating retained source would remove its licence obligations.
## Removed interfaces and identifiers
Removed integer functions are `crpx_hashint_rrmixer64`, `crpx_hashint_nasam64`, `crpx_hashint_moremur64` and `crpx_hashint_entropy`. Removed generator functions are `crpx_rng_rrmixer_seed64` and `crpx_rng_moremur_seed64`. Context generator IDs 2 (Moremur) and 13 (RRMixer) now return `CRPX_ERROR_INVALID_ARGUMENT`; previously retired ID 19 remains invalid. No replacement algorithm silently takes their names or IDs. Existing clients using these functions must select an explicitly named retained alternative and rebuild.
