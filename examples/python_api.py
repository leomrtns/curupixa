"""Call each Curupixa shared library using ctypes; run with a CMake build directory."""
import ctypes as c
from pathlib import Path
import sys


def main():
  directory = Path(sys.argv[1] if len(sys.argv) > 1 else "build/shared").resolve()

  def load(component):
    for suffix in ("so", "dylib", "dll"):
      path = directory / f"libcurupixa_{component}.{suffix}"
      if path.exists():
        return c.CDLL(str(path), mode=c.RTLD_GLOBAL)
    raise FileNotFoundError(f"Shared library for {component} not found in {directory}")

  phylo = load("phylogenetics")
  mapped = load("gene_species")
  inference = load("inference")
  phylo.crpx_status_string.argtypes = [c.c_int]
  phylo.crpx_status_string.restype = c.c_char_p

  def check(status):
    if status:
      raise RuntimeError(phylo.crpx_status_string(status).decode())

  phylo.crpx_tree_from_newick.argtypes = [c.c_char_p, c.POINTER(c.c_void_p), c.POINTER(c.c_size_t)]
  phylo.crpx_tree_from_newick.restype = c.c_int
  phylo.crpx_tree_free.argtypes = [c.c_void_p]
  phylo.crpx_tree_free.restype = None
  phylo.crpx_tree_node_count.argtypes = [c.c_void_p]
  phylo.crpx_tree_node_count.restype = c.c_size_t
  mapped.crpx_gene_species_map_labels.argtypes = [c.c_void_p, c.c_void_p, c.c_char, c.c_size_t,
                                                 c.POINTER(c.c_uint32), c.c_size_t]
  mapped.crpx_gene_species_map_labels.restype = c.c_int
  inference.crpx_metropolis_accept.argtypes = [c.c_double, c.c_double, c.c_double, c.POINTER(c.c_bool)]
  inference.crpx_metropolis_accept.restype = c.c_int

  phylo.crpx_tree_label.argtypes = [c.c_void_p, c.c_uint32]
  phylo.crpx_tree_label.restype = c.c_char_p
  inference.crpx_amalgamation_new.argtypes = [c.POINTER(c.c_void_p), c.POINTER(c.POINTER(c.c_uint32)),
                                             c.c_size_t, c.c_size_t, c.POINTER(c.c_char_p),
                                             c.POINTER(c.c_uint64), c.POINTER(c.c_void_p)]
  inference.crpx_amalgamation_new.restype = c.c_int
  inference.crpx_amalgamation_workspace_new.argtypes = [c.POINTER(c.c_void_p)]
  inference.crpx_amalgamation_workspace_new.restype = c.c_int
  inference.crpx_amalgamation_log_probability.argtypes = [c.c_void_p, c.c_void_p, c.POINTER(c.c_uint32),
                                                         c.c_void_p, c.POINTER(c.c_double)]
  inference.crpx_amalgamation_log_probability.restype = c.c_int
  for name in ("crpx_amalgamation_free", "crpx_amalgamation_workspace_free"):
    function = getattr(inference, name)
    function.argtypes = [c.c_void_p]
    function.restype = None

  gene, species = c.c_void_p(), c.c_void_p()
  model, workspace = c.c_void_p(), c.c_void_p()
  try:
    check(phylo.crpx_tree_from_newick(b"((A|1,A|2),B|1);", c.byref(gene), None))
    check(phylo.crpx_tree_from_newick(b"(A,B);", c.byref(species), None))
    nodes = phylo.crpx_tree_node_count(gene)
    mapping = (c.c_uint32 * nodes)()
    check(mapped.crpx_gene_species_map_labels(gene, species, b"|", 0, mapping, nodes))
    accepted = c.c_bool()
    check(inference.crpx_metropolis_accept(0.0, 0.0, -0.5, c.byref(accepted)))
    assert accepted.value
    # One sampled topology has conditional-clade probability one.
    labels = (c.c_char_p * 3)(b"A|1", b"A|2", b"B|1")
    tip_ids = {phylo.crpx_tree_label(gene, node): node for node in range(nodes)}
    tips = (c.c_uint32 * 3)(*(tip_ids[label] for label in labels))
    trees = (c.c_void_p * 1)(gene.value)
    tip_arrays = (c.POINTER(c.c_uint32) * 1)(tips)
    check(inference.crpx_amalgamation_new(trees, tip_arrays, 1, 3, labels, None, c.byref(model)))
    check(inference.crpx_amalgamation_workspace_new(c.byref(workspace)))
    log_probability = c.c_double()
    check(inference.crpx_amalgamation_log_probability(model, gene, tips, workspace, c.byref(log_probability)))
    assert log_probability.value == 0.0
    print(f"Mapped gene tree with {nodes} nodes; all three C libraries called successfully.")
  finally:
    inference.crpx_amalgamation_workspace_free(workspace)
    inference.crpx_amalgamation_free(model)
    phylo.crpx_tree_free(gene)
    phylo.crpx_tree_free(species)


if __name__ == "__main__":
  main()
