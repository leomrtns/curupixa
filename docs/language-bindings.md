# Calling Curupixa from Python and other languages
The three libraries provide C functions with the same interfaces used by C and C++ programs. A language binding can call them through `ctypes`, CFFI or a compiled extension. Curupixa does not currently provide an installable Python package.
## Build shared libraries
```sh
cmake -S . -B build/shared -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build/shared --parallel
python3 examples/python_api.py build/shared
```
The example loads all three libraries, constructs trees, maps two gene copies to one species and evaluates a conditional-clade tree probability. It declares every function signature before calling C. On Linux the files are named `libcurupixa_phylogenetics.so`, `libcurupixa_gene_species.so` and `libcurupixa_inference.so`; macOS shared-library names end in `.dylib`. The build records the libraries' dependencies so that the upper libraries can find the lower ones. For installed packages, the system loader must also be able to find their installation directory.
## Declare types and check errors
With `ctypes`, use `c_void_p` for a pointer to an opaque object, `c_size_t` for `size_t`, and `c_uint32` for a node ID. A pointer returned through an output argument needs `byref()`. Set `argtypes` and `restype` before calling a function; Python's default integer return type cannot safely represent a pointer.
Most calculation functions return `CRPX_OK` (zero) on success. Otherwise, `crpx_status_string()` gives an error description that your wrapper can turn into an exception. Some numerical results can legitimately be negative infinity, for example the log probability of an unsupported tree in the conditional-clade model. Follow the individual function's documented contract rather than interpreting every non-finite value as a failed call.
## Release memory and retain inputs
An object created by `crpx_tree_new()` or `crpx_tree_from_newick()` must eventually be released with `crpx_tree_free()`. Use the corresponding destructor for each other object type. Do not release Curupixa memory through Python's allocator. Put destruction in `finally`, or write a wrapper class with an explicit `close()` method and context-manager support.
Some accessors return a view of arrays that belong to a workspace. Copy these arrays into Python if they must survive another evaluation or destruction of that workspace. Keep Python arrays and callbacks alive for as long as C may use them. In particular, an amalgamation sampler's callback must remain alive throughout the sampling call. The individual headers state whether an input is copied or retained.
## Repeat calculations and use threads
A workspace is an object storing temporary arrays for a calculation. Create it once, reuse it for successive evaluations and destroy it when finished. Prepared indices can similarly avoid repeating tree setup. These are separate from the input tree, so constructing a tree alone does not allocate reconciliation or likelihood arrays.
Use a different mutable workspace for each simultaneous calculation. Share a read-only model or index only where its interface allows this. Python's `ctypes.CDLL` calls can release the interpreter lock, so the Python interpreter does not protect a shared C workspace from concurrent access. Curupixa does not start worker threads itself.
