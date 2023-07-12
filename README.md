# mosim
Molecular Orbital Solver

## Dependencies
Currently, only linux is known to be supported.

### Compiler
On linux: `clang`/`gcc` should work.

### Libraries
- `fmtlib`
- `eigen3`
- `glfw3`
- `glm`
- `glew`
- `openmp`
- `libarchive`

## Build 
```sh
meson setup build --buildtype release 
cd build 
ninja -j12
```

## Known Issues
Some examples have half-filled orbitals, which is bad since i'm not sure if the implementation of hartree-fock supports it
