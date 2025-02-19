![example simulation](config.gif)

# Celadro_3D: Cells as active droplets

Three-dimensional phase-field modelling of epithelial cells using finite-difference integrator. 

## Building 

```
mkdir build
cd build
cmake ..
make
```

We rely on the `boost::program_options` which must be installed prior to
building the program. We also use modern C++ features, such that you will
require a modern compiler.

## Running

The code is run from the command line and a simulation card `simCard.dat` and an input file `input_str.dat` must always be provided. Specifically, the `simCard.dat` should be given as an argument 

As an exampke, from directory `example` execute:
`../build/celadro simCard.dat [options]`

The `simCard.dat` specifies simulation details such as number of time steps, boundary conditions. The `input_str.dat` contains initial positions for cells.


Examples can be found in the `example` directory. Every option can also be
given to the program using the command line as `../build/celadro simCard.dat --option=arg`.
A complete list of available options can be obtained by typing `../celadro -h`.

By default the program writes output files in the current directory. This can be
changed using `--output=dir/` or `-o dir/`, where `dir/` is the target
directory. The program also supports compressed (using zip) output with the option
flag `-compress-full`. Typical usage is `../buil/celadro -compress-full simCard.dat`

Type `../build/celadro -h` for a list of available options.

## Examples

Examples runs and ploting scripts can be found in the `example` directory. 

## Visualization

VTK library is needed for visualization. A code is provided in `example/vtk_VolRender_05012021.py`. You can view the .vtk files using paraview or by running `example/vtk_VolRender_05012021.py`.


## References

The relevant publications are: 

```
Monfared et al. eLife 2023;12:e82435. DOI: https://doi.org/10.7554/eLife.82435
Monfared et al. arXiv:2210.08112 [cond-mat.soft]. DOI: https://doi.org/10.48550/arXiv.2210.08112
Mueller et al., PRL 122, 048004. DOI: https://doi.org/10.1103/PhysRevLett.122.048004
```


