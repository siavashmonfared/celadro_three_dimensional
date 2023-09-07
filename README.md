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

The code is run from the command line and a simCard must always be given as the
first argument, from directory `example/`:

`../build/celadro simCard.dat [options]`

A simCard.dat is a simple file providing the parameters for the run. Example
simCard can be found in the `example` directory. Every option can also be
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


