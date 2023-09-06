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

The code is run from the command line and a runcard must always be given as the
first argument:

`./celadro runcard.dat [options]`

A runcard is a simple file providing the parameters for the run. Example
runcards can be found in the `simulation_examples/` directory. Every option can also be
given to the program using the command line as `./celadro runcard.dat --option=arg`.
A complete list of available options can be obtained by typing `./celadro -h`.

By default the program writes output files in the current directory. This can be
changed using `--output=dir/` or `-o dir/`, where `dir/` is the target
directory. The program also supports compressed (using zip) output with the option
flag `--compression` or `-c`. Typical usage is `./celadro runcard.dat -fco output/`

Type `./celadro -h` for a list of available options.

## Examples

Examples runs and ploting scripts can be found in the `simulation_examples` directory. 

## Visualization

VTK library is needed for visualization. Go to `VTK_VolRender/vtk_VolRender_01012021.cpp` to change camera position. The code can be compile with: 

```
cmake .
make 
```
You also need to create a .vtk file from .json. This can be done with the following python code: 

```
/scripts/write_vtk_from_JSON_03012021.py 
```

you can run it by simply change what frame(s) to use for creating .vtk file(s). Then simply write: 

```
python write_vtk_from_JSON_03012021.py [output directory containing .json files]
```

