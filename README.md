# DNAOAD

Dual Numbers for Arbitrary Order Automatic Differentiation

This library performs automatic forward differentiation of code, by
setting the independent variable to be of a specific type and
initializing it's value and first derivative, then performing standard
operations on it.

## Compilation

The code uses *cmake* and *make* to compile. Simple compilation, on
unix-like systems, from the top source directory compile in the
selected *_build* directory via

> cmake -B_build .

> cd _build

> make

For windows compilation, either use a framework like *MinGW* and the
above compilation, or *cmake* is integrated with visual studio.

If locally installed, where to find Fortuno source files and the fypp
preprocessor can be specified at configutation time with

> cmake -DFETCHCONTENT_SOURCE_DIR_FYPP=/path/to/fypp

or

> cmake -DFETCHCONTENT_SOURCE_DIR_FORTUNO=/path/to/Fortuno

where in each cases, the path to the top level of the their
repositories (or also directly to the binary location for fypp) should
be specified.

Alternativelty, these packages will be fetched as subprojects from
their respectitve repositories.

### Preprocessing

The source code uses the [fypp](https://github.com/aradi/fypp)
preprocessor for *.F90* source files. This tool is downloaded during
the compilation phase if not already installed. As a result, building
the code has a requirement for Python3.5 or above.

### Testing

In the build directory, after compilation type:

> ctest

to run unit tests for a number of functions of the code with the
current build.

The [fortuno](https://github.com/fortuno-repos) framework is used for
testing, and is also downloaded during the compilation phase if not
installed.

The github repository also runs these tests on push or PR for code
changes, testing is currently performed for

| OS | Compiler | Version |
| ----------- | ----------- | ----------- |
| Linux | gfortran | 13 |
| Linus | intel ifx | one-api latest |
| Windows | gfortran | mingw latest |


See

> .github/workflows

for the actions that are invoked.

## Examples

See the `examples` folder for demonstrations of the library. Again,
these are compiled by the cmake framework.

## Code and algorithms

For more details, see *Section 4* of the
[article](https://doi.org/10.48550/arXiv.2501.04159)

## License

DNAOAD is licensed under the MIT license. The SPDX license identifier
for this project is MIT.
