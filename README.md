## Invoking FDFB and Homomorphic Decomposition
### Building
First build and install OpenFHE using the following commands. Replace `YOUR_CMAKE_INSTALL_PREFIX`  with your own path.
```shell
git clone https://github.com/openfheorg/openfhe-development.git
cd openfhe-development
git checkout 745a492
git apply 20230504.patch
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DWITH_NATIVEOPT=ON -DCMAKE_INSTALL_PREFIX=YOUR_CMAKE_INSTALL_PREFIX
cmake --build build -j8 && cmake --install build
```
[Google Benchmark](https://github.com/google/benchmark.git) is also required if you want to meaure the performance of FDFB and homomorphic decomposition algorithms.

Modify the `CMakeLists.txt` in the current directory (not OpenFHE directory) by replacing `YOUR_CMAKE_INSTALL_PATH` with the path your just passed to CMake. If Google Benchmark is unavailable, comment the corresponding lines in `CMakeLists.txt`  and remove `benchfdfb` and `benchsign` from target executables.

Now run the following commands.
```shell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```
The executables can be found in the `build` directory.
### Running
`benchfdfb` benchmarks all FDFB algorithms. `benchsign` benchmarks all homomorphic decomposition algorithms. Their CLI is managed by Google Benchmark.

`evalfunc`  evaluates FDFB. `evalsign` performs homomorphic decomposition. `evalrelu` evaluates large-precision ReLU with $\textbf{HomDecomp-Reduce}$. The parameter set used by the three executables can be specified using command line arguements. A basic usage is printed when invoking these executables without any argument. Refer to their source code for a full list of predefined parameters sets.

## Parameter Selection and Noise Analysis
`param.py` provides convenient functions to estimate the noise of FDFB and homomorphic decomposition algorithms. To use these functions, type `from param import *` in Python command line. Refer to the python file for details on the usage of these functions.

`sage_code.txt` provides a simple template to measure the concrete security of an LWE (or RLWE) scheme. Both [SageMath](https://www.sagemath.org/) and [lattice-estimator](https://github.com/malb/lattice-estimator.git) are required to use it.