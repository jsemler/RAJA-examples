# RAJA Proxy Applications

This repository contains three proxy applications using the RAJA programming
framework. The CMakeLists file is set up to fetch the latest version of RAJA
from GitHub and build it along with the applcations. All you need to do is
ensure CMake is using a compiler that supports C++11 and OpenMP:

    mkdir build && cd build
    cmake ../
    make

If you have a version of RAJA built already and you'd like to use that, just set
`RAJA_DIR` to point to it before running CMake.

## Testing RAJA PRs

To have CMake download and build a branch other than `develop`, you can
specify the branch name in the `RAJA_GIT_TAG` variable. Likewise, to change the
repository URL, simply set the variable `RAJA_GIT_REPO`.
