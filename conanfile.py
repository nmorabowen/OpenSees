from conan import ConanFile
from conan.tools.cmake import cmake_layout

class OpenSeesDependencies(ConanFile):
    name            = "OpenSeesDependencies"
    version         = "1.0.0"
    description     = "Provides Software Packages needed to build OpenSees"
    license         = "BSD-3-Clause"
    author          = "Ladruño"

    settings        = "os", "compiler", "build_type", "arch"
    options         = {"shared": [True, False]}
    default_options = {
        "shared": False,
        "mkl-static/*:threaded": False,
        "ipp-static/*:simcenter_backend": True
    }

    generators      = "CMakeDeps", "CMakeToolchain"

    def requirements(self):
        self.requires("hdf5/1.14.0")
        self.requires("tcl/8.6.11")
        self.requires("zlib/1.2.13")
        self.requires("eigen/3.4.0")

    def layout(self):
        cmake_layout(self)
