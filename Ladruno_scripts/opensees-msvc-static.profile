[settings]
os=Windows
arch=x86_64
build_type=Release
compiler=msvc
# Ladruno: MSVC toolset version is machine-dependent. Default 194 (VS 2022 /
# v17, the canonical build box); override to 195 for VS 2026 / v18, etc. via
# the LADRUNO_MSVC_VERSION env var so the committed profile stays portable.
compiler.version={{ os.getenv("LADRUNO_MSVC_VERSION", "194") }}
compiler.cppstd=17
compiler.runtime=static
compiler.runtime_type=Release

[conf]
tools.cmake.cmaketoolchain:generator=Ninja
