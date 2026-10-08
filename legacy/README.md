# Historical C Prototype

`c-prototype/` preserves the unfinished LBFlow C rewrite that preceded the
restored C++ simulator. Its `main.c` is a Hello World entry point; it does not
perform fluid simulation. Its JSON examples, GoogleTest download logic, GitLab
runner setup and developer scripts belong to that historical implementation.

These files are retained for source history and possible future porting. They
are excluded from the supported root CMake build and CI. Use the repository-root
README and generated `.prm` examples for current simulations. Third-party source
license notices are retained in place.
