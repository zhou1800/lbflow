# Source History And Baseline

The restored solver has two distinct source revisions in this repository:

- `66923961c911684cff8b703914118e431f8ce445` imports the recovered DRLab
  0.2.2.1 numerical implementation before the repairs. This is a 2026 import,
  not an assertion that the code was published in 2017. The original archive's
  revision is `1bef9bbe7047ca74280f5870f4787be4df12fe20`.
- The following implementation commit integrates the numerical repairs,
  supported root build, regression checks and publication sources.

The import excludes IDE state, editor undo files, large raw sample geometry
and the generated configuration header. The original source hashes remain in
`../source-manifest.json`. Original numerical source bytes are preserved in the
import; later changes can be inspected with `git diff 6692396 -- solver`.
The incomplete C prototype from the prior LBFlow tree remains in
`../../../legacy/c-prototype/`.

For an isolated original-solver counterexample, run from the repository root:

```sh
mkdir -p results/original-source
git archive 66923961c911684cff8b703914118e431f8ce445 solver | tar -x -C results/original-source
cmake -S results/original-source/solver -B results/original-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -include cstdint"
cmake --build results/original-build --target LBSim --parallel 4
python3 scripts/benchmarks.py verification --solver results/original-build/apps/LBSim --root results/original-verification
```

The historical CMake file overwrites `CMAKE_CXX_FLAGS`, so the compatibility
include is supplied in `CMAKE_CXX_FLAGS_RELEASE`. These commands were checked
with CMake 3.28.3. The forced `<cstdint>` include is a compiler compatibility adjustment; it does
not change numerical operators. The original executable is expected to fail
the analytical force checks. The maintained solver passes them. Old checkpoint
files do not have the repaired stage contract, so use separate result roots.
