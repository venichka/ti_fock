# ti_fock

Computational tools for time-interface quantum optics with Fock states.

## Repository structure

- `src/TI_quantum.jl`: core package code.
- `test/runtests.jl`: automated smoke/regression tests.
- `notebooks/`: Pluto notebooks for experiments/figures.
- `.github/workflows/ci.yml`: GitHub Actions CI for package checks.

## Package setup (portable)

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Output folders (`figs/`, `data/`)

`TI_quantum.path()` auto-creates portable output directories.

```julia
using TI_quantum
PATH_FIGS, PATH_DATA = path()
```

By default they are created in the repository root:

- `./figs/`
- `./data/`

Optional override:

```bash
export TI_QUANTUM_OUTPUT_DIR=/absolute/output/path
```

## Notebooks

Notebook dependencies are isolated in `notebooks/Project.toml`.
Notebook environment targets Julia `1.11+` (validated on Julia `1.12.4` on March 2, 2026).

If you use `juliaup`, install/select that channel first:

```bash
juliaup add release
julia +release --project=notebooks -e 'using Pkg; Pkg.instantiate()'
```

```bash
julia --project=notebooks -e 'using Pkg; Pkg.instantiate()'
```

The notebooks load the local package source directly from `../src/TI_quantum.jl`.
