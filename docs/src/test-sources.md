# Test source map

This test suite is organized by numerical behavior type. The goal is to make
it easy to trace each group of checks to mathematical or software references.

## File layout

- `test_elementary_identities.jl`
  - Reductions to elementary/special functions (`exp`, `sin`, `cos`, `besselk`)
  - Sources:
    - DLMF §16.17 (definition and standard identities context)
    - mpmath Meijer G docs examples for trigonometric and Bessel reductions
- `test_api_and_symmetry.jl`
  - API-level consistency checks: split vs index form, vector vs tuple inputs,
    inversion symmetry
  - Sources:
    - Standard Meijer G inversion identity
    - MATLAB symbolic `meijerG` algorithm notes on transformed-parameter identity
- `test_domain_and_reduction.jl`
  - Input-contract errors (`m`, `n`, `z=0`), simple-pole guards, order-reduction
    equivalence, and near-`|z|=1` branch-switch smoke tests
  - Sources:
    - DLMF §16.17 contour/series regime split
    - Implementation contract in `src/core.jl`
- `test_high_precision_references.jl`
  - BigFloat checks and fixed numerical references
  - Sources:
    - mpmath `meijerg` docs and `mpmath/tests/test_functions2.py`
    - Internal high-precision self-convergence reference (1024-bit)

## External references used for fixed values

1. mpmath documentation, Hypergeometric functions → Meijer G function:
   https://mpmath.org/doc/current/functions/hypergeometric.html#meijer-g-function
2. mpmath regression tests (`test_meijerg`):
   https://github.com/mpmath/mpmath/blob/main/mpmath/tests/test_functions2.py
3. DLMF §16.17 Definition of Meijer G:
   https://dlmf.nist.gov/16.17
4. MATLAB symbolic `meijerG` algorithm page:
   https://www.mathworks.com/help/symbolic/sym.meijerg.html

## Notes on scope

- These tests intentionally avoid confluent-pole logarithmic cases except where
  rejection behavior is being tested.
