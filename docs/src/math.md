# Mathematical Notes

`MeijerG.jl` evaluates the Meijer G-function via Slater-style residue expansions. See [DLMF §16.17](https://dlmf.nist.gov/16.17) for the Mellin–Barnes integral representation and residue (Slater) expansions, [Slater (1966)](https://archive.org/details/generalizedhyper0000unse) for the classical treatment, and the [Mellin–Barnes / Barnes integral](https://en.wikipedia.org/wiki/Mellin%E2%80%93Barnes_integral) overview for background on contour choices and residue calculus.

## Definition

```math
G_{p,q}^{m,n}\!\left(z\;\middle|\;\begin{matrix}a_1,\dots,a_p\\b_1,\dots,b_q\end{matrix}\right)
= \frac{1}{2\pi i}\int_{\mathcal{L}}
\frac{\prod_{j=1}^{m}\Gamma(b_j-s)\prod_{j=1}^{n}\Gamma(1-a_j+s)}
{\prod_{j=m+1}^{q}\Gamma(1-b_j+s)\prod_{j=n+1}^{p}\Gamma(a_j-s)}z^s\,ds.
```

## Expansion selection

- If `p < q`, use lower expansion.
- If `p > q`, use upper expansion.
- If `p = q`, use lower for `|z| ≤ 1` and upper for `|z| > 1`.

## Lower expansion

For simple poles among the first `m` lower parameters:

```math
G_{p,q}^{m,n}(z) = \sum_{k=1}^{m} A_k \, z^{b_k} \, {}_pF_{q-1}(\cdots; (-1)^{p-m-n}z)
```

with gamma-ratio coefficient `A_k`.

## Upper expansion

For simple poles among the first `n` upper parameters:

```math
G_{p,q}^{m,n}(z) = \sum_{h=1}^{n} B_h \, z^{a_h-1} \, {}_qF_{p-1}(\cdots; (-1)^{q-m-n}z^{-1})
```

with gamma-ratio coefficient `B_h`.

## References

- DLMF §16.17: https://dlmf.nist.gov/16.17
- mpmath Meijer G documentation: https://mpmath.org/doc/current/functions/hypergeometric.html#meijer-g-function
- Meijer G (Wikipedia): https://en.wikipedia.org/wiki/Meijer_G-function
- Mellin–Barnes / Barnes integral (Wikipedia): https://en.wikipedia.org/wiki/Mellin%E2%80%93Barnes_integral
- Slater, L.J. (1966), Generalized Hypergeometric Functions: https://archive.org/details/generalizedhyper0000unse
