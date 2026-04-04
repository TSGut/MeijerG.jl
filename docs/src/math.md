# Mathematical Notes

Where no reductions to simpler functions are known `MeijerG.jl` evaluates the Meijer G-function via residue expansions derived
from the Mellin-Barnes contour integral representation. See
[DLMF §16.17](https://dlmf.nist.gov/16.17).

## Definition

The Meijer G-function is defined by the Mellin-Barnes contour integral:

```math
G_{p,q}^{m,n}\left(z\;\middle|\;\begin{matrix}a_1,\dots,a_p\\b_1,\dots,b_q\end{matrix}\right)
= \frac{1}{2\pi i}\int_{\mathcal{L}}
\frac{\displaystyle\prod_{j=1}^{m}\Gamma(b_j-s)\;\prod_{j=1}^{n}\Gamma(1-a_j+s)}
     {\displaystyle\prod_{j=m+1}^{q}\Gamma(1-b_j+s)\;\prod_{j=n+1}^{p}\Gamma(a_j-s)}
\,z^s\,ds
```

with $0 \le m \le q$, $0 \le n \le p$, and $z \ne 0$. The contour $\mathcal{L}$
separates poles of $\Gamma(b_j - s)$ on the right from poles of
$\Gamma(1 - a_j + s)$ on the left. [DLMF §16.17](https://dlmf.nist.gov/16.17) shows a nice representative diagram of the different contour scenarios.

While the contour integral approach is implemented in this package and exported using the `meijerg_contour` API, the recommended default evaluator `meijerg` uses the Slater expansion via the residue theorem to obtain finite sums of
generalized hypergeometric functions instead because accurate contour integration is temperamental and often slow.

When poles in the selected residue family are confluent, meaning integer-separated
parameters in the active set, Slater's hypergeometric series expansion fails.
`MeijerG.jl` evaluates these cases by a parameter-perturbation limit of the
same residue representation instead. If all else has failed, the package resorts to contour integration.

## Slater's Lower Expansion

When $p < q$, or when $p = q$ and $|z| \le 1$, the contour is closed around
the poles of $\Gamma(b_j - s)$ for $j = 1, \dots, m$. This yields
([DLMF §16.17.2](https://dlmf.nist.gov/16.17.2)):

```math
\boxed{G_{p,q}^{m,n}\left(z\,\middle|\,\mathbf{a};\,\mathbf{b}\right)
= \sum_{k=1}^{m}
    \underbrace{
        \frac{\displaystyle\prod_{\substack{j=1\\j\ne k}}^{m}\!\!\Gamma(b_j-b_k)
                 \;\prod_{j=1}^{n}\Gamma(1+b_k-a_j)}
                 {\displaystyle\prod_{j=m+1}^{q}\!\!\Gamma(1+b_k-b_j)
                 \;\prod_{j=n+1}^{p}\!\!\Gamma(a_j-b_k)}
    }_{=:\;A_k}
    \;z^{b_k}\;
    {}_{p}F_{q-1}\left(\begin{matrix}1+b_k-a_1,\dots,1+b_k-a_p\\
        (1+b_k-b_j)_{j\ne k}\end{matrix}
        \;\middle|\;(-1)^{p-m-n}z\right)}
```

## Slater's Upper Expansion

When $p > q$, or when $p = q$ and $|z| > 1$, the contour is closed around
the poles of $\Gamma(1-a_j+s)$ for $j = 1, \dots, n$, yielding a series in $z^{-1}$:

```math
\boxed{G_{p,q}^{m,n}\left(z\,\middle|\,\mathbf{a};\,\mathbf{b}\right)
= \sum_{h=1}^{n}
    \underbrace{
        \frac{\displaystyle\prod_{\substack{j=1\\j\ne h}}^{n}\!\!\Gamma(a_h-a_j)
                 \;\prod_{j=1}^{m}\Gamma(1-a_h+b_j)}
                 {\displaystyle\prod_{j=n+1}^{p}\!\!\Gamma(1-a_h+a_j)
                 \;\prod_{j=m+1}^{q}\!\!\Gamma(a_h-b_j)}
    }_{=:\;B_h}
    \;z^{a_h-1}\;
    {}_{q}F_{p-1}\left(\begin{matrix}1-a_h+b_1,\dots,1-a_h+b_q\\
        (1-a_h+a_j)_{j\ne h}\end{matrix}
        \;\middle|\;(-1)^{q-m-n}z^{-1}\right)}
```

## Choosing the Expansion

- If $p < q$, use the lower expansion.
- If $p > q$, use the upper expansion.
- If $p = q$ and $|z| \le 1$, use the lower expansion.
- If $p = q$ and $|z| > 1$, use the upper expansion.

For balanced cases ($p = q$) near $|z| = 1$, `MeijerG.jl` promotes the Slater evaluation internally to `BigFloat` to maintain accuracy. Output type remains unaffected by this.

## Order Reduction

Before computation, the implementation removes cancelling $\Gamma$ factors:

- If $a_k = b_j$ for some $k \le n$ and $j > m$, remove both parameters and decrement $(p, q, n)$ by one.
- If $a_k = b_j$ for some $k > n$ and $j \le m$, remove both parameters and decrement $(p, q, m)$ by one.

## Adaptive Perturbation for Confluent Poles

When residue-based expansions (Slater lower/upper) encounter confluent poles the associated $\Gamma$ denominators vanish and the standard expansion formula becomes indeterminate. 

As mentioned above `MeijerG.jl` handles this by **adaptive perturbation**: a parameter-regularization technique that computes
the desired G-function as a limit. For confluent cases, a small perturbation $\delta \ll 1$ is applied to the parameter set (e.g., $b_j \to b_j + \delta\, \epsilon_j$ where $\epsilon_j$ are independent random variables), and the expansion is evaluated at increasing refinement levels. Richardson extrapolation on the results drives convergence to the true value without explicit limit taking or numerical differentiation.

The perturbation approach:

1. Computes the Slater expansion at a perturbed parameter set.
2. Refines by halving the perturbation magnitude.
3. Extrapolates successive results using Richardson coefficients to remove leading-order error terms.
4. Terminates when convergence reaches the user-requested tolerance (relative and absolute).

If perturbation fails to converge after a configurable number of refinement levels (default: 8), the package falls back to direct Mellin-Barnes contour integration using [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl).

## References

- [DLMF §16.17](https://dlmf.nist.gov/16.17)
- Slater, L.J. (1966), *Generalized Hypergeometric Functions*: [archive link](https://archive.org/details/generalizedhyper0000unse)
