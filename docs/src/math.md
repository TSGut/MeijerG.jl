# Mathematical Notes

`MeijerG.jl` evaluates the Meijer G-function via residue expansions derived
from the Mellin-Barnes contour integral representation. See
[DLMF §16.17](https://dlmf.nist.gov/16.17),
[Slater (1966)](https://archive.org/details/generalizedhyper0000unse), and
[Mellin-Barnes integral background](https://en.wikipedia.org/wiki/Mellin%E2%80%93Barnes_integral).

## Definition

The Meijer G-function is defined by the Mellin-Barnes contour integral ([DLMF §16.17](https://dlmf.nist.gov/16.17)):

```math
G_{p,q}^{m,n}\left(z\;\middle|\;\begin{matrix}a_1,\dots,a_p\\b_1,\dots,b_q\end{matrix}\right)
= \frac{1}{2\pi i}\int_{\mathcal{L}}
\frac{\displaystyle\prod_{j=1}^{m}\Gamma(b_j-s)\;\prod_{j=1}^{n}\Gamma(1-a_j+s)}
     {\displaystyle\prod_{j=m+1}^{q}\Gamma(1-b_j+s)\;\prod_{j=n+1}^{p}\Gamma(a_j-s)}
\,z^s\,ds
```

with $0 \le m \le q$, $0 \le n \le p$, and $z \ne 0$. The contour $\mathcal{L}$
separates poles of $\Gamma(b_j - s)$ on the right from poles of
$\Gamma(1 - a_j + s)$ on the left.

Rather than evaluating this integral numerically, the package uses the residue
theorem to obtain finite sums of generalized hypergeometric functions.

When poles in the selected residue family are confluent, meaning integer-separated
parameters in the active set, Slater's expansion acquires logarithmic terms.
`MeijerG.jl` evaluates these cases by a parameter-perturbation limit of the
same residue representation, which recovers the finite logarithmic contribution.

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

## Order Reduction

Before computing, the implementation removes cancelling $\Gamma$ factors:

- If $a_k = b_j$ for some $k \le n$ and $j > m$, remove both parameters and decrement $(p, q, n)$ by one.
- If $a_k = b_j$ for some $k > n$ and $j \le m$, remove both parameters and decrement $(p, q, m)$ by one.

## Special Case Verification

Many classical special functions are Meijer G special cases. Examples:

- $e^x$: $G_{0,1}^{1,0}\left(-x\,\middle|\,\begin{matrix}-\\0\end{matrix}\right)$
- $\sin x$: $\sqrt{\pi}\,G_{0,2}^{1,0}\left(\tfrac{x^2}{4}\,\middle|\,\begin{matrix}-\\\tfrac{1}{2},0\end{matrix}\right)$
- $\cos x$: $\sqrt{\pi}\,G_{0,2}^{1,0}\left(\tfrac{x^2}{4}\,\middle|\,\begin{matrix}-\\0,\tfrac{1}{2}\end{matrix}\right)$
- $J_\nu(x)$: $\left(\tfrac{x}{2}\right)^\nu G_{0,2}^{1,0}\left(\tfrac{x^2}{4}\,\middle|\,\begin{matrix}-\\\tfrac{\nu}{2},-\tfrac{\nu}{2}\end{matrix}\right)$
- $K_\nu(x)$: $\tfrac{1}{2}G_{0,2}^{2,0}\left(\tfrac{x^2}{4}\,\middle|\,\begin{matrix}-\\\tfrac{\nu}{2},-\tfrac{\nu}{2}\end{matrix}\right)$
- $\gamma(\alpha, x)$: $G_{1,2}^{1,1}\left(x\,\middle|\,\begin{matrix}1\\\alpha,0\end{matrix}\right)$
- $\Gamma(\alpha, x)$: $G_{1,2}^{2,0}\left(x\,\middle|\,\begin{matrix}1\\\alpha,0\end{matrix}\right)$

## References

- [DLMF §16.17](https://dlmf.nist.gov/16.17)
- [mpmath Meijer G documentation](https://mpmath.org/doc/current/functions/hypergeometric.html#meijer-g-function)
- [Meijer G (Wikipedia)](https://en.wikipedia.org/wiki/Meijer_G-function)
- [Mellin-Barnes / Barnes integral (Wikipedia)](https://en.wikipedia.org/wiki/Mellin%E2%80%93Barnes_integral)
- Slater, L.J. (1966), *Generalized Hypergeometric Functions*: [archive link](https://archive.org/details/generalizedhyper0000unse)
