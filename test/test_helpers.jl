relerr(x, y) = abs(x - y) / max(abs(y), eps(typeof(abs(y))))

isfinitevalue(x::Real) = isfinite(x)
isfinitevalue(x::Complex) = isfinite(real(x)) && isfinite(imag(x))
