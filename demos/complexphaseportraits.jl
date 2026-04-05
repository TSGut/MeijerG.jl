using MeijerG
using ComplexPhasePortrait
using Images
using ColorSchemes

# Grid resolution
xs = range(-5.0, 5.0; length=1000)
ys = range(-5.0, 5.0; length=1000)
Z  = [x + im*y for y in ys, x in xs]

# Build a periodic Viridis-style colormap so the phase wrap closes continuously.
function periodic_colormap(colors)
	base = [RGB(c.r, c.g, c.b) for c in colors]
	return vcat(base, base[end-1:-1:2])
end

viridis_colormap = periodic_colormap(ColorSchemes.viridis.colors)

# Ensure output folder exists (relative to this script)
figures_dir = joinpath(@__DIR__, "figures")
mkpath(figures_dir)

# 1. Bessel‑related Meijer G
f_bessel(z) = meijerg((), (0.0,0.0), 2, 0, z)
img1 = portrait(f_bessel.(Z), PTcgrid, ctype=:nist);
img2 = portrait(f_bessel.(Z), PTcgrid, colormap=viridis_colormap);
save(joinpath(figures_dir, "bessel_phase.png"), img1);
save(joinpath(figures_dir, "bessel_phase_viridis.png"), img2);

# 2. Confluent-pole exotic Meijer G
f_confluent(z) = meijerg((1/2,1/2,-1/4,3.0), (0.0,1/2,5,-3/2), 2, 1, z)
img3 = portrait(f_confluent.(Z), PTcgrid);
img4 = portrait(f_confluent.(Z), PTcgrid, colormap=viridis_colormap);
save(joinpath(figures_dir, "confluent_phase.png"), map(clamp01nan,img3));
save(joinpath(figures_dir, "confluent_phase_viridis.png"), map(clamp01nan,img4));

# 3. Strongly asymmetric spiral Meijer G
f_asym1(z) = meijerg((1/4,), (0.0, -2.2, 7/5), 1, 0, z)
img5 = portrait(f_asym1.(Z), PTcgrid);
img6 = portrait(f_asym1.(Z), PTcgrid, colormap=viridis_colormap);
save(joinpath(figures_dir, "asym1_phase.png"), img5);
save(joinpath(figures_dir, "asym1_phase_viridis.png"), img6);

# 4. Rational‑symmetry Spiral Meijer G
f_rational(z) = meijerg((1/7,2/7,4/7), (3/7,5/7,6/7), 2, 1, z)
img7 = portrait(f_rational.(Z), PTcgrid);
img8 = portrait(f_rational.(Z), PTcgrid, colormap=viridis_colormap);
save(joinpath(figures_dir, "rational_phase.png"), img7);
save(joinpath(figures_dir, "rational_phase_viridis.png"), img8);
