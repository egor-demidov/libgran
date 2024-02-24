# Unary force models

libran comes with some build-in binary force models. Here, we provide class fererences
for reach built-in binary model template:

- [Triangular facet](Van-der-Waals-attraction.md) - a wrapper that simulates rigid triangular facets with particle-surface interaction
  forces described by surface force models:
  - [Van der Waals attraction](Van-der-Waals-attraction-surface.md) - simulates Van der Waals attraction between a
    particle and a plane based on [Hamaker theory](https://doi.org/10.1016/S0031-8914(37)80203-7)
  - [Frictional contact](Frictional-contact-surface.md) - simulates frictional forces between a particle
    and a plane in contact based on the approach described in [Luding 2008](https://doi.org/10.1007/s10035-008-0099-x)

<seealso>
<category ref="related">
   <a href="Overview.md">Overview</a>
    <a href="Installation.md">Installation</a>
    <a href="Tutorials.md">Tutorials</a>
    <a href="Class-reference.md">Class reference</a>
</category>
<category ref="external">
    <a href="https://github.com/egor-demidov/libgran">libgran on GitHub</a>
    <a href="https://github.com/egor-demidov/libtimestep">libtimestep on GitHub</a>
</category>
</seealso>