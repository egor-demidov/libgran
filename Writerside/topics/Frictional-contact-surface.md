# Frictional contact (surface)

<tldr>
<p>A surface force model that simulates frictional forces between a particle and a plane in contact</p>
<p>Name: <code>surface_contact_force_functor</code></p>
<p>Defined in: <code>&lt;libgran/contact_force/surface_contact_force.h&gt;</code></p>
</tldr>

## Template parameters

{type="wide"}
`field_value_t`
: primary field type

`real_t`
: scalar field type

## Constructor arguments

{type="wide"}
`size_t n_part`
: number of particles in the system

`real_t k`
: normal stiffness

`real_t gamma_n`
: normal damping coefficient

`real_t k_t`
: tangential stiffness

`real_t gamma_t`
: tangential damping coefficient

`real_t mu_s`
: static friction coefficient for sticking/sliding

`real_t phi_d`
: dynamic-to-static friction ratio for sticking/sliding

`real_t k_r`
: rolling resistance stiffness

`real_t gamma_r`
: rolling resistance damping coefficient

`real_t mu_r`
: static friction coefficient for rolling resistance

`real_t phi_r`
: dynamic-to-static friction ratio for rolling resistance

`real_t k_o`
: torsion resistance stiffness

`real_t gamma_o`
: torsion resistance damping coefficient

`real_t mu_o`
: static friction coefficient for torsion resistance

`real_t phi_o`
: dynamic-to-static friction ratio for torsion resistance

`real_t r_part`
: particle radius

`real_t mass`
: particle mass

`real_t inertia`
: particle moment of inertia

`real_t dt`
: integration time step

`field_value_t field_zero`
: zero-valued primary field

`real_t real_zero`
: zero-valued scalar field

## Public member functions

### operator ()

Synopsis:

Called by the [triangular facet](Triangular-facet.md) to compute the acceleration
of particle i due to its frictional contact with the triangular facet

Arguments:

{type="wide"}
`size_t i`
: index of the particle that is being accelerated

`field_value_t const & x_facet`
: point on the facet closest to the center of particle i

`std::vector<fielv_value_t> const & x`
: const reference to the position buffer

`std::vector<fielv_value_t> const & v`
: const reference to the velocity buffer

`std::vector<fielv_value_t> const & theta`
: const reference to the orientation buffer

`std::vector<fielv_value_t> const & omega`
: const reference to the angular velocity buffer

`real_t t`
: simulation time

Return value:

{type="wide"}
`std::pair<field_value_t, field_value_t>`
: translational acceleration and angular acceleration of particle i due to its
frictional contact with particle j

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
