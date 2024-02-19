# Van der Waals attraction

<tldr>
<p>A binary force model that simulates Van der Waals attraction</p>
<p>Name: <code>hamaker_functor</code></p>
<p>Defined in: <code>&lt;libgran/hamaker_force/hamaker_force.h&gt;</code></p>
</tldr>

## Template parameters

{type="wide"}
`field_value_t`
: primary field type

`real_t`
: scalar field type

## Constructor arguments

{type="wide"}
`real_t A`
: Hamaker constant of the material

`real_t h0`
: saturation distance

`real_t r_part`
: particle radius

`real_t mass`
: particle mass

`field_value_t field_zero`
: zero-valued primary field

`real_t real_zero`
: zero-valued scalar field

## Public member functions

### operator ()

Synopsis:

Called by the [binary force container](Binary-force-container.md) to compute the acceleration
of particle i due to its Van der Waals attraction to particle j

Arguments:

{type="wide"}
`size_t i`
: index of the particle that is being accelerated

`size_t j`
: index of the particle that is acting on particle i

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
Van der Waals attraction to particle j

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
