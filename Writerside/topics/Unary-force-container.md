# Unary force container

<tldr>
<p>A record type that stores a variable number of references to unary force model
instances</p>
<p>Name: <code>unary_force_container_t</code></p>
<p>Defined in: <code>&lt;libgran/granular_system/granular_system.h&gt;</code></p>
</tldr>

## Template parameters

{type="wide"}
`field_value_t`
: primary field type

`real_t`
: scalar field type

`unary_force_functors_t ...`
: list of unary force model types

## Constructor arguments

{type="wide"}
`unary_force_functors_t & ...`
: list of unary force model references

## Public member functions

### operator ()

Synopsis:

Called by the [granular system](Granular-system.md) to compute the acceleration of particle i due to unary forces acting on it

Arguments:

{type="wide"}
`size_t i`
: index of the particle that is being accelerated

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
: translational acceleration and angular acceleration of particle i due to unary forces acting on it

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