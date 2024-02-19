# Binary force container

<tldr>
<p>Name: <code>binary_force_container_t</code></p>
<p>Defined in: <code>&lt;libgran/granular_system/granular_system.h&gt;</code></p>
<p>A record type that stores a variable number of references to binary force model
instances</p>
</tldr>

## Template parameters

{type="wide"}
`field_value_t`
: primary field type

`real_t`
: scalar field type

`binary_force_functors_t ...`
: list of binary force model types

## Member functions

### Constructor

Arguments:

{type="wide"}
`binary_force_functors_t & ...`
: list of binary force model references

### operator ()

Synopsis:

Called by the granular system to compute the acceleration of particle i due to binary interactions
with other particles

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
: translational acceleration and angular acceleration of particle i due to its interaction
with particle j
