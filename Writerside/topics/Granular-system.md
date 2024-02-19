# Granular system

<tldr>
<p>Name: <code>granular_system</code></p>
<p>Defined in: <code>&lt;libgran/granular_system/granular_system.h&gt;</code></p>
<p>An object that represents a granular system with force models</p>
</tldr>

## Template parameters

{type="wide"}
`field_value_t`
: primary field type

`real_t`
: scalar field type

`integrator_t`
: integrator type

`step_handler_t`
: step handler type

`binary_force_contatiner_t`
: binary force container type

`unary_force_container_t`
: unary force container type

## Member functions

### Constructor

Arguments:

{type="wide"}
`std::vector<primary_field_t> x0`
: Buffer with initial positions

`std::vector<primary_field_t> v0`
: Buffer with initial velocities

`std::vector<primary_field_t> theta0`
: Buffer with initial orientations

`std::vector<primary_field_t> omega0`
: Buffer with initial angular velocities

`real_t t0`
: Time at which the simulation starts

`field_value_t field_zero`
: Zero-valued primary field

`real_t real_zero`
: Zero-valued scalar field

`step_handler_t & step_handler`
: Reference to step handler

`binary_force_contatiner_t binary_force_functors`
: Binary force model container

`unary_force_contatiner_t unary_force_functors`
: Unary force model container

### do_step

Synopsis:

Advances the granular system by one time step

Arguments:

{type="narrow"}
`realt_t dt`
: time step to advance the simulation by

Return value:

{type="narrow"}
`void`
: no return value

### compute_accelerations (binary)

Synopsis:

Computes the translational and angular accelerations of particle i due to its
binary interactions with particle j

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

### compute_accelerations (unary)

Synopsis:

Computes the translational and angular accelerations of particle i due to 
the unary forces acting on it

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
: translational acceleration and angular acceleration of particle i due to the unary force
models acting on it

### get_x (inherited)

Synopsis:

Returns a const reference to the particle position buffer

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`std::vector<field_value_t> const &`
: const reference to the particle position buffer

### get_v (inherited)

Synopsis:

Returns a const reference to the particle velocity buffer

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`std::vector<field_value_t> const &`
: const reference to the particle velocity buffer

### get_theta (inherited)

Synopsis:

Returns a const reference to the particle orientation buffer

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`std::vector<field_value_t> const &`
: const reference to the particle orientation buffer

### get_omega (inherited)

Synopsis:

Returns a const reference to the particle angular velocity buffer

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`std::vector<field_value_t> const &`
: const reference to the particle angular velocity buffer

### reset_acceleration_buffers (inherited)

Synopsis:

Sets all the values in translational and angular acceleration buffers to zero

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`void`
: no return value

<seealso>
<!--<category ref="related">
           <a href="Links.topic">Topic about links</a>
           <a href="Some_other.topic"/>
</category>-->
</seealso>
