# Triangular facet

<tldr>
<p>A type that stores a variable number of references to surface force models</p>
<p>Name: <code>triangular_facet</code></p>
<p>Defined in: <code>&lt;libgran/surface_force/triangular_facet.h&gt;</code></p>
</tldr>

## Template parameters

{type="wide"}
`field_value_t`
: primary field type

`real_t`
: scalar field type

`surface_force_functors_t ...`
: list of surface force model types

## Constructor arguments

{type="wide"}
`field_value_t v_facet`
: initial translational velocity of the facet

`std::tuple<field_value_t, field_value_t, field_value_t> vertices`
: vertices of the triangular facet

`surface_force_functors_t & ...`
: list of surface force model references

## Public member functions

### operator ()

Synopsis:

Called by the [unary force container](Unary-force-container.md) to compute the acceleration of particle i due to
surface forces exerted on it because of the triangular facet

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
: translational acceleration and angular acceleration of particle i due to surface forces exerted on it by the
  triangular facet

### get_unit_normal()

Synopsis:

Returns the normal unit vector of the plane in which the triangular facet lies

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`field_value_t const &`
: const reference to the normal unit vector of the plane in which the triangular facet lies

### get_vertices()

Synopsis:

Returns the tuple with current positions of the vertices of the triangular facet

Arguments:

{type="wide"}
`void`
: no arguments

Return value:

{type="wide"}
`std::tuple<field_value_t, field_value_t, field_value_t> const &`
: const reference to the tuple containing the positions of the vertices of the triangular facet

### update_positions() 

Synopsis:

Needs to be called at every integration time step to update the positions of the vertices based on the
facet velocity, `v_facet`, unless the facet is static

Arguments:

{type="wide"}
`real_t dt`
: time step by which the vertex positions need to be advanced

{type="wide"}
`void`
: no return value

## Public member variables

{type="wide"}
`field_value_t v_facet`
: current translational velocity of the triangular facet

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
