# libgran ![Unit tests workflow badge](https://github.com/egor-demidov/libgran/actions/workflows/test.yml/badge.svg) ![Documentation deployment workflow badge](https://github.com/egor-demidov/libgran/actions/workflows/build-docs.yml/badge.svg)

libgran is a Discrete Element Method (DEM) framework for simulating the mechanical behavior of soot aggregates. DEM is a
technique for simulation of granular media consisting of rigid spherical particles. The resultant force and torque acting
on each particle are computed and used with Newton's second law to compute the motion of particles:
```math
m\ddot{\mathbf{x}}=\mathbf{f}
```
```math
I\ddot{\boldsymbol{\omega}}=\boldsymbol{\tau}
```
The forces that particles experience arise from friction at inter-particle contacts, bonding between particles, 
inter-particle attraction, field forces, etc. libgran contains a bonded and a non-bonded contact model, a Van der Waals
attraction model and is designed to be easily extensible with custom models. A simulation is set up in the driver program, which needs to initialize three components:
- Force functor container
- Step handler
- Granular system

### Force functor container

A force functor container is an object that contains instances of all force models that are used in the simulation.
The choice of which force models to use is made statically and the list of force models is provided in template arguments
to the force functor container. It is a good practice to create aliases to the force models and the force container types
that will be used in the simulation at the start of the driver program to improve code readability later on. For example,
the code snippet below creates an alias to a force functor container with a frictional contact force model and Van der Waals
attraction force model called `binary_force_container_t`:
```c++
#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>; // Contact force
using vdw_force_functor_t = hamaker_functor<Eigen::Vector3d, double>; // Van der Waals force
using binary_force_container_t = binary_force_functor_container<Eigen::Vector3d, double, contact_force_functor_t, vdw_force_functor_t>; // Binary force container
```
Then, in the driver program, each force model and the force functor container need to be instantiated:
```c++
int main() {
    
    /* Declare the constants ... */
    
    // Create an instance of contact force model
    contact_force_functor_t contact_force_model(x0.size(),
        k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of Van der Waals attraction model
    vdw_force_functor_t hamaker_model(A, h0,
        r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    binary_force_container_t
            binary_force_functors{contact_force_model, hamaker_model};
    
    /* Do other simulation stuff ... */
    
    return 0;
}
```
`binary_force_container_t` is an alias to `binary_force_functor_container` - an object meant to contain force models that describe
binary interactions between particles. It is also possible to define unary force models for forces that affect each particle
independently. Examples of such forces could be gravity or an electric field. If we do not want to use any unary forces in the
simulation, we will need to create an en empty placeholder unary force container:
```c++
/* ... */
using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>; // Unary force container (empty)

int main () {
    
    /* ... */

    unary_force_container_t unary_force_functors;
    
    /* ... */
    
    return 0;
}
```
### Step handler
A step handler is an object that provides a layer of abstraction between the integrator that computes the position and velocity
increments, and the representation of positions / velocities of individual particles. A trivial step handler simply takes 
position or velocity increments as arguments and increments the position and velocity of the specified particle. But in some
simulations this default behavior needs to be overriden or additional manipulations need to be performed with position / velocity
increments computed for a particle. In those scenarios, it is convenient to define a step handler. For example, if we want to keep
one particle in the simulation fixed in space, we can define a step handler that never increments velocity of that particle.

The default step handler that simply increments positions and velocities is provided by the implementation and can be created
in the following manner:
```c++
int main() {
    
    /* ... */
    
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;
    
    /* ... */
    
    return 0;
}
```
### Granular system
A granular system is an object that contains positions and velocities of the particles and encapsulates the integrator.
It makes calls to the force functor container when it needs accelerations computed and to the step handler when it needs
positions incremented. Initial positions and velocities of the particles need to be provided to the granular system constructor.
Granular system contains a method, `do_step(dt)`, that should be called by the driver program to advance the system by
time step dt. Example of usage of granular system:
```c++
/* ... */

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>; // Granular system representation
    
int main() {
    
    /* ... */

    // Initialize two particles
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    x0.emplace_back(0.0, 0.0, 0.0);
    x0.emplace_back(0.0, 2.5*r_part, 0.0);
    
    // Initialize the remaining buffers
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
`   std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Buffers with initial positions anv velocities are copied into the binary system object
    // and can be safely deleted after initialization of granular system
    // Step handler and force functor containers are stored in the granular system as references
    // and need to exist for the duration of use of the granular system object
    granular_system_t system(x0,
        v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(),
        0.0, step_handler_instance, binary_force_functors, unary_force_functors);

    /* ... */
    
    return 0;
}
```
Now that the granular system is initialized, all that's left is to perform the time stepping, output data periodically,
and post-process the results as necessary:
```c++
/* ... */

int main() {

    /* ... */

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            /* Perform a data dump */
        }
        system.do_step(dt); // Advance the system
    }
    
    /* Post-process the data */
    
    return 0;
}
```
### libgran work flow
The work flow involved in a simulation and the relationships between components are presented in the diagram below:
```mermaid
sequenceDiagram
    participant Driver program
    Driver program->>Force functor container: Initialize the force functors
    Driver program->>Step handler: Initialize the step handler
    Driver program->>Granular system: Initialize the granular system
    loop Perform time stepping
        Driver program->>Granular system: call do_step(dt)
        loop Integrate each particle
            loop Add up contributions from other particles
                Granular system->>Force functor container: call operator ()
                Force functor container->>Granular system: Return computed acceleration
            end
            Granular system->>Step handler: Increment position and velocity
        end
    end
    participant Granular system
    participant Step handler
    participant Force functor container
```
## Included force models

The force models that are provided by the implementation are described in this section.

### Contact force

The contact model is based on constraining the four degrees of freedom of motion of two particles relative to each other:
normal translation, tangential translation, torsion, and rolling. Let us consider two particles i and j. We can begin by 
defining a unit normal vector:
```math
\mathbf{n}=\frac{\mathbf{x}_{j}-\mathbf{x}_{i}}{\lVert\mathbf{x}_{j}-\mathbf{x}_{i}\rVert}
```
The relative velocity at the point of contact is:
```math
\mathbf{v}_{ij}=\mathbf{v}_j-\mathbf{v}_i+\boldsymbol{\omega}_{j}\times a\mathbf{n}+\boldsymbol{\omega}_{i}\times a\mathbf{n}
```
where $a$ is particle radius corrected for inter-particle overlap/separation, $\mathbf{v}$ is particle translational velocity, and
$\boldsymbol{\omega}$ is particle angular velocity. Relative velocity can be decomposed into normal and residual (tangential) components:
```math
\mathbf{v}_{ij,\rm n}=\left(\mathbf{v}_{ij}\cdot \mathbf{n}\right)\mathbf{n}
```
```math
\mathbf{v}_{ij,\rm t}=\mathbf{v}_{ij}-\mathbf{v}_{ij,\rm n}
```
Similarly, relative angular velocity:
```math
\boldsymbol{\omega}_{ij}=\boldsymbol{\omega}_j-\boldsymbol{\omega}_i
```
can be decomposed into normal (torsional) and residual (rolling) components:
```math
\boldsymbol{\omega}_{ij,\rm o}=\left(\boldsymbol{\omega}_{ij}\cdot\mathbf{n}\right)\mathbf{n}
```
```math
\boldsymbol{\omega}_{ij,\rm r}=\boldsymbol{\omega}_{ij}-\boldsymbol{\omega}_{ij,\rm o}
```
To constrain the four degrees of freedom, we insert four springs, as illustrated in the figure below: 
![Illustration of the four degrees of freedom that we would like
to constrain and the springs that are inserted for each DOF](resources/images/degrees-of-freedom.svg)

The length of the normal spring, $\delta$, can be computed directly at any point in the simulation from positions of the particles,
$\mathbf{x}$, and their radius, $r$:
```math
\delta=\lVert\mathbf{x}_j-\mathbf{x}_i\rVert-2r
```
The remaining three springs have zero length at the time the contact is formed and have their lengths accumulated
throughout the duration of the contact. Let a spring vector be $\boldsymbol\xi$. Then the rate of stretching /
contraction of a spring, $\dot{\boldsymbol\xi}$, is given by:
```math
\dot{\boldsymbol\xi}_{\rm t}=\mathbf{v}_{ij,\rm t}
```
```math
\dot{\boldsymbol\xi}_{\rm o}=r\boldsymbol{\omega}_{ij,\rm o}
```
```math
\dot{\boldsymbol\xi}_{\rm r}=\boldsymbol{\omega}_{ij,\rm r}\times a\mathbf{n}
```
Then, as long as the contact lasts, spring $\boldsymbol{\xi}$ is incremented at each time step to
obtain a new spring, $\boldsymbol{\xi}'$, to be used at the next time step:
```math
\boldsymbol{\xi}'=\boldsymbol{\xi}+\dot{\boldsymbol{\xi}}\Delta t
```

#### Bonded contact force

For a pair of particles that is connected with a rigid bond, we would like to approximate a rigid-body motion.
In other words, the common reference frame of particles i and j can rotate and translate, but any translation or
rotation of particle i relative to particle j should be restricted. The distance between all points in the pair of
particles should be approximately preserved over time. It can be shown that when using the springs defined in the
[contact force section](#contact-force) to restrict the motion of particles i and j, then the union of particles i and j will undergo
rigid body motion as the stiffness of inserted springs approaches infinity. In the simulation we need to use a finite
stiffness value, but as long as the amplitude of oscillations is much smaller than the length scale of particles in the
simulation, the motion will, approximately, be rigid.

To stabilize the system over time and dissipate any vibrational kinetic energy in the bonds, each spring is supplemented by a
dashpot element. The force, $\mathbf{f}$, exerted on particle i by each spring in the i-j bond is given by:
```math
\mathbf{f}=k\boldsymbol{\xi}+\gamma\dot{\boldsymbol{\xi}}
```
where $k$ is stiffness and $\gamma$ is the damping coefficient of the respective spring.
And forces exerted on particle j are equal in magnitude and opposite in direction.
Forces arising from the normal and tangential springs are applied to the particles.
The force associated with the tangential spring will also give rise to torques because the tangential force is not
collinear with $\mathbf{n}$. Forces computed from the torsion and rolling resistance springs are quasi-forces
that are not applied to particles i and j, but are only used to compute torques that will be applied to the particles.

#### Frictional contact force

The model described in [Luding 2008](https://doi.org/10.1007/s10035-008-0099-x) is used to simulate
frictional contacts between non-bonded particles. A brief description
is provided here and the reader is referred to [Luding 2008](https://doi.org/10.1007/s10035-008-0099-x)
for a more detailed description. Luding's model uses the same four springs described
in an earlier section to compute normal and tangential forces, rolling
and torsion resistance torques. Instead of directly setting force proportional
to spring elongation, a certain degree of slip is allowed between particles in contact.
That is done by computing a test force, $\mathbf{f}_{0}$:
```math
\mathbf{f}_{0}=k\boldsymbol{\xi}+\gamma\dot{\boldsymbol{\xi}}
```
and deciding whether static or dynamic friction should be used based on Coulomb's law of friction:
```math
f_{C,s}=\mu_s f_{\rm n}
```
```math
f_{C,d}=\mu_d f_{\rm n}
```
where $\mu_s$ is the static friction coefficient, $\mu_d$ is the dynamic friction coefficient,
and $f_n$ is the magnitude of the normal force between the two particles.
Then a choice is made whether static or dynamic friction should be used based
on the following condition:
```math
\text{if}\ \lVert\mathbf{f}_0\rVert\leq f_{C,s}\ \text{use static friction}
```
```math
\text{if}\ \lVert\mathbf{f}_0\rVert> f_{C,s}\ \text{use dynamic friction}
```

In case the contact is determined to be in the state of static friction,
the tangential spring $\boldsymbol\upxi$ is incremented as described in
[contact force section](#contact-force) and the test force is used to compute torques and,
in the case of the tangential force, is also applied to the contacting particles.
In case the contact is determined to be in the state of dynamic friction,
the spring is not allowed to stretch anymore. Instead, a certain extent of slipping is allowed.
The spring to be used at the next iteration, $\boldsymbol\upxi'$, is set to
```math
\boldsymbol\upxi'=-\frac{1}{k}\left(f_{C,d}\frac{\mathbf{f}_0}{\lVert \mathbf{f}_0\rVert}+\gamma\dot{\boldsymbol\upxi}\right)
```
and the magnitude of the Coulomb's force, $f_{C,d}$, is used to compute torques and, if applicable, accelerate the particles in contact.

The model is only enabled when the normal force is repulsive. 
Once particles are not overlapping, the frictional force is set
to zero and accumulated springs are reset.

### Van der Waals attraction force

[Hamaker 1937](https://doi.org/10.1016/S0031-8914(37)80203-7) derived the potential energy, $U$, due
to Van der Waals attraction between two particles of radius $r$ whose surfaces are separated by distance
$\delta$ from each other to be:
```math
U=-\frac{A}{6}\left[\frac{2r^2}{(4r+\delta)\delta}+\frac{2r^2}{(2r+\delta)^2}+\ln\frac{(4r+\delta)\delta}{(2r+\delta)^2}\right]
```
where $A$ is the Hamaker constant - a material property. The magnitude of force acting on the particles can be derived by
differentiating potential energy, $U$, with respect to separation distance, $\delta$, and the direction of the force will
coincide with the normal unit vector $\bf n$ defined earlier:
```math
\mathbf{f}=-\frac{A}{6}\left[\frac{(4r+2\delta)}{(4r+\delta)\delta}-\frac{2}{(2r+\delta)}-\frac{4r^2}{(2r+\delta)^3}-\frac{2r^2(4r+2\delta)}{(4r+\delta)^2\delta^2}\right]\mathbf{n}
```
Since in the limit as $\delta$ approaches $0$ the magnitude of force (and potential energy) becomes infinite,
a saturation distance $\delta_0$ is introduced. Attractive force does not increase past the saturation distance. The magnitude of $\delta_0$ varies between 0.4 and 1 nm ([Ranade 1987](
https://doi.org/10.1080/02786828708959155)).

## Simulation example

In this section, we set up a simulation where two cubes made up of particles held together by friction and
Van der Waals attraction are assembled and collided with each other. We will use libeigen for linear algebra
functionality. If your project is a git repository, the dependencies for this example can be installed as
git submodules:
```shell
git submodule add https://gitlab.com/libeigen/eigen.git deps/eigen
git submodule add https://github.com/egor-demidov/libgran deps/libgran
git submodule add https://github.com/egor-demidov/libtimestep deps/libtimestep
```
Then add the include directories in your CMakeLists.txt file:
```cmake
include_directories(deps/libtimestep/include)
include_directories(deps/libgran/include)
include_directories(deps/eigen)
```
See [installation](#installation) section for more information about compiler flags and parallelization capabilities.

Now that the dependencies are met, we can start writing the driver program. In your main.cpp, include
the necessary headers:
```c++
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>

// Using Eigen for linear algebra
#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>
```
In addition to some STL headers and libeigen, we included three libgran headers. Two of them define the force models
that we want to use in this simulation, contact force and Van der Waals force, and the thord header defines
the granular system object.

Now, let us "assemble" the templates that we will use and create shorter aliases for the produced types:
```c++
using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>; // Contact force
using vdw_force_functor_t = hamaker_functor<Eigen::Vector3d, double>; // Van der Waals force
using binary_force_container_t
    = binary_force_functor_container<Eigen::Vector3d, double,
    contact_force_functor_t, vdw_force_functor_t>; // Binary force container

using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>; // Unary force container (empty)

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>; // Granular system representation
```
Okay, now we can proceed to implement the main function. After the alias definitions, write:
```c++
int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 1.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;
    const size_t n_thermo_dumps = 10000;
    const size_t thermo_dump_period = n_steps / n_thermo_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 1.0e-9;
    
    return 0;
}
```
We just declared and initialized the simulation constants. They include general parameters, such as time step size,
number of time steps, number of data dumps, etc. The constants in this section also contain parameters
for the force models that we are going to use. The list of parameters that each force models requires is provided
in the [class reference](#class-reference) section. If `M_PI` is undefined, add the `_USE_MATH_DEFINES` compile definition.

Now, after the definition of constants, we can initialize the initial conditions for the simulation. 
We need to provide four buffers populated with initial values, each of size $N$, where $N$ is the number of particles
in the system. These buffers correspond to positions, translational velocities, orientations, and angular velocities
respectively. We will declare four corresponding `std::vector`s: `x0`, `v0`, `theta0`, and `omega0`:
```c++
std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
```
Let us build two cubes, with each side of a cube containing three particles. The cubes will have initial
translational velocities directed towards each other.
```c++
const double x_offset = 10.0 * r_part;
const double y_offset = 3.0 * r_part;

for (size_t i = 0; i < 3; i ++) {
    auto x = double(i) * 2.0 * r_part;
    for (size_t j = 0; j < 3; j ++) {
        auto y = double(j) * 2.0 * r_part;
        for(size_t n = 0; n < 3; n ++) {
            auto z = double(n) * 2.0 * r_part;
            
            // Add particle to cube 1
            x0.emplace_back(x, y, z);
            // Set velocity of the newly created particle
            v0.emplace_back(0.5, 0.0, 0.0);
            
            // Add particle to cube 2
            x0.emplace_back(x + x_offset, y + y_offset, z);
            // Set velocity of the newly created particle
            v0.emplace_back(-0.5, 0.0, 0.0);
        }
    }
}
```
Now, let us set the remaining buffers - orientations and angular velocities - to zero. It is important
that **all initial condition buffers are initialized** and are of the same size.
```c++
 // Initialize the remaining buffers
theta0.resize(x0.size());
omega0.resize(x0.size());
std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());
```
The next step is to create instances of the force models, force model containers, and the step handler:
```c++
// Create an instance of contact force model
contact_force_functor_t contact_force_model(x0.size(),
    k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
    mu_o, phi, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

// Create an instance of Hamaker model
vdw_force_functor_t hamaker_model(A, h0,
    r_part, mass, Eigen::Vector3d::Zero(), 0.0);

// Create an instance of binary force container
binary_force_container_t
    binary_force_functors{contact_force_model, hamaker_model};

// Create an instance of unary force container (empty)
unary_force_container_t
    unary_force_functors;

// Create an instance of step_handler
rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
    step_handler_instance;
```
At this point, we have everything we need to assemble our granular system object. Add the following line
after instantiation of force models:
```c++
granular_system_t system(x0,
    v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
    step_handler_instance, binary_force_functors, unary_force_functors);
```
Now, all that is left is to perform the time-stepping, perform periodic data dumps, and post-process
the results. Add an integration loop:
```c++
for (size_t n = 0; n < n_steps; n ++) {
    if (n % dump_period == 0) {
        /* TODO: add some code to dump data */
    }
    system.do_step(dt);
}
```
To advance the system by on time step of size dt, we simply need to call the `do_step(dt)` method
of the granular system object. If you run the program now, the simulation will be performed but no output will
be produced. We need to add some code to write out the state of the system periodically. Above your
main function, add a new declaration:
```c++
void dump_particle_positions(std::string const & dir, size_t count,
    std::vector<Eigen::Vector3d> const & x, double r_part) {
    
    std::stringstream out_file_name;
    out_file_name << dir << "/particles_" << count << ".csv";
    std::ofstream ofs(out_file_name.str());
    
    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    ofs << "x, y, z\n";
    for (auto const & point : x) {
        ofs << point[0] / r_part << ", " << point[1] / r_part << ", " << point[2] / r_part << "\n";
    }
}
```
Now back in your integration loop, where you left space for data dumping code, add:
```c++
dump_particle_positions("run", n / dump_period, system.get_x(), r_part);
```
Now, create a directory named run in the same directory as the executable and run it. You should see
the run directory get populated with 300 .csv files, each containing the coordinates of all particles at the
corresponding time step.

### Visualizing the results

Generated data files can be visualized using [ParaView](https://www.paraview.org/download/).
To create a visualization, open the files in ParaView. Press `File->Open` and navigate to
the directory with generated data files. ParaView will automatically recognize
sequentially numbered files a time series:

![Screenshot of ParaView file dialog](resources/images/01_pv_open.png)

Select the entry and press OK. ParaView will ask to select the appropriate data reader.
Select the CSV reader and press OK:

![Prompt to select the data reader](resources/images/02_pv_reader.png)

ParaView will load the data as a table. Press the green Apply button in the left panel.
A table previewing the data should appear on the right side of the screen. This preview 
can be closed.

![Previewing the loaded data](resources/images/03_pv_applied.png)

Now, we need to tell ParaView to interpret the rows in the loaded table as x-y-z points.
Select the loaded data set in the left panel and got to `Filters->Seacrh`. Start typing
`table to points` and select the corresponding filter when it appears among search results.

![Adding the table to points filter](resources/images/05_filter.png)

Now in the left panel there should be a filter that is a child of the table. Click
on it and in the dropdowns, select x column, y column, and z column to be x, y, and z
respectively.

![Selecting the x, y, and z columns](resources/images/06_converted.png)

Press the green Apply button. With the filter selected in the left panel, click the
`Glyph` button in the ribbon. In the left panel, set Glyph Type to sphere, scale array 
and orientation array to none, scale factor to 2, and Glyph Mode to all points. Click the 
green Apply button. Cubes made up of spherical particles should appear:

![Rendered cubes made up of spherical particles](resources/images/07_ready.png)

In the top ribbon, you can see that the current frame is 0 out of 300. You can press the
green play button to preview the animation, or directly type in 300 in the text box to jump 
to the last frame. If you do, you will see that at the end of the simulation, the cubes had collided
and are now stuck together.

![Collided cubes at the end of the simulation](resources/images/08_final_state.png)

## Implementing custom binary force models

Let us create a custom binary force model - linear pairwise attraction.
Particle i will be attracted to particle j as:
```math
\mathbf{f}=k_{\rm a}\delta\mathbf{n}
```
In your project, create a file `linear_attraction.h` and add include guards to it:
```c++
#ifndef LINEAR_ATTRACTION_H
#define LINEAR_ATTRACTION_H

/* Code will go here */

#endif //LINEAR_ATTRACTION_H
```
Now, between the include guards, add the template declaration:
```c++
template<typename field_value_t, typename real_t>
struct linear_attraction_functor {
    
};
```
The template parameter `field_value_t` is the type of the primary field (like `Eigen::Vector3d`
or `Eigen::Vector3f`) and the parameter `real_t` is the type of a real-valued scalar value (like `double` or `float`).

Now, let us implement a constructor. In order compute accelerations, we need to know the
$r$, $k_{\rm a}$, and $m$. We will also need to know the zero value of the primary field.
Inside the struct, declare the private member constants and the 
constructor. The constructor only initializes the member constants:

```c++
linear_attraction_functor(real_t r, real_t k_a, real_t mass, field_value_t zero_field) :
    r(r), k_a(k_a), mass(mass), zero_field(zero_field) {}

private:
const real_t r, k_a, mass;
const field_value_t zero_field;
```

The algorithm is then as follows:
```math
\mathbf{n}=\frac{\mathbf{x}_j-\mathbf{x}_i}{\lVert\mathbf{x}_j-\mathbf{x}_i\rVert}
```
```math
\delta=\lVert\mathbf{x}_j-\mathbf{x}_i\rVert - 2r
```
```math
\mathbf{a}=\frac{k_{\rm a}\delta}{m}\mathbf{n}
```

In your struct, after the constructor and before the `private` keyword, add:
```c++
std::pair<field_value_t, field_value_t> operator () 
   (size_t i, size_t j,
    std::vector<field_value_t> const & x,
    std::vector<field_value_t> const & v [[maybe_unused]],
    std::vector<field_value_t> const & theta [[maybe_unused]],
    std::vector<field_value_t> const & omega [[maybe_unused]],
    real_t t [[maybe_unused]]) const {
    
    /* Compute and return the accelerations */
}
```
This lengthy statement is the signature of `operator ()` - the function that will be called by the 
integrator to compute the acceleration acing on particle i due to its linear attractive interaction
with particle j. The arguments are: index of particle i, index of particle j, position buffer,
velocity buffer, orientation buffer, angular velocity buffer, and time. All these arguments must be accepted
by the `operator ()` function even if they are not used. In case an argument is not used, it should be
marked with the `[[maybe_unused]]` attribute to avoid causing compiler warnings.

Now, all that is left to do is to compute the acceleration and return it. Note that the return type
is an `std::pair`, because we need to return two accelerations: translational and angular. However, we assume that
this linear attractive model is a body force acting on the particles and does not directly produce torques. Therefore,
the angular acceleration that we return will be zero. Inside your `operator ()`, write:
```c++
field_value_t n = (x[j] - x[i]).normalized(); // Compute the normal unit vector
real_t delta = (x[j] - x[i]).norm() - 2.0 * r; // COmpute the separation between particles
field_value_t a = k_a * delta * n / mass; // Compute the translational acceleration
return std::make_pair(a, zero_field); // Return the translational and angular accelerations
```

Now, create a simple driver program to test this newly created binary force model
([driver program tutorial](#simulation-example)). We will use a file with geometry of a fractal
aggregates, [aggregate.h](https://github.com/egor-demidov/libgran/blob/cbd7c7d3423d52fd55a95de6b81c56775a852624/resources/aggregate.h), to initialize this simulation. It needs to
be downloaded and added to your project.

```c++
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>

// Using Eigen for linear algebra
#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "linear_attraction.h" // Include your custom model

// Include the fractal aggregate used as initial conditions
#include "aggregate.h"

using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>; // Contact force
using vdw_force_functor_t = hamaker_functor<Eigen::Vector3d, double>; // Van der Waals force
using lienar_force_functor_t = linear_attraction_functor<Eigen::Vector3d, double>; // Your custom model
using binary_force_container_t =
    binary_force_functor_container<Eigen::Vector3d, double,
    contact_force_functor_t, vdw_force_functor_t, lienar_force_functor_t>; // Binary force container

using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>; // Unary force container (empty)

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>; // Granular system representation

void dump_particle_positions(std::string const & dir, size_t count,
    std::vector<Eigen::Vector3d> const & x, double r_part) {

    std::stringstream out_file_name;
    out_file_name << dir << "/particles_" << count << ".csv";
    std::ofstream ofs(out_file_name.str());

    if (!ofs.good()) {
        std::cerr << "Unable to create a dump file at " << out_file_name.str() << std::endl;
        exit(EXIT_FAILURE);
    }

    ofs << "x, y, z\n";
    for (auto const & point : x) {
        ofs << point[0] / r_part << ", " << point[1] / r_part << ", " << point[2] / r_part << "\n";
    }
}

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 1.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;
    const size_t n_thermo_dumps = 10000;
    const size_t thermo_dump_period = n_steps / n_thermo_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 1.0e-9;
    
    // Parameter for your custom model
    const double k_attr = 4.0e-8;

    // Declare the initial condition buffers
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;

    // Copy the particle coordinates into the initial position buffer
    x0.resize(std::size(aggregate));
    std::transform(std::begin(aggregate), std::end(aggregate), x0.begin(), [r_part] (auto const & particle) -> Eigen::Vector3d {
        return Eigen::Vector3d{particle[0], particle[1], particle[2]} * r_part;
    });

    // Fill the remaining buffers with zeros
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Create an instance of contact force model
    contact_force_functor_t contact_force_model(x0.size(),
        k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o,
        mu_o, phi, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of Hamaker model
    vdw_force_functor_t hamaker_model(A, h0,
        r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    // Create an isntance of your custom model
    lienar_force_functor_t linear_model(r_part, k_attr, mass, Eigen::Vector3d::Zero());

    // Create an instance of binary force container
    binary_force_container_t
        binary_force_functors{contact_force_model, hamaker_model, linear_model};

    // Create an instance of unary force container (empty)
    unary_force_container_t
        unary_force_functors;

    // Create an instance of step_handler
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
        step_handler_instance;

    granular_system_t system(x0,
        v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
        step_handler_instance, binary_force_functors, unary_force_functors);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump " << n / dump_period << " out of " << n_dumps << std::endl;
            dump_particle_positions("run", n / dump_period, system.get_x(), r_part);
        }
        system.do_step(dt);
    }
    
    return 0;
}
```

If you run this simulation and use the steps described in the [driver program tutorial](#simulation-example)
to set up a visualization,
you will see a fractal aggregate that becomes more compact due to the linear attraction
between particles:

![Aggregate undergoing compaction](resources/images/aggregate.gif)

## Implementing custom unary force models

## Implementing custom step handlers

## Class reference

## Installation

If your project is a git repository, libgran can be added to it as a git submodule. For example, to install libgran and
its dependency libtimestep in a directory named deps under your project root, run:
```shell
git submodule add https://github.com/egor-demidov/libgran deps/libgran
git submodule add https://github.com/egor-demidov/libtimestep deps/libtimestep
```
Since libgran and libtimestep are header-only libraries, they can be included in your project without linking against
an object. Simply, add the include directories of libgran and libtimestep in your CMakeLists.txt file:
```cmake
include_directories(deps/libgran/include)
include_directories(deps/libtimestep/include)
```
By default, libgran uses C++ 17 parallel algorithms for computation of binary interactions. In case you would like to
use OpenMP instead, `LIBGRAN_USE_OMP` compile definition needs to be added. In your CMakeLists.txt, add the following
line:
```cmake
add_compile_definitions(LIBGRAN_USE_OMP)
```
Note that for improved performance / parallelization capabilities, additional compiler flags might be required on your
system. Platform- and toolchain-specific instructions are provided below. Also, the optimal number of threads depends on
the size of the granular system. Using too many threads can be detrimental to performance.

### Linux & GNU C++ compiler

The following compiler flags are recommended for best performance:
```
-O3 -flto=auto -march=native
```

#### Parallelization with C++ 17 algorithms

If you opted for C++ 17 parallel algorithms (default), Intel TBB library needs to be installed and linked against.
On Ubuntu, TBB can be installed with:
```shell
sudo apt install libtbb2-dev
```
Then, in your CMakeLists.txt file, add:
```cmake
find_package(TBB REQUIRED)

# Set up your targets...

target_link_libraries(<your target name> PRIVATE TBB::tbb)
```

#### Parallelization with OpenMP

If you opted for OpenMP, the following flag needs to be added:
```
-fopenmp
```
Then, every time you run your simulation, the number of threads can be set prior to execution with:
```shell
export OMP_NUM_THREADS=4
```

### Windows & MSVC compiler

The following compiler flags are recommended for best performance:
```
/O2 /EHsc /GL /fp:except
```
It is important that you **build configuration is set to Release**. Otherwise, MSVC produces
very slow binaries.

#### Parallelization with C++ 17 algorithms

C++ 17 parallel algorithms work out of the box under MSVC without the need for additional compiler
flags or libraries.

#### Parallelization with OpenMP

If you opted for OpenMP, the following flag needs to be added:
```
/openmp
```
Then, every time you run your simulation, the number of threads can be set prior to execution with:
```shell
set OMP_NUM_THREADS=4
```

### macOS & LLVM/clang compiler

LLVM/clang compiler has to be used with OpenMP parallelization.
LLVM/clang and OpenMP need to be installed:
```shell
brew install llvm libomp
```
To use LLVM/clang instead of AppleClang with CMake, run:
```shell
export CC=/usr/local/opt/llvm/bin/clang
export CXX=/usr/local/opt/llvm/bin/clang++
cmake <path to CMakeLists.txt>
```
The following compiler flags are recommended for best performance:
```
-O3 -fopenmp -flto -march=native
```
The following linker flags are recommended for best performance:
```
-lomp -flto
```
Every time you run your simulation, the number of threads can be set prior to execution with:
```shell
export OMP_NUM_THREADS=4
```

### macOS & AppleClang compiler

The default AppleClang compiler bundled with XCode **is not supported**, because it does not support
C++ 17 parallel algorithms or OpenMP.

## Acknowledgement

- Project funded by U.S. National Science Foundation, award #AGS-2222104

