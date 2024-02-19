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
attraction model and is designed to be easily extensible with custom models.

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

## Documentation

Comprehensive documentation, including tutorials and class reference,
[is available externally](https://egor-demidov.github.io/libgran).

## Acknowledgement

- Project funded by U.S. National Science Foundation, award #AGS-2222104

