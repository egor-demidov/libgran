# Overview

libgran is a Discrete Element Method (DEM) framework for simulating the mechanical behavior of soot aggregates. DEM is a
technique for simulation of granular media consisting of rigid spherical particles. The resultant force and torque acting
on each particle are computed and used with Newton's second law to compute the motion of particles:
```tex
m\ddot{\mathbf{x}}=\mathbf{f}
```
```tex
I\ddot{\omega}=\bm{\tau}
```
The forces that particles experience arise from friction at inter-particle contacts, bonding between particles,
inter-particle attraction, field forces, etc. libgran contains a bonded and a non-bonded contact model, a Van der Waals
attraction model and is designed to be easily extensible with custom models.