# Custom unary force model

This tutorial will teach you how to create a custom unary force model to be used with libgran.
A unary force is a force acting on a particle independently from other particles. In this tutorial,
we will create a unary linear central force. A particle will be attracted to the center of mass
of the aggregate as:

```tex
\mathbf{f}=k_{\rm a}\left(\mathbf{x}_i-\mathbf{x}_0\right)
```

where **f** is force, _k<sub>a</sub>_ is the proportionality constant, **x**<sub>i</sub> is the position
of the particle, and **x**<sub>0</sub> is the center of mass of the aggregate.

## Before you start

Make sure you have the following:

- A cmake-based C++ project
- Proper compiler flags are set in the **CMakeLists.txt** file (see [installation](Installation.md))
- You have completed the [example simulation tutorial](Example-simulation.md)

## Defining a custom unary force

In this section, we will define a custom unary force model in a header file inside your project.
This model will later be used by a driver program in a simulation.

1. Using the same environment as in the [example simulation tutorial](Example-simulation.md), create a header file named
   **central_attraction.h** and add include guards to it:

    ```C++
    #ifndef CENTRAL_ATTRACTION_H
    #define CENTRAL_ATTRACTION_H
    
    /* Code will go here */
    
    #endif //CENTRAL_ATTRACTION_H
    ```

2. Now, between the include guards, add the template declaration:

    ```C++
    template<typename field_value_t, typename real_t>
    struct central_attraction_functor {
        
    };
    ```

   The template parameter `field_value_t` is the type of the primary field (like `Eigen::Vector3d`
   or `Eigen::Vector3f`) and the parameter `real_t` is the type of a real-valued scalar value (like `double`
   or `float`).

3. Now, we will implement a constructor. In order compute accelerations, we need to know
   _r_, _k<sub>a</sub>_, and _m_. We will also need to know the zero value of the primary field.
   Inside the struct, declare the private member constants and the
   constructor. The constructor only initializes the member constants:

    ```C++
   central_attraction_functor(real_t r, real_t k_a, real_t mass, field_value_t zero_field) :
        r(r), k_a(k_a), mass(mass), zero_field(zero_field) {}

    private:
    const real_t r, k_a, mass;
    const field_value_t zero_field;
    field_value_t x0;
    ```

   In your struct, after the constructor and before the `private` keyword, add:

   ```c++
   std::pair<field_value_t, field_value_t> operator () 
      (size_t i,
       std::vector<field_value_t> const & x,
       std::vector<field_value_t> const & v [[maybe_unused]],
       std::vector<field_value_t> const & theta [[maybe_unused]],
       std::vector<field_value_t> const & omega [[maybe_unused]],
       real_t t [[maybe_unused]]) const {
       
       /* Compute and return the accelerations */
   }
   ```

   This lengthy statement is the signature of `operator ()` - the function that will be called by the
   integrator to compute the acceleration acing on particle i due to its linear attraction to the center of
   mass of the granular system. The arguments are: index of particle i,, position buffer,
   velocity buffer, orientation buffer, angular velocity buffer, and time. All these arguments must be accepted
   by the `operator ()` function even if they are not used. In case an argument is not used, it should be
   marked with the `[[maybe_unused]]` attribute to avoid causing compiler warnings.

4. Acceleration can be computed with the following equation:

    ```tex
    \mathbf{a}=\frac{k_{\rm a}}{m}\left(\mathbf{x}_i-\mathbf{x}_0\right)
    ```
   
    To implement this computation, add inside your `operator ()`:

    ```C++
    field_value_t a = k_a / mass * (x[i] - x0); // Compute the central acceleration
    return std::make_pair(a, zero_field); // Return the translational and angular accelerations
    ```

   > Note that the return type
   > is an `std::pair`, because we need to return two accelerations: translational and angular. However, we assume that
   > this central attractive model is a body force acting on the particles and does not directly produce torques. Therefore,
   > the angular acceleration that we return is be zero.

5. At every time step, the center of mass needs to be re-computed as it may shift with the evolution of the system.
   We do not want to compute the center of mass directly in the body of `operator ()` because that computation would be
   done for every particle. We only want the center of mass to be re-computed once per time step.
   We will define an `update_center_of_mass()` method, which will need to be called by the driver program at every time step. 
   Inside your struct before the `private` keyword, add:

    ```C++
    void update_center_of_mass(std::vector<field_value_t> const & x) {
        field_value_t x0_new = zero_field;
        for (auto const & point : x) {
            x0_new += point;
        }
        x0 = x0_new / real_t(x.size());
    }
    ```
    
    Since we want to ensure that the central force model instance that we store in the driver program is the 
    same instance that is being used in the granular system object
    (otherwise updating the center of mass will have no effect on the simulation), it is a good idea to `delete` the copy
    constructor of `central_attraction_functor` to prevent unintentional copies. In your struct, add:

    ```C++
    central_attraction_functor(central_attraction_functor const &) = delete;
    ```

## Driver program for this model

In this section of the tutorial, we provide a driver program that can be used to test the custom force model defined in the
previous section. For more details about driver programs, refer to the [example simulation tutorial](Example-simulation.md).

1. In this simulation, we will be using a fractal aggregate as the initial geometry for the simulation.
   Copy the **aggregate.h** file, as defined below, into your project. This file contains initial coordinates of particles
   normalized by particle radius and will be included in the driver program to set up the initial position buffer.

   ```C++
   #ifndef AGGREGATE_H
   #define AGGREGATE_H
   
   constexpr double aggregate[][3] = {
         5.381048,   18.483654,    5.143242,
         5.949202,   16.941475,    6.282927,
         6.522552,   19.743778,    4.090125,
         4.254210,   20.535814,    4.447582,
         3.061821,   20.474491,    2.843072,
         4.059029,   19.564165,    6.184764,
         1.961491,   22.027952,    3.456290,
         4.100143,   16.437365,    6.854641,
         8.205635,   20.156048,    3.091492,
         8.266253,   20.484835,    5.063351,
         9.585005,   19.294472,    1.927438,
         9.182722,   21.896002,    2.957845,
        10.888521,   21.026497,    3.535957,
         0.298204,   13.117410,   -3.956084,
        -1.850106,   10.994704,   -0.392680,
        -0.157987,   10.679225,    0.625754,
        -2.152361,   10.906271,   -2.367730,
        -0.440803,   11.384176,   -3.285434,
        -2.186500,   11.783336,   -4.176061,
         1.443839,   14.256097,   -2.776716,
         0.293551,   12.450349,    1.437687,
        -2.130326,   16.880444,    2.636410,
        -0.960546,   15.313308,    3.055574,
         0.426633,   16.545813,    3.801675,
         1.009729,   16.275637,    5.695614,
         2.422916,   17.210443,    4.633050,
         2.095561,   18.965055,    3.730735,
         1.240510,   20.379923,    4.856364,
         2.820676,   15.253491,    4.522942,
         3.315668,   14.109584,    6.087059,
         2.162122,   12.846807,    7.123750,
        -1.178428,   13.457395,    2.342753,
        -1.954622,   14.948363,    1.258995,
        -2.105695,   11.862226,    3.114516,
        -3.326299,   14.136928,    0.050659,
        -4.645139,   11.413218,    2.226980,
        -6.326225,   12.281743,    1.579208,
        -7.691842,   11.562547,    0.307262,
        -3.592058,   11.337117,    0.528383,
         3.433121,    1.539043,   -1.650279,
         3.214089,   -0.063356,   -2.826862,
         2.581532,    2.291854,   -3.295901,
         1.107985,    3.338170,   -2.439232,
        -0.651946,    2.841426,   -3.249099,
         3.702398,    2.006700,   -7.527482,
         2.479157,    4.073030,   -4.199743,
         3.837380,    5.696042,   -6.078282,
         2.305369,    4.413163,   -6.162930,
         2.603780,    3.617363,   -7.973360,
         2.595674,    2.119801,   -9.298969,
         4.004651,   -0.764775,    6.374567,
         4.431592,   -1.264468,    4.485642,
         4.255418,    0.571209,    3.711526,
         3.228247,    1.314791,    2.164913,
         1.164471,    0.713804,    4.137331,
         1.859790,    2.339574,    3.202768,
        -0.435625,   -0.397126,    3.683978,
         3.454807,    1.993493,    0.297285,
        -3.619100,    1.481467,   -1.152247,
        -2.618960,    1.530801,    0.579020,
        -0.225961,   -1.293818,    0.275802,
        -2.165546,   -1.756202,    0.120204,
         0.382767,    0.640474,    2.182770,
        -0.865854,    0.569151,    0.622046,
         2.572918,   -7.939448,    1.429085,
         1.021404,   -8.169539,    2.669998,
         4.336183,   -8.941447,    2.506176,
         4.491859,   -7.414514,    1.223896,
         9.602220,   -7.352553,    2.891611,
         6.242034,   -9.460876,    2.819118,
         7.920323,   -8.041297,    4.601250,
         7.872787,   -8.320031,    2.621339,
         7.614533,  -11.089431,    5.461537,
         6.889518,  -11.184531,    3.600003,
         5.142938,  -11.849626,    4.312113,
         8.887940,   -5.485045,    2.938691,
         7.860165,   -3.887006,    4.323622,
         6.701442,   -2.581828,    2.096128,
         8.379007,   -3.624541,    2.410009,
         5.009143,   -2.646241,    3.160069,
         0.546897,   -7.253179,   -0.813710,
         2.084806,   -6.452823,    0.183420,
        11.189775,   -3.693578,    4.221244,
         9.665753,   -2.591838,    3.540412,
         0.240105,   -8.962543,    0.575863,
         1.044943,   -9.141250,   -1.246308,
        -8.864514,    6.251542,    1.875130,
        -8.735852,    8.157282,    1.282172,
        -9.732181,    4.710372,    0.941343,
        -5.316808,    5.346671,    1.328840,
        -7.201322,    5.817655,    0.852616,
        -6.969355,    8.757685,    2.002597,
        -5.052598,    9.232308,    1.685147,
       -10.212317,   10.201171,   -0.612693,
        -8.928578,   10.114470,    0.918480,
        -1.060849,    5.208689,    1.357853,
         0.308861,    5.613111,    2.757973,
        -1.075000,    5.751388,   -0.567057,
        -0.253025,    4.021375,   -1.142740,
        -3.377661,    5.562436,    3.342318,
        -2.314340,    4.112841,    2.465950,
        -5.970988,    5.047156,   -1.298219,
        -5.273189,    3.628005,   -0.073847,
        -4.030825,    3.822539,    1.481367,
        -3.126646,  -14.718781,   -6.370821,
        -3.336624,  -15.911931,   -4.779500,
        -4.828164,  -12.135728,   -5.731589,
        -3.721381,  -11.820348,   -7.367306,
        -4.545673,  -11.434608,   -3.879934,
        -2.197195,  -13.032529,   -6.911812,
        -0.998573,   -9.687873,   -5.965328,
        -1.987809,   -8.058442,   -5.360044,
        -2.164439,   -7.910609,   -8.875383,
        -1.257186,   -7.508889,   -7.138861,
        -1.272729,  -10.431705,   -7.801506,
        -0.543833,  -12.293704,   -7.760655,
        -1.640402,   -7.176370,  -10.660398,
        -4.383595,   -4.705949,   -3.012371,
        -3.663079,   -6.502570,   -2.509365,
        -2.311150,   -7.805171,   -1.452906,
        -1.256305,   -6.845490,   -0.050651,
        -2.250354,   -7.570956,   -3.438215,
        -5.456369,   -5.816705,   -4.283346,
        -4.878307,  -15.496701,   -6.942243,
        -0.552877,  -17.577950,    0.008364,
         0.189493,  -16.927700,    1.747923,
         0.937788,  -14.188067,   -1.174223,
         0.861116,  -12.785125,   -2.597557,
         1.244249,  -14.455119,    0.784032,
         2.833199,  -13.464568,    0.081148,
        -0.221441,  -14.991487,    2.034657,
        -1.131078,  -12.925971,   -2.491109,
         4.009218,  -14.982114,   -0.479234,
         5.564094,  -14.644542,   -1.691013,
        -3.145387,  -13.399182,    0.264053,
        -1.355115,  -12.946044,   -0.503797,
        -4.013848,  -14.934554,    1.206608,
        -2.738651,  -15.922336,    0.024169,
        -6.909004,  -15.684735,   -0.757439,
        -7.107459,  -14.953906,   -2.608522,
        -6.167473,  -13.942424,   -0.113636,
        -5.983214,  -15.205661,    1.425936,
        -5.263436,  -15.812198,   -4.252769,
        -7.157816,  -16.397881,   -3.991421,
        -6.993089,  -20.749481,   -3.291241,
        -6.177739,  -20.545391,   -5.106055,
        -8.282816,  -19.537605,   -2.359597,
        -4.882554,  -17.485018,   -3.224855,
        -4.975591,  -19.141996,   -4.340994,
        -6.324211,  -19.777241,   -6.946840
   };
   
   #endif //AGGREGATE_H
   ```
   {collapsible="true" collapsed-title="aggregate.h"}

2. Copy the driver program source code, **tutorial_03.cpp**, as defined below, into your project. In this program,
   we use the aggregate file to initialize the geometry and run the simulation with three force models:
   frictional contact, Van der Waals attraction, and central attraction.

   ```C++
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
   
   #include "central_attraction.h" // Include your custom model
   
   // Include the fractal aggregate used as initial conditions
   #include "aggregate.h"
   
   using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>; // Contact force
   using vdw_force_functor_t = hamaker_functor<Eigen::Vector3d, double>; // Van der Waals force
   using central_force_functor_t = central_attraction_functor<Eigen::Vector3d, double>; // Your custom model
   using binary_force_container_t =
       binary_force_functor_container<Eigen::Vector3d, double,
       contact_force_functor_t, vdw_force_functor_t>; // Binary force container
   
   using unary_force_container_t =
       unary_force_functor_container<Eigen::Vector3d, double, central_force_functor_t>; // Unary force container (empty)
   
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
       const double k_attr = 8.0e-6;
   
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
       central_force_functor_t central_model(r_part, k_attr, mass, Eigen::Vector3d::Zero());
   
       // Create an instance of binary force container
       binary_force_container_t
           binary_force_functors{contact_force_model, hamaker_model};
   
       // Create an instance of unary force container (empty)
       unary_force_container_t
           unary_force_functors{central_model};
   
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
           central_model.update_center_of_mass(x0);
           system.do_step(dt);
       }
   
       return 0;
   }
   ```
   {collapsible="true" collapsed-title="tutorial_03.cpp"}

3. Create a directory named **run** in the same folder as the executable and run the simulation. Set up a visualization,
   as shown in the [example simulation tutorial](Example-simulation.md). You should see an initially fractal aggregate
   that undergoes compaction due to the central attractive force:

   ![Fractal aggregate undergoing compaction](aggregate_central.gif)

## What you've learned

You have learned how to:

- Implement a custom force model for independent forces acting on particles
- Use the custom model in combination with built-in force models

## Complete source file

```C++
#ifndef CENTRAL_ATTRACTION_H
#define CENTRAL_ATTRACTION_H

template<typename field_value_t, typename real_t>
struct central_attraction_functor {
    central_attraction_functor(central_attraction_functor const &) = delete;

    central_attraction_functor(real_t r, real_t k_a, real_t mass, field_value_t zero_field) :
    r(r), k_a(k_a), mass(mass), zero_field(zero_field) {}

    std::pair<field_value_t, field_value_t> operator ()
       (size_t i,
        std::vector<field_value_t> const & x,
        std::vector<field_value_t> const & v [[maybe_unused]],
        std::vector<field_value_t> const & theta [[maybe_unused]],
        std::vector<field_value_t> const & omega [[maybe_unused]],
        real_t t [[maybe_unused]]) const {

        field_value_t a = k_a / mass * (x0 - x[i]); // Compute the central acceleration
        return std::make_pair(a, zero_field); // Return the translational and angular accelerations
    }

    void update_center_of_mass(std::vector<field_value_t> const & x) {
        field_value_t x0_new = zero_field;
        for (auto const & point : x) {
            x0_new += point;
        }
        x0 = x0_new / real_t(x.size());
    }

private:
    const real_t r, k_a, mass;
    const field_value_t zero_field;
    field_value_t x0;
};

#endif //CENTRAL_ATTRACTION_H
```
{collapsible="true" collapsed-title="central_attraction.h"}