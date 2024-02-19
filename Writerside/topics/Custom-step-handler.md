# Custom step handler

This tutorial wil teach you how to implement a custom step handler to be used with libgran.
A step handler is an object that takes position and velocity increments for each particle
as inputs from the integrator and updates the positions and accelerations of particles in
the granular system. So far, we have been using the built-in trivial step handler that
merely increments the positions and velocities. However, it may sometimes be desirable 
to override the default behavior. For example, if you need one particle in the simulation
to be fixed in space or undergo uniform motion,
you can define a custom step handler that will never increment the velocity
of that particle. This is what we will do in this tutorial. A custom step handler will
be created that takes as input the index of a particle that shall be fixed and does not
allow it to accelerate, while the other particles are allowed to accelerate normally.

## Before you start

Make sure you have the following:

- A cmake-based C++ project
- Proper compiler flags are set in the **CMakeLists.txt** file (see [installation](Installation.md))
- You have completed the [example simulation tutorial](Example-simulation.md)
- You have completed the [custom unary force model tutorial](Custom-unary-force-model.md)

## Defining a custom step handler

In this section, we will define a custom step handler in a header file inside your
project. This step handler will later be used by a driver program in a simulation.

1. Using the same environment as in the [example simulation tutorial](Example-simulation.md), create a header file named
   **fixed_particle_step_handler.h** and add include guards to it:

    ```C++
    #ifndef FIXED_PARTICLE_STEP_HANDLER
    #define FIXED_PARTICLE_STEP_HANDLER
    
    /* Code will go here */
    
    #endif //FIXED_PARTICLE_STEP_HANDLER
    ```
   
2. Now, between the include guards, add the template declaration:

    ```C++
    #include <memory>
   
    #include <libtimestep/step_handler/step_handler.h>
    
    template <
        typename field_container_t,
        typename field_value_t,
        template <
            typename _field_container_t,
            typename _field_value_t>
        typename base_step_handler_t>
    struct fixed_particle_step_handler {
    
    };
    ```
   
    The template parameter `field_contatiner_t` is the type of the container used to store positions
    and accelerations (`std::vector<Eigen::Vector3d>` in our case), `field_value_t` is the type of the primary field, and `base_step_handler_t` is the type of
    step handler that this step handler extends.

3. Now, implement the constructor. The constructor initializes takes as arguments the index of the particle that
   shall be fixed and the step handler instance that is being extended. Both these arguments are copied into member
   variables.

    ```C++
    fixed_particle_step_handler(
        size_t fixed_particle,
        base_step_handler_t<field_container_t, field_value_t> base_step_handler) :
        fixed_particle(fixed_particle), base_step_handler(std::move(base_step_handler)) {}
        
    private:
    const size_t fixed_particle;
    base_step_handler_t<field_container_t, field_value_t> base_step_handler;
    ```

    In your struct, after the constructor and before the `private` keyword, add:
    
    ```C++
    void increment_v(size_t n,                                                      // index of the value to increment
                     field_value_t const & dv,                                      // value of the velocity increment
                     typename field_container_t::const_iterator x_begin_itr,        // iterator pointing to the start of the x buffer
                     typename field_container_t::iterator v_begin_itr,              // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr,    // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr,    // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer
    
        // If this is the fixed particle, do not increment velocity
        if (n == fixed_particle)
            return;
    
        // Otherwise, let the base step handler handle in
        base_step_handler.increment_v(n, dv, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }
    ```
   
    This is the function that will be called by the integrator to increment velocities. Inside the function,
    we check if the index of the target particle matches the index of the particle that should be fixed. If it does,
    we do not increment the velocity. Otherwise, we let the base handler handle the increment.

4. For position, orientation, and angular velocity, simply pass on the arguments to the base step handler:

    ```C++
    void increment_x(size_t n,                                                      // index of the value to increment
                     field_value_t const & dx,                                      // value of the position increment
                     typename field_container_t::iterator x_begin_itr,              // iterator pointing to the start of the x buffer
                     typename field_container_t::const_iterator v_begin_itr,        // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr,    // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr,    // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer
    
        base_step_handler.increment_x(n, dx, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }
    ```
    
    ```C++
    void increment_theta(size_t n,                                                  // index of the value to increment
                     field_value_t const & dtheta,                                  // value of the orientation increment
                     typename field_container_t::const_iterator x_begin_itr,        // iterator pointing to the start of the x buffer
                     typename field_container_t::const_iterator v_begin_itr,        // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                     typename field_container_t::iterator theta_begin_itr,          // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr,    // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer
    
        base_step_handler.increment_theta(n, dtheta, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }
    ```
    
    ```C++
    void increment_omega(size_t n,                                                  // index of the value to increment
                     field_value_t const & domega,                                  // value of the angular velocity increment
                     typename field_container_t::const_iterator x_begin_itr,        // iterator pointing to the start of the x buffer
                     typename field_container_t::const_iterator v_begin_itr,        // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr,    // iterator pointing to the start of the theta buffer
                     typename field_container_t::iterator omega_begin_itr,          // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer
    
        base_step_handler.increment_omega(n, domega, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }
    ```
   
## Driver program for this model

In this section of the tutorial, we provide a driver program that can be sued to test the custom 
step handler defined in the previous section. For mode details about driver programs, refer to the
[example simulation tutorial](Example-simulation.md). We will be using the custom unary force model
created in the [custom unary force model](Custom-unary-force-model.md) tutorial, **central_attraction.h**, and the initial
geometry, **aggregate.h**. Either copy them into this project or create a new driver program in the same project.
The simulation setup will be similar to what was done in the [custom unary force model](Custom-unary-force-model.md)
tutorial, except one of the particles in the system will be fixed in space.

1. Copy the driver program source code, **tutorial_04.cpp**, as defined below, into your project.
   The program is similar to the one in the previous tutorial, except we also create an instance of
   `fixed_particle_step_handler` and pass it to the granular system constructor. We use the default step
   handler as the base handler for our custom handler. We also increase the Hamaker constant, `A`, by a
   factor of 10 to make adhesion between neighboring monomers stronger (so the aggregate does not fall apart immediately).

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
    #include "fixed_particle_step_handler.h" // Include your custom step handler
    
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
    
    // Create an alias for your custom step handler using rotational_step_handler as the base
    template <typename field_container_t, typename field_value_t>
    using step_handler_t = fixed_particle_step_handler<field_container_t, field_value_t, rotational_step_handler>;
    
    // Make sure that granular system is nofigured to use the new step handler
    using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
        step_handler_t, binary_force_container_t, unary_force_container_t>; // Granular system representation
    
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
        const double A = 1.0e-19;
        const double h0 = 1.0e-9;
    
        // Parameter for your custom model
        const double k_attr = 8.0e-6;
    
        // Parameter for your custom step handler
        const size_t fixed_particle = 10;
    
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
    
        // Create an instance of base step handler
        rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            base_step_handler_instance;
    
        step_handler_t<std::vector<Eigen::Vector3d>, Eigen::Vector3d>
            step_handler_instance(fixed_particle, base_step_handler_instance);
    
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
    {collapsible="true" collapsed-title="tutorial_04.cpp"}

2. Create a **run** directory in the same directory as the executable. Run the simulation.
   After the **run** directory has been populated with 300 **.csv** files,
   follow the steps from the [example simulation tutorial](Example-simulation.md) to set up
   a visualization in ParaView. You will see a fractal aggregate undergoing compaction do to a central
   force acting on each primary particle. One of the particles near the top of the aggregate is held stationary.
   That particle resists compaction and the aggregate ends up undergoing a complicated motion
   about the fixed particle.

   ![Aggregate undergoing compaction due to a central force](aggregate_fixed.gif)

## What you've learned

You have leaned how to:

- Implement a custom step handler
- Use the custom step handler in a simulation

## Complete source file

```C++
#ifndef FIXED_PARTICLE_STEP_HANDLER_H
#define FIXED_PARTICLE_STEP_HANDLER_H

#include <memory>

#include <libtimestep/step_handler/step_handler.h>

template <
    typename field_container_t,
    typename field_value_t,
    template <
        typename _field_container_t,
        typename _field_value_t>
    typename base_step_handler_t>
struct fixed_particle_step_handler {
    fixed_particle_step_handler(
    size_t fixed_particle,
    base_step_handler_t<field_container_t, field_value_t> base_step_handler) :
    fixed_particle(fixed_particle), base_step_handler(std::move(base_step_handler)) {}

    void increment_v(size_t n,                                                  // index of the value to increment
                 field_value_t const & dv,                                      // value of the velocity increment
                 typename field_container_t::const_iterator x_begin_itr,        // iterator pointing to the start of the x buffer
                 typename field_container_t::iterator v_begin_itr,              // iterator pointing to the start of the v buffer
                 typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                 typename field_container_t::const_iterator theta_begin_itr,    // iterator pointing to the start of the theta buffer
                 typename field_container_t::const_iterator omega_begin_itr,    // iterator pointing to the start of the omega buffer
                 typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer

        // If this is the fixed particle, do not increment velocity
        if (n == fixed_particle)
            return;

        // Otherwise, let the base step handler handle in
        base_step_handler.increment_v(n, dv, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    void increment_x(size_t n,                                                  // index of the value to increment
                 field_value_t const & dx,                                      // value of the position increment
                 typename field_container_t::iterator x_begin_itr,              // iterator pointing to the start of the x buffer
                 typename field_container_t::const_iterator v_begin_itr,        // iterator pointing to the start of the v buffer
                 typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                 typename field_container_t::const_iterator theta_begin_itr,    // iterator pointing to the start of the theta buffer
                 typename field_container_t::const_iterator omega_begin_itr,    // iterator pointing to the start of the omega buffer
                 typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_x(n, dx, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    void increment_theta(size_t n,                                              // index of the value to increment
                 field_value_t const & dtheta,                                  // value of the orientation increment
                 typename field_container_t::const_iterator x_begin_itr,        // iterator pointing to the start of the x buffer
                 typename field_container_t::const_iterator v_begin_itr,        // iterator pointing to the start of the v buffer
                 typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                 typename field_container_t::iterator theta_begin_itr,          // iterator pointing to the start of the theta buffer
                 typename field_container_t::const_iterator omega_begin_itr,    // iterator pointing to the start of the omega buffer
                 typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_theta(n, dtheta, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    void increment_omega(size_t n,                                              // index of the value to increment
                 field_value_t const & domega,                                  // value of the angular velocity increment
                 typename field_container_t::const_iterator x_begin_itr,        // iterator pointing to the start of the x buffer
                 typename field_container_t::const_iterator v_begin_itr,        // iterator pointing to the start of the v buffer
                 typename field_container_t::const_iterator a_begin_itr,        // iterator pointing to the start of the a buffer
                 typename field_container_t::const_iterator theta_begin_itr,    // iterator pointing to the start of the theta buffer
                 typename field_container_t::iterator omega_begin_itr,          // iterator pointing to the start of the omega buffer
                 typename field_container_t::const_iterator alpha_begin_itr) {  // iterator pointing to the start of the alpha buffer

    base_step_handler.increment_omega(n, domega, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
}

private:
    const size_t fixed_particle;
    base_step_handler_t<field_container_t, field_value_t> base_step_handler;
};

#endif //FIXED_PARTICLE_STEP_HANDLER_H
```
{collapsible="true" collapsed-title="fixed_particle_step_handler.h"}

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
