//
// Created by egor on 2/22/24.
//

#ifndef TRIANGULAR_FACET_H
#define TRIANGULAR_FACET_H

// triangular_facet is a unary force model that is a container
// for surface unary force models
template <
    typename field_value_t,
    typename real_t,
    typename... surface_force_functors_t>
struct triangular_facet {

    typedef std::vector<field_value_t> field_container_t;

    explicit triangular_facet(field_value_t v_facet, std::tuple<field_value_t, field_value_t, field_value_t> vertices,
        surface_force_functors_t & ... force_functors) :
            v_facet(std::move(v_facet)), vertices(std::move(vertices)), surface_force_functors{force_functors...} {

        // Compute the area of the triangle
        field_value_t u = std::get<1>(vertices) - std::get<0>(vertices);
        field_value_t v = std::get<2>(vertices) - std::get<0>(vertices);

        facet_area = 1.0 / 2.0 * u.cross(v).norm();

        // Initiazlize the unit normal
        plane_unit_normal = u.cross(v).normalized();
    }

    std::pair<field_value_t, field_value_t> operator () (size_t i, field_container_t const & x,
                                                                  field_container_t const & v,
                                                                  field_container_t const & theta,
                                                                  field_container_t const & omega,
                                                                  real_t t) {
        static_assert(std::tuple_size<decltype(surface_force_functors)>::value > 0, "at least one surface force functor must be provided "
                                                                            "as a template parameter");

        // Compute the nearest point on the facet to the center of the particle
        field_value_t x_facet = closest_point_on_triangle(x[i]);

        // Create a copy of facet velocity that can be captured by lambdas
        field_value_t v_facet_copy = v_facet;

        // Call the surface force functors passing the x_facet to each of them as an argument
        std::pair<field_value_t, field_value_t> accelerations = std::apply([i, &x_facet, &v_facet_copy, &x, &v, &theta, &omega, t] (auto & ... e) -> std::pair<field_value_t, field_value_t> {
            return (e(i, x_facet, v_facet_copy, x, v, theta, omega, t) + ...);
        }, surface_force_functors);
        return accelerations;
    }

    field_value_t const & get_unit_normal() const {
        return plane_unit_normal;
    }

    // Velocity of the facet
    field_value_t v_facet;

private:

    field_value_t closest_point_on_plane(field_value_t const & test_point) const {
        field_value_t offstet_from_origin = std::get<0>(vertices).dot(plane_unit_normal) * plane_unit_normal;
        field_value_t projection = test_point - test_point.dot(plane_unit_normal) * plane_unit_normal;

        return projection + offstet_from_origin;
    }

    bool is_point_on_triangle(field_value_t const & test_point) const {
        real_t area = 0.0;
        field_value_t pa = std::get<0>(vertices) - test_point;
        field_value_t pb = std::get<1>(vertices) - test_point;
        field_value_t pc = std::get<2>(vertices) - test_point;
        area += 1.0 / 2.0 * pa.cross(pb).norm();
        area += 1.0 / 2.0 * pb.cross(pc).norm();
        area += 1.0 / 2.0 * pc.cross(pa).norm();

        return !(area > facet_area);
    }

    template <size_t i, size_t j>  // Vertex indices to define the edge
    field_value_t closest_point_on_edge(field_value_t const & test_point) const {
        field_value_t line_vector = std::get<j>(vertices) - std::get<i>(vertices);
        field_value_t unit_line_vector = line_vector.normalized();
        field_value_t test_point_vector = test_point - std::get<i>(vertices);
        real_t line_length = line_vector.norm();
        real_t test_point_vector_projection = test_point_vector.dot(unit_line_vector);

        if (test_point_vector_projection < 0.0)
            return std::get<i>(vertices);

        if (test_point_vector_projection > line_length)
            return std::get<j>(vertices);

        return std::get<i>(vertices) + test_point_vector_projection * unit_line_vector;
    }

    field_value_t closest_point_on_triangle(field_value_t const & test_point) const {
        field_value_t closest_point = closest_point_on_plane(test_point);

        if (is_point_on_triangle(closest_point))
            return closest_point;

        field_value_t closest_point_edge_1 = closest_point_on_edge<0, 1>(test_point);
        field_value_t closest_point_edge_2 = closest_point_on_edge<1, 2>(test_point);
        field_value_t closest_point_edge_3 = closest_point_on_edge<2, 0>(test_point);

        field_value_t dist1 = test_point - closest_point_edge_1;
        field_value_t dist2 = test_point - closest_point_edge_2;
        field_value_t dist3 = test_point - closest_point_edge_3;

        double dist1_squared = dist1.dot(dist1);
        double dist2_squared = dist2.dot(dist2);
        double dist3_squared = dist3.dot(dist3);

        if (dist1_squared < dist2_squared) {
            if (dist1_squared < dist3_squared) {
                return closest_point_edge_1;
            }
        } else {
            if (dist2_squared < dist3_squared) {
                return closest_point_edge_2;
            }
        }
        return closest_point_edge_3;
    }

    real_t facet_area;
    std::tuple<field_value_t, field_value_t, field_value_t> vertices;
    std::tuple<surface_force_functors_t & ...> surface_force_functors;
    field_value_t plane_unit_normal;
};

#endif //TRIANGULAR_FACET_H
