#ifndef NUMERICAL_H
#define NUMERICAL_H

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef NUMDEF
#ifdef NUMERICAL_STATIC
#define NUMDEF static
#else
#define NUMDEF extern
#endif
#endif

    typedef struct
    {
        double x;
        double y;
        double z;
    } vec3_t;

    NUMDEF bool assert_eq(double a, double b);

    NUMDEF vec3_t vec3_subtract(vec3_t *origin, vec3_t *end);
    NUMDEF double max(double a, double b);
    NUMDEF double min(double a, double b);
    NUMDEF double vec3_dot(vec3_t *a, vec3_t *b);
    NUMDEF vec3_t vec3_cross(vec3_t *a, vec3_t *b);
    NUMDEF double vec3_l2_norm(vec3_t *a);
    NUMDEF void vec3_normalise_in_place(vec3_t *a);
    NUMDEF void vec3_invert_in_place(vec3_t *a);
    NUMDEF void vec3_coord_trans_xz_in_place(vec3_t *a, double radian);
    NUMDEF void vec3_coord_trans_xy_in_place(vec3_t *a, double radian);

    NUMDEF bool is_simple_polygon(vec3_t *vertices, size_t nvtx);

#ifdef __cplusplus
}
#endif

#endif // NUMERICAL_H

#ifdef NUMERICAL_IMPLEMENTATION

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#define NUM_EPSILON 1.0e-9

bool assert_eq(double a, double b)
{
    return fabs(a - b) < NUM_EPSILON;
}

double max(double a, double b)
{
    return (a > b) ? a : b;
}

double min(double a, double b)
{
    return (a < b) ? a : b;
}

NUMDEF vec3_t vec3_subtract(vec3_t *origin, vec3_t *end)
{
    vec3_t a = {
        .x = end->x - origin->x,
        .y = end->y - origin->y,
        .z = end->z - origin->z,
    };

    return a;
}

NUMDEF double vec3_dot(vec3_t *a, vec3_t *b)
{
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

NUMDEF vec3_t vec3_cross(vec3_t *a, vec3_t *b)
{
    vec3_t prod = {
        .x = a->y * b->z - a->z * b->y,
        .y = a->z * b->x - a->x * b->z,
        .z = a->x * b->y - a->y * b->x,
    };

    return prod;
}

NUMDEF double vec3_l2_norm(vec3_t *a)
{
    return sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
}

NUMDEF void vec3_normalise_in_place(vec3_t *a)
{
    double norm = vec3_l2_norm(a);

    a->x /= norm;
    a->y /= norm;
    a->z /= norm;
}

NUMDEF void vec3_invert_in_place(vec3_t *a)
{
    a->x *= -1.0;
    a->y *= -1.0;
    a->z *= -1.0;
}

NUMDEF void vec3_coord_trans_xz_in_place(vec3_t *a, double radian)
{
    double x = a->x;
    double z = a->z;

    a->x = x * cos(radian) - z * sin(radian);
    a->z = x * sin(radian) + z * cos(radian);
}

NUMDEF void vec3_coord_trans_xy_in_place(vec3_t *a, double radian)
{
    double x = a->x;
    double y = a->y;

    a->x = x * cos(radian) + y * sin(radian);
    a->y = -1.0 * x * sin(radian) + y * cos(radian);
}

static bool is_coplanar(vec3_t *vertices, size_t nvtx)
{
    // define new origin p0
    vec3_t p0 = vertices[0];

    // find first distinct vertex p1
    vec3_t p1 = {0};
    size_t i;
    bool found_p1 = false;

    for (i = 1; i < nvtx; ++i)
    {
        vec3_t diff = vec3_subtract(&vertices[i], &p0);

        if (vec3_l2_norm(&diff) > NUM_EPSILON)
        {
            p1 = vertices[i];
            found_p1 = true;
            break;
        }
    }

    if (!found_p1)
    {
        fprintf(stderr, "%-8sAll vertices are identical\n", "[ERROR]");
        return false;
    }

    // find a third vertex p2 that is non-collinear with p0 and p1
    vec3_t p2 = {0};
    vec3_t u = vec3_subtract(&p1, &p0);
    bool found_p2 = false;

    for (size_t j = i + 1; j < nvtx; ++j)
    {
        vec3_t v_tmp = vec3_subtract(&vertices[j], &p0);
        vec3_t cross = vec3_cross(&u, &v_tmp);

        if (vec3_l2_norm(&cross) > NUM_EPSILON)
        {
            p2 = vertices[j];
            found_p2 = true;
            break;
        }
    }

    if (!found_p2)
    {
        fprintf(stderr, "%-8sAll vertices are collinear; please use non-collinear vertices to define a plane\n", "[ERROR]");
        return false;
    }

    // determine the plane normal using p0, p1, and p2
    vec3_t v = vec3_subtract(&p2, &p0);
    vec3_t normal = vec3_cross(&u, &v);

    // check coplanarity
    for (size_t k = 0; k < nvtx; ++k)
    {
        vec3_t diff = vec3_subtract(&vertices[k], &p0);
        double dist = vec3_dot(&diff, &normal);

        if (fabs(dist) > NUM_EPSILON)
        {
            fprintf(stderr, "%-8sVertex %zu deviates from the plane\n", "[ERROR]", k);
            return false;
        }
    }

    // coplanarity test passed
    return true;
}

static void vertices_proj_to_2d_in_place(vec3_t *vertices, size_t nvtx)
{
    // select first vertex as origin
    vec3_t p0 = vertices[0];

    // create a basis vector u
    vec3_t u = vec3_subtract(&vertices[1], &p0);
    vec3_normalise_in_place(&u);

    // determine plane normal, using first three non-collinear vertices
    // (assuming all vertices are non-coplanar)
    vec3_t normal = {0};

    for (size_t i = 2; i < nvtx; ++i)
    {
        vec3_t v_cand = vec3_subtract(&vertices[i], &p0);
        vec3_t prod = vec3_cross(&u, &v_cand);

        if (vec3_l2_norm(&prod) > NUM_EPSILON)
        {
            normal = prod;
            break;
        }
    }

    // define second basis vector v, perpendicular to u
    vec3_t v = vec3_cross(&normal, &u);
    vec3_normalise_in_place(&v);

    // for each vertex, mutate it to corresponding 2D projection
    for (size_t i = 0; i < nvtx; ++i)
    {
        vec3_t diff = vec3_subtract(&vertices[i], &p0);
        vertices[i].x = vec3_dot(&diff, &u);
        vertices[i].y = vec3_dot(&diff, &v);
        vertices[i].z = 0.0;
    }
}

static int vertices_orientation(vec3_t *a, vec3_t *b, vec3_t *c)
{
    double res = (b->y - a->y) * (c->x - b->x) - (b->x - a->x) * (c->y - b->y);

    // return 0 if collinear, 1 if clockwise, -1 if counterclockwise
    if (fabs(res) < NUM_EPSILON)
        return 0;

    return (res > 0) ? 1 : -1;
}

static bool vertex_on_segment(vec3_t *a, vec3_t *b, vec3_t *c)
{
    // check if b lies on segment ac (assuming collinearity)
    return (b->x >= min(a->x, c->x) && b->x <= max(a->x, c->x) &&
            b->y >= min(a->y, c->y) && b->y <= max(a->y, c->y));
}

static bool is_intersect(vec3_t *p1, vec3_t *p2, vec3_t *q1, vec3_t *q2)
{
    int o1 = vertices_orientation(p1, p2, q1);
    int o2 = vertices_orientation(p1, p2, q2);
    int o3 = vertices_orientation(q1, q2, p1);
    int o4 = vertices_orientation(q1, q2, p2);

    // return true if segments intersect each other
    if (o1 * o2 < 0 && o3 * o4 < 0)
        return true;

    // special cases
    if (o1 == 0 && vertex_on_segment(p1, q1, p2))
        return true;
    if (o2 == 0 && vertex_on_segment(p1, q2, p2))
        return true;
    if (o3 == 0 && vertex_on_segment(q1, p1, q2))
        return true;
    if (o4 == 0 && vertex_on_segment(q1, p2, q2))
        return true;

    // intersect test passed
    return false;
}

NUMDEF bool is_simple_polygon(vec3_t *vertices, size_t nvtx)
{
    // check if enough vertices to form a polygon
    if (nvtx < 3)
    {
        fprintf(stderr, "%-8sNumber of vertices = %zu < 3. At least 3 vertices are needed to form a polygon\n", "[ERROR]", nvtx);
        return false;
    }

    // chech coplanarity
    if (!is_coplanar(vertices, nvtx))
    {
        fprintf(stderr, "%-8sCoplanarity test failed\n", "[ERROR]");
        return false;
    }

    // project vertices to 2D
    vertices_proj_to_2d_in_place(vertices, nvtx);

    // check intersect of every pair of non-adjacent edges
    for (size_t i = 0; i < nvtx; ++i)
    {
        vec3_t p1 = vertices[i];
        vec3_t p2 = vertices[(i + 1) % nvtx];

        for (size_t j = i + 1; j < nvtx; ++j)
        {
            // skip edges shearing a vertex
            if (j == i || j == (i + 1) % nvtx || (i == 0 && j == nvtx - 1))
                continue;

            vec3_t q1 = vertices[j];
            vec3_t q2 = vertices[(j + 1) % nvtx];

            if (is_intersect(&p1, &p2, &q1, &q2))
            {
                fprintf(stderr, "%-8sIntersection found between segments (%zu-%zu) and (%zu-%zu)\n", "[ERROR]", i, (i + 1) % nvtx, j, (j + 1) % nvtx);
                return false;
            }
        }
    }

    // polygon simplicity test passed
    return true;
}

#endif // NUMERICAL_IMPLEMENTATION
