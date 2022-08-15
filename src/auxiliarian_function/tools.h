#ifndef _TOOLS_ORIENTATION_
#define _TOOLS_ORIENTATION_

#include <vector>
#include <algorithm>
//#include "../../Indirect_Predicates/implicit_point.h"
#include <iostream>
#include <cmath>
#include <tuple>

#include "debug.h"
#include <new>

// template<typename T>
// void fusion_sort(std::vector<T>& a, bool(*compare)(T, T));

typedef std::vector<double> points;

std::vector<std::vector<double>> matrix_product(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);
points matrix_product(const std::vector<std::vector<double>> A,const points X);
points intersection_line_plane(points a, points b, points x, points y, points z);
double euclidean_distance(std::vector<double>& u, std::vector<double>& v);
std::vector<double> barycenter_polygon2D(std::vector<std::vector<double>> vertices);
std::vector<double> barycenter_polygon3D(std::vector<std::vector<double>> vertices);
std::vector<double> centroid_3D(std::vector<std::vector<double>>);
std::vector<double> barycenter_points(std::vector<std::vector<double>> vertices, std::vector<double> weight);
std::vector<std::vector<double>> generate_basis(std::vector<std::vector<double>> vertices);
std::vector<double> min_vect(std::vector<double> a, std::vector<double> b);
std::vector<double> max_vect(std::vector<double> a, std::vector<double> b);
std::vector<double> min_vect(std::vector<std::vector<double>> a);
std::vector<int> min_vect_ind(std::vector<std::vector<double>> a);
double area_triangle(points a, points b, points c);
double volume_tet(points a, points b, points c, points d);

std::vector<std::vector<double>> rotation_matrix(double teta1, double teta2, double teta3);
std::vector<std::vector<double>> inverse_rotation_matrix(double teta1, double teta2, double teta3);


double random_float();


struct Face_add{
    double area = -1;
    std::vector<double> barycenter = {};
    std::vector<std::vector<double>> vertex_coord = {};
    std::vector<uint32_t> vertex_ind = {};
    std::vector<std::vector<double>> basis = {};
    std::vector<double> min = {};
    std::vector<double> max = {};
};



#endif