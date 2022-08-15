#ifndef COLLISION_H
#define COLLISION_H

#include <vector>
#include <iostream>
#include <chrono>
#include "../delaunay.h"

int triBoxOverlap(float boxcenter[3],float boxhalfsize[3],float triverts[3][3]);

typedef std::vector<double> points;

// inline double orient3d(const double* p, const double* q, const double* r, const double* s);
// //int orient3d(const double* p, const double* q, const double* r, const double* s);

bool intersection_segment_box(points a, points b, points min_bound, points max_bound);

class Triangle{
    public:
        points a, b, c;
        float vertices_f[3][3];
        points barycenter;
        bool touched = false;
        Triangle(points x, points y, points z): a(x), b(y), c(z){
            vertices_f[0][0] = x[0]; vertices_f[0][1] = x[1], vertices_f[0][2] = x[2];
            vertices_f[1][0] = y[0]; vertices_f[1][1] = y[1], vertices_f[1][2] = y[2];
            vertices_f[2][0] = z[0]; vertices_f[2][1] = z[1], vertices_f[2][2] = z[2];
        }
};

bool intersection_box_tri(Triangle T, const points min_bound, const points max_bound);

class Octree{
    bool end;
    bool start;
    points min_bound, max_bound;
    std::vector<Triangle*> in_triangle;
    Octree* son[8];
    Octree(std::vector<Triangle*>& intput, points min_input, points max_input);
        
    public:
        static long int time_box;
        static long int time_tri;
        Octree(std::vector<std::vector<points>> triangles);
        int number_intersection(points a, points b);
        ~Octree();
};



#endif