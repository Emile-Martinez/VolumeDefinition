#include "tools.h"


double scalar_product(double* u, double* v, int size){
    double res = 0;
    for (int i = 0; i < size; i++)
    {
        res += u[i]*v[i];
    }
    return res;
}

double scalar_product(points u, points v){
    double res = 0;
    int size = std::min(u.size(), v.size());
    for (int i = 0; i < size; i++)
    {
        res += u[i]*v[i];
    }
    return res;
}

// Do the sum of two vectors, truncated at the minimum of their sizes
points sum(points a, points b){
    int size = std::min(a.size(), b.size());
    points res(size, 0);
    for (int  i = 0; i < size; i++)
    {
        res[i] = a[i] + b[i];
    }
    return res;
}

points minus(points& a, points& b){
    int size = std::min(a.size(), b.size());
    points res(size, 0);
    for (int  i = 0; i < size; i++)
    {
        res[i] = a[i] - b[i];
    }
    return res;
}

// Multiply a vector by a scalar
points scalar_mult(double lambda, points v){
    for (int i = 0; i < v.size(); i++)
    {
        v[i] = v[i] * lambda;
    }
    return v;
}

// Return the euclidean distance between the two vectors, truncated at the minimum of their two sizes
// Note that if one on the vertex is empty, the output is 0
double euclidean_distance(points& u, points& v){
    int size = std::min(u.size(), v.size());
    double res = 0;
    for (int i = 0; i < size; i++)
    {
        double x = u[i] - v[i];
        res += x*x;
    }
    return std::sqrt(res);
}

std::vector<std::vector<double>> matrix_product(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B){
    std::vector<std::vector<double>> res(A.size(), std::vector<double>(B[0].size(), 0));

    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < B[0].size(); j++)
        {
            for (int k = 0; k < B.size(); k++)
            {
                res[i][j] += A[i][k] * B[k][j];
            }   
        }
    }
    return res;
}

points matrix_product(const std::vector<std::vector<double>> A,const points X){
    points res(A.size(), 0);
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < X.size(); j++)
        {
            res[i] += A[i][j] * X[j];
        }   
    }
    return res;
}

// Output the norm of a vector. Note that empty vector has 0 as norm
double norm(points v){
    double res = 0;
    for (int i = 0; i < v.size(); i++)
    {
        res += v[i]*v[i];
    }
    return std::sqrt(res);
}

double sign(double a){
    if (a== 0) return 0;
    else if (a<0) return -1;
    else if (a>0) return 1;
    else return a; // if it is nan for example
}

points cross_product(points a, points b){
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}


double area_triangle(points a, points b, points c){
    return norm(cross_product(minus(a, b), minus(a, c)))/2;
}


double volume_tet(points a, points b, points c, points d){
    return std::fabs(scalar_product(minus(a,d), cross_product(minus(b, d), minus(c, d))))/6;
}

// Inverse only 3*3 matrix
std::vector<std::vector<double>> inversion_matrix(std::vector<std::vector<double>> m){
    std::vector<std::vector<double>> res(3, std::vector<double>(3,0));

    double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
             m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
             m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if (det == 0) return res;
    double invdet = 1 / det;

    res[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    res[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
    res[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
    res[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    res[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
    res[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
    res[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
    res[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
    res[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

    return res;
}


// We calculate the barycenter of a 2D polygon using the wikipedia formula
// https://en.wikipedia.org/wiki/Centroid#Of_a_polygon

points barycenter_polygon2D(std::vector<points> vertices2D){
    int N = vertices2D.size();

    double A = 0, Cx = 0, Cy = 0;
    
    for (int i = 0; i < N; i++)
    {
        double x_i = vertices2D[i][0], y_i = vertices2D[i][1], x_j = vertices2D[(i+1)%N][0], y_j = vertices2D[(i+1)%N][1];
        double area = x_i*y_j - x_j*y_i;
        A += area;
        Cx += (x_i+x_j)*area;
        Cy += (y_i+y_j)*area;
    }
    
    Cx /= 3*A;
    Cy /= 3*A;

    return {Cx, Cy};

}


// Compute the centroid of a 3D convex polygon
points barycenter_polygon3D(std::vector<points> vertices){
    int N = vertices.size();

    double total_area = 0;
    points res = {0,0,0};

    points v0 = vertices[0];

    for (int i = 1; i < N-1; i++)
    {
        points v1 = vertices[i], v2 = vertices[i+1];
        auto u = minus(v1, v0), v = minus(v2, v1);
        points area_vertex = cross_product(u,v);
        double area = norm(area_vertex);
        points local_centroid = sum(v0, v1);
        local_centroid = sum(local_centroid, v2);
        local_centroid = scalar_mult(area/3., local_centroid);
        total_area += area;
        res = sum(local_centroid, res);
    }
    
    return scalar_mult(1./total_area, res);

}

//Not the true centroid but still a good enough approximation
points centroid_3D(std::vector<points> vertices){
    points res = {0,0,0};
    for (int i = 0; i < vertices.size(); i++)
    {
        res = sum(res, vertices[i]);
    }
    return scalar_mult(1./((double) vertices.size()), res);
}


points centroid_tetrahedron(points a, points b, points c, points d, double &total_area){
    points u = minus(a, b), v = minus(a, c), w = minus(a, d);
    points intermediary = cross_product(u, v);
    double area = std::fabs(1./6. * scalar_product(intermediary, w));
    total_area += area;
    intermediary = sum(a, b);
    points res = sum(intermediary, c);
    return scalar_mult(area, res);
}


points centroid_pyramid(std::vector<points> face, points apex, double &total_area){
    points res={0,0,0};
    for (int i = 1; i < face.size()-1; i++)
    {
        res = sum(res, centroid_tetrahedron(face[0], face[i], face[i+1], apex, total_area));
    }
    return res;
}

std::vector<double> barycenter_points(std::vector<points> vertices, std::vector<double> weight){

    points res = {0,0,0};
    double total_coeff = 0;
    for (int i = 0; i < vertices.size(); i++)
    {
        res = sum(res, scalar_mult(weight[i], vertices[i]));
        total_coeff += weight[i];
    }
    return scalar_mult(1./total_coeff, res);
}

// Return a orthonormal basis of the plan generated by the input points, which are assumed colinear
// It returns also a a normal vector
std::vector<points> generate_basis(std::vector<points> vertices){
    
    int N = vertices.size();

    // We will work on 2D, to apply the right formula, and thus project eveything on a 2D plan
    std::vector<double> v1(N), v2(N); //It will be the two vectors of the 2D base

    v1 = minus(vertices[0], vertices[1]);
    v1 = scalar_mult(1./norm(v1), v1);

    // Now we have the first one. We now just need an other non coplanar and turns all this in an orthonormal basis
    // In order to reduce rounding errors, we will take the "less parallel" of the vectors made of the first vertex
    // of the polygon and anyother one

    std::vector<double> normalized_scalar_product(N-1);
    for (int i = 1; i < N; i++)
    {
        auto intermediate = minus(vertices[0], vertices[i]);
        normalized_scalar_product[i-1] = std::fabs(scalar_product(v1, intermediate)) / norm(intermediate);
    }

    int minimal_direction_index = std::min_element(normalized_scalar_product.begin(), normalized_scalar_product.end()) - normalized_scalar_product.begin() +1;

    v2 = minus(vertices[0], vertices[minimal_direction_index]);
    v2 = scalar_mult(1./norm(v2), v2);
    auto intermediate2 = scalar_mult(scalar_product(v1,v2), v2);
    v2 = minus(v2, intermediate2);
    v2 = scalar_mult(1./norm(v2), v2);
    
    return {v1, v2, cross_product(v1, v2)};
}



std::vector<double> min_vect(std::vector<double> a, std::vector<double> b){
    int N = std::min(a.size(), b.size());
    std::vector<double> res(a.size());
    for (int i = 0; i < N; i++)
    {
        res[i] = std::min(a[i], b[i]);
    }
    return res;
}

std::vector<double> min_vect(std::vector<std::vector<double>> a){
    auto res = a[0];
    for (auto &&x : a)
    {
        res = min_vect(x, res);
    }
    return res;
}

std::vector<int> min_vect_ind(std::vector<std::vector<double>> a){
    auto min = a[0];
    std::vector<int> res(a[0].size(), 0);
    for (int i = 1; i < a.size(); i++)
    {
        for (int j = 0; j < a[0].size(); j++)
        {
            if (min[j] > a[i][j]) {
                min[j] = a[i][j];
                res[j] = i;
            }
        }   
    }
    return res;
}

std::vector<double> max_vect(std::vector<double> a, std::vector<double> b){
    int N = std::min(a.size(), b.size());
    std::vector<double> res(a.size());
    for (int i = 0; i < N; i++)
    {
        res[i] = std::max(a[i], b[i]);
    }
    return res;
}


std::vector<std::vector<double>> rotation_matrix(double teta1, double teta2, double teta3){

    std::vector<std::vector<double>> rot1, rot2, rot3;
    rot1 = {{1 ,0 ,0}, {0, std::cos(teta1), -std::sin(teta1)}, {0, std::sin(teta1), std::cos(teta1)}};
    rot2 = {{std::cos(teta2), 0, std::sin(teta2)}, {0, 1, 0}, {-std::sin(teta2), 0, std::cos(teta2)}};
    rot3 = {{std::cos(teta3), -std::sin(teta3), 0}, {std::sin(teta3), std::cos(teta3), 0}, {0, 0, 1}};
    return matrix_product(rot1, matrix_product(rot2, rot3));
}


std::vector<std::vector<double>> inverse_rotation_matrix(double teta1, double teta2, double teta3){
    std::vector<std::vector<double>> rot1, rot2, rot3;
    rot1 = {{1 ,0 ,0}, {0, std::cos(-teta1), -std::sin(-teta1)}, {0, std::sin(-teta1), std::cos(-teta1)}};
    rot2 = {{std::cos(-teta2), 0, std::sin(-teta2)}, {0, 1, 0}, {-std::sin(-teta2), 0, std::cos(-teta2)}};
    rot3 = {{std::cos(-teta3), -std::sin(-teta3), 0}, {std::sin(-teta3), std::cos(-teta3), 0}, {0, 0, 1}};
    return matrix_product(rot3, matrix_product(rot2, rot1));
}



double random_float(){
    return ((float) std::rand())/ ((float) RAND_MAX);
}

// Gives the intersection between the line (a,b) and the plane (x,y,z)
points intersection_line_plane(points a, points b, points x, points y, points z){
    std::vector<std::vector<double>> m = {{a[0]-b[0], a[1]-b[1], a[2]-b[2]},
                                            {x[0]-y[0], x[1]-y[1], x[2]-y[2]},
                                            {x[0]-z[0], x[1]-z[1], x[2]-z[2]}
                                           };
    

    double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if (det == 0) return a;
    double invdet = 1 / det;

    
    double inv00 = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    double inv10 = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    double inv20 = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;

    double t = inv00*(a[0]-x[0]) + inv10*(a[1]-x[1]) + inv20*(a[2]-x[2]);

    return {a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1]), a[2] + t*(b[2]-a[2])};
    
}
