#include "collision.h"


bool intersection_segment_box(points a, points b, points min_bound, points max_bound){
    double t_min = 0, t_max = 1;
    for (int i = 0; i < 3; i++)
    {
        if (b[i] - a[i] > 0){
            t_min = std::max(t_min, (min_bound[i]-a[i])/(b[i]-a[i]));
            t_max = std::min(t_max, (max_bound[i]-a[i])/(b[i]-a[i]));
        }
        else{
            t_min = std::max(t_min, (max_bound[i]-a[i])/(b[i]-a[i]));
            t_max = std::min(t_max, (min_bound[i]-a[i])/(b[i]-a[i]));
        }
    }

    return t_min < t_max;
    
}

bool intersection_box_tri(Triangle T, const points min_bound, const points max_bound){
    float boxcenter[3], boxhalfsize[3];
    boxcenter[0] = (min_bound[0] + max_bound[0])/2.;
    boxcenter[1] = (min_bound[1] + max_bound[1])/2.;
    boxcenter[2] = (min_bound[2] + max_bound[2])/2.;
    boxhalfsize[0] = (max_bound[0] - min_bound[0])/2.;
    boxhalfsize[1] = (max_bound[1] - min_bound[1])/2.;
    boxhalfsize[2] = (max_bound[2] - min_bound[2])/2.;

    return triBoxOverlap(boxcenter, boxhalfsize, T.vertices_f);
}

Octree::Octree(std::vector<Triangle*>& input, points min_input, points max_input):min_bound(min_input), max_bound(max_input), in_triangle(input){
    start = false;
    if (input.size() < 5){
        end = true;
    }
    else{
    end = false;
    
    for (int i = 0; i < 8; i++)
    {
        std::vector<Triangle*> son_triangles(0);
        points new_min(3), new_max(3);
        if((i+1)%4 < 2){
            new_min[0] = min_input[0];
            new_max[0] = (max_input[0] + min_input[0])/2;
        }
        else{
            new_min[0] = (max_input[0] + min_input[0])/2;
            new_max[0] = max_input[0];
        }

        if(i < 4){
            new_min[1] = min_input[1];
            new_max[1] = (max_input[1] + min_input[1])/2;
        }
        else{
            new_min[1] = (max_input[1] + min_input[1])/2;
            new_max[1] = max_input[1];
        }

        if(i%4 < 2){
            new_min[2] = min_input[2];
            new_max[2] = (max_input[2] + min_input[2])/2;
        }
        else{
            new_min[2] = (max_input[2] + min_input[2])/2;
            new_max[2] = max_input[2];
        }

        for(Triangle* tri: input){
            if(intersection_box_tri(*tri, new_min, new_max)){
                son_triangles.push_back(tri);
            }
        }

        if(son_triangles.size() == in_triangle.size()) {
            end = true; 
            for (int j = 0; j < i; j++)
            {
                delete son[j];
            }
            break;
        }

        son[i] = new Octree(son_triangles, new_min, new_max);
        
    }

    }
}

Octree::Octree(std::vector<std::vector<points>> triangles){
    time_tri = 0; time_box = 0;
    start = true;
    end = false;
    points min = triangles[0][0];
    points max = min;
    in_triangle = {};

    for (int i = 0; i < triangles.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            min[0] = std::min(min[0], triangles[i][j][0]);
            min[1] = std::min(min[1], triangles[i][j][1]);
            min[2] = std::min(min[2], triangles[i][j][2]);
            max[0] = std::max(max[0], triangles[i][j][0]);
            max[1] = std::max(max[1], triangles[i][j][1]);
            max[2] = std::max(max[2], triangles[i][j][2]);
        }
        in_triangle.push_back(new Triangle(triangles[i][0], triangles[i][1], triangles[i][2]));
    }

    min_bound = min;
    max_bound = max;

    for (int i = 0; i < 8; i++)
    {
        std::vector<Triangle*> son_triangles(0);
        points new_min(3), new_max(3);
        if((i+1)%4 < 2){
            new_min[0] = min_bound[0];
            new_max[0] = (max_bound[0] + min_bound[0])/2;
        }
        else{
            new_min[0] = (max_bound[0] + min_bound[0])/2;
            new_max[0] = max_bound[0];
        }

        if(i < 4){
            new_min[1] = min_bound[1];
            new_max[1] = (max_bound[1] + min_bound[1])/2;
        }
        else{
            new_min[1] = (max_bound[1] + min_bound[1])/2;
            new_max[1] = max_bound[1];
        }

        if(i%4 < 2){
            new_min[2] = min_bound[2];
            new_max[2] = (max_bound[2] + min_bound[2])/2;
        }
        else{
            new_min[2] = (max_bound[2] + min_bound[2])/2;
            new_max[2] = max_bound[2];
        }

        for(Triangle* tri: in_triangle){
            if(intersection_box_tri(*tri, new_min, new_max)){
                son_triangles.push_back(tri);
            }
        }
        son[i] = new Octree(son_triangles, new_min, new_max);
    }


}

Octree::~Octree(){
    if(!end){
        for (int i = 0; i < 8; i++)
        {
            delete son[i];
        }
    }
    if(start){
        for (auto t: in_triangle)
        {
            delete t;
        }
        std::cout << "Time for triangle intersection : " << time_tri/1000 << "us\n";
        std::cout << "Time for box intersection : " << time_box/1000 << "us\n";
    }
}

inline bool segment_IntersectTriangle(points x, points y, Triangle *T){
    points a = T->a, b = T->b, c = T->c;
    //std::cerr << orient3d(x.data(), y.data(), a.data(), b.data()) << orient3d(x.data(), y.data(), b.data(), c.data()) << orient3d(x.data(), y.data(), c.data(), a.data()) << orient3d(x.data(), a.data(), b.data(), c.data()) << orient3d(y.data(), a.data(), b.data(), c.data());
    return (orient3d(x.data(), y.data(), a.data(), b.data()) == orient3d(x.data(), y.data(), b.data(), c.data()))
        && (orient3d(x.data(), y.data(), b.data(), c.data()) == orient3d(x.data(), y.data(), c.data(), a.data()))
        && (orient3d(x.data(), a.data(), b.data(), c.data()) != orient3d(y.data(), a.data(), b.data(), c.data()));
}


int Octree::number_intersection(points a, points b){
    if (end){
        for(Triangle* T:in_triangle){
            auto start = std::chrono::high_resolution_clock::now();
            T->touched = segment_IntersectTriangle(a, b, T);
            auto end = std::chrono::high_resolution_clock::now();
            time_tri += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        }
    }
    else{
        for (int i = 0; i < 8; i++)
        {
            auto start = std::chrono::high_resolution_clock::now();
            bool keep_going = intersection_segment_box(a, b, son[i]->min_bound, son[i]->max_bound);
            auto end = std::chrono::high_resolution_clock::now();
            time_box += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
            if(keep_going)
                    son[i]->number_intersection(a, b);
            }
    }
    if (start){
        int res = 0;
        for(Triangle* T : in_triangle){
            res += T->touched;
            T->touched = false;
        }
        return res;
    }
    return 0;
}


long int Octree::time_box = 0;
long int Octree::time_tri = 0;