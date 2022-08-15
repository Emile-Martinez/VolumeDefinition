#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <ctime>
#include <chrono>
#include <cmath>
#include <map>
#include "BSP.h"
#include <math.h>
#include <cassert>

#include "graph_cut/GCoptimization.h"
#include "graph_cut/LinkedBlockList.cpp"
#include "graph_cut/GCoptimization.cpp"

//#define Time_info

#define OUTSIDE UINT32_MAX
#define IS_GHOST_CELL(c) (c==UINT64_MAX)
#define GHOST_CELL UINT64_MAX
#define IS_NOT_OPPOSITE(a, b) (a!=-b)
#define OTHER_ONE(arr, value) (arr[0] == value) ? arr[1] : arr[0]
#define IS_IN(arr, value) ((arr[0] == value) || (arr[1] == value))
#define ROTATE_ARG1(a, b, c) c, a, b
#define ROTATE_ARG2(a, b, c) b, c, a


#define DETERMINANT3X3(a11, a12, a13, a21, a22, a23, a31, a32, a33) ((a11)*((a22)*(a33) - (a23)*(a32)) - (a12)*((a21)*(a33) - (a23)*(a31)) + (a13)*((a21)*(a32) - (a22)*(a31)))

// Works only for convex faces, which is supposed to be always the case
double approxFaceArea(BSPface& face, BSPcomplex* cpx, const std::vector<double>& approxCoords)
{
    std::vector<uint32_t> vs(face.edges.size(), 0);
    cpx->list_faceVertices(face, vs);

    const double* acp = approxCoords.data();

    const double* tv0, * tv1, * tv2;
    double a = 0.0;
    tv0 = acp + vs[0] * 3;
    for (size_t i = 2; i < vs.size(); i++)
    {
        tv1 = acp + vs[i - 1] * 3;
        tv2 = acp + vs[i    ] * 3;
        double x, y, z;
        x = (tv0[1]-tv1[1])*(tv0[2]-tv2[2]) - (tv0[2]-tv1[2])*(tv0[1]-tv2[1]);
        y = (tv0[2]-tv1[2])*(tv0[0]-tv2[0]) - (tv0[0]-tv1[0])*(tv0[2]-tv2[2]);
        z = (tv0[0]-tv1[0])*(tv0[1]-tv2[1]) - (tv0[1]-tv1[1])*(tv0[0]-tv2[0]);
        a += std::sqrt(x*x+y*y+z*z);
    }

    return fabs(a)/2;
}

double approxFaceArea(BSPface& face, BSPcomplex* cpx, const std::vector<points>& approxCoords)
{
    std::vector<uint32_t> vs(face.edges.size(), 0);
    cpx->list_faceVertices(face, vs);


    points tv0, tv1, tv2;
    double a = 0.0;
    tv0 = approxCoords[vs[0]];
    for (size_t i = 2; i < vs.size(); i++)
    {
        tv1 = approxCoords[vs[i-1]];
        tv2 = approxCoords[vs[i]];
        a += area_triangle(tv0, tv1, tv2);
    }

    return fabs(a);
}



// Returns TRUE if face is part of the skin according to skin_colour
inline bool isSkinFace(const BSPface& face, uint32_t skin_colour)
{
    return face.colour & skin_colour;
}

inline void setInternalCell(BSPcell& c, uint32_t internal_label)
{
    c.place |= internal_label;
}

inline void setExternalCell(BSPcell& c)
{
    c.place = EXTERNAL;
}

// ########################################### DEBUG ###################################################
#pragma region 
void save_faces(BSPcomplex cpx, std::vector<int> faces){
    for (BSPface &face : cpx.faces) face.colour = WHITE;
    for (int face : faces)
    {
        cpx.faces[face].colour = BLACK_A;
    }
    cpx.saveBlackFaces("specified_faces.off");
    exit(0);
}

void save_faces(BSPcomplex cpx){
    for (BSPface &face : cpx.faces) face.colour = BLACK_A;
    cpx.saveBlackFaces("specified_faces.off");
    exit(0);
}

void save_tetrahedron(BSPcomplex cpx, std::vector<int> cells){
    for (BSPcell& cell : cpx.cells) cell.place = EXTERNAL;
    for (int i : cells)
    {
        cpx.cells[i].place = INTERNAL_A;
    }
    cpx.saveMesh("Speficied_mesh.msh", '0');
    exit(0);
}

void save_tetrahedron(BSPcomplex cpx){
    for (BSPcell& cell : cpx.cells) cell.place = INTERNAL_A;
    cpx.saveMesh("Specified_mesh.msh", '0');
    exit(0);
}
#pragma endregion
// ######################################## END OF DEBUG ###############################################

#pragma region 

// Label = 1 means INTERNAL
// Label = 0 means EXTERNAL

// Take first coplanar constraint associated to this face
// and return TRUE if the cell vertices are 'below' such constraint.
bool isFirstConnCellBelowFace(BSPface& f, BSPcomplex* cpx)
{
    const uint32_t* cid = cpx->constraints_vrts.data() + f.coplanar_constraints[0] * 3;
    const genericPoint* pv1 = cpx->vertices[cid[0]];
    const genericPoint* pv2 = cpx->vertices[cid[1]];
    const genericPoint* pv3 = cpx->vertices[cid[2]];

    BSPcell& cell = cpx->cells[f.conn_cells[0]];
    uint64_t num_cellEdges = UINT64_MAX;
    uint32_t num_cellVrts = cpx->count_cellVertices(cell, &num_cellEdges);
    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
    cpx->list_cellVertices(cell, num_cellEdges, cell_vrts);

    for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 1;

    for (uint32_t vi : cell_vrts) if (!cpx->vrts_visit[vi])
    {
        const genericPoint* cv = cpx->vertices[vi];
        const int o = genericPoint::orient3D(*cv, *pv1, *pv2, *pv3);
        if (o)
        {
            for (uint64_t ei : f.edges) cpx->vrts_visit[cpx->edges[ei].vertices[0]] = cpx->vrts_visit[cpx->edges[ei].vertices[1]] = 0;
            return (o > 0);
        }
    }

    ip_error("Degenerate cell\n");
    return false;
}


void BSPcomplex::markInternalCells(uint32_t skin_colour, uint32_t internal_label, const std::vector<double>& face_areas)
{

    // auto graph = get_graph(this, skin_colour);
    // auto distance_to_outside = graph.djikstra(OUTSIDE);


    // Allocate dual graph: num cells + 1 to account for the external "ghost" cell
    GCoptimizationGeneralGraph gc((GCoptimization::SiteID)cells.size() + 1, 2);

    // gc is the dual graph of the cell complex
    // - a node in gc corresponds to a cell in the complex
    // - an arc in gc exists if two cells share a WHITE face

    // In gc, cells that share a white face are connected by an arc weighted on face area
    // The 'data cost' associated to each cell is:
    //  - total area of BLACK faces 'consistently oriented' with the cell, if label is EXTERNAL
    //  - 0, if label is INTERNAL
    //
    // The 'smooth cost' associated to each WHITE face is:
    //  - the face area, if label_A and label_B are different
    //  - 0 otherwise
    //
    // Note: a BLACK face is considered to be consistently oriented with one of its incident cells
    // if the cell is 'below' the first of its coplanar constraints.

    // evs == 1 if edge is on boundary of skin
    std::vector<uint8_t> evs(edges.size(), 0);
    for (BSPface& f : faces) if (isSkinFace(f, skin_colour))
      for (uint64_t eid : f.edges) if (evs[eid] < 2) evs[eid]++;

    // vvs == 1 if vertex is on boundary of skin
    std::vector<uint8_t> vvs(vertices.size(), 0);
    for (size_t i = 0; i < edges.size(); i++) if (evs[i] == 1)
    {
        const BSPedge& e = edges[i];
        vvs[e.vertices[0]] = vvs[e.vertices[1]] = 1;
    }

    std::vector<double> cell_costs_external(cells.size() + 1, 0.0);
    std::vector<double> cell_costs_internal(cells.size() + 1, 0.0);

    int F = faces.size();

    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);

        if (isSkinFace(f, skin_colour))
        {
            if (isFirstConnCellBelowFace(f, this))
            {
                cell_costs_external[cell1] += face_areas[i];
                cell_costs_internal[cell2] += face_areas[i];
            }
            else
            {
                cell_costs_external[cell2] += face_areas[i];
                cell_costs_internal[cell1] += face_areas[i];
            }
        }
    }

    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        if (!isSkinFace(f, skin_colour))
        {
            // 'w' is an additional weight for arcs that promotes the cut of arcs corresp. to faces having
            // all their vertices on the boundary of the input surface (i.e. promote hole filling)
            double w = 0.1;
            for (uint64_t eid : f.edges) if (vvs[edges[eid].vertices[0]] == 0 || vvs[edges[eid].vertices[1]] == 0) { w = 1.0; break; }
            const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);
            gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, face_areas[i]*w);
        }
    }

    const double int_weight = 0.1; // Internal cell penalization less than external to avoid artifacts at intersections
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 0, cell_costs_external[i]);
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 1, cell_costs_internal[i] * int_weight);
    gc.setDataCost((GCoptimization::SiteID)cells.size(), 1, 1.0); // Ghost cell must be external

    // Run graph cut algorithm
    // I.e., label all the cells so that the total data cost + smooth cost is minimized
    gc.swap();

    for (size_t i = 0; i < cells.size(); i++)
        if (gc.whatLabel((GCoptimization::SiteID)i)) setInternalCell(cells[i], internal_label);
}

void BSPcomplex::constraintsSurface_complexPartition(bool two_files)
{


    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);

    // Clear vrts_visit for use in isFirstConnCellBelowFace()
    for (size_t i = 0; i < vertices.size(); i++) vrts_visit[i] = 0;

    // Precalculate approximate vertex coordinates for use in approxFaceArea()
    std::vector<double> approxCoords(vertices.size() * 3);
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i * 3], approxCoords[i * 3 + 1], approxCoords[i * 3 + 2]);

    // Precalculate approximate face areas for use in markInternalCells()
    std::vector<double> face_areas(faces.size(), 0.0);
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_areas[i] = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_areas[i];
    }


    // Normalize all areas to avoid overflows in graphcut
    for (size_t i = 0; i < faces.size(); i++) face_areas[i] /= tot_face_area;

    if (two_files)
    {
        markInternalCells(BLACK_A, INTERNAL_A, face_areas);
        markInternalCells(BLACK_B, INTERNAL_B, face_areas);
    }
    else markInternalCells(BLACK_A, INTERNAL_A, face_areas);

}

#pragma endregion


// Smoothings
#pragma region 

void BSPcomplex::smoothing1(std::vector<int> inside, std::vector<int> outside, std::vector<Face_add>& face_info, uint32_t internal_label, uint32_t skin_colour){
    // We are going to turn the complex in its dual graph

    std::vector<int> is_inside(cells.size()+1, 0);

    // is_inside[i] is true if and only if we had more inside ray than outside ray for the cell i (i.e. we want to classify it inside)
    for (size_t i = 0; i < cells.size(); i++) is_inside[i] = (inside[i] > outside[i]);

    // We then initialize the graph
    Minimize_area graph(is_inside);

    // cells.size() is the index of the outside cell (i.e. R^3 minus the convex hull of the input)
    // And this cell must of course remain outside
    graph.set_unchangeable(cells.size());

    for (int ind_f = 0; ind_f < faces.size(); ind_f++)
    {
        BSPface f = faces[ind_f];
        int u = f.conn_cells[0];
        int v = f.conn_cells[1];
        if(IS_GHOST_CELL(u)) u = cells.size();
        if(IS_GHOST_CELL(v)) v = cells.size();

        // We connect two points if they represent two cells which are connected
        // We take the area has weight, with a negative weight if the facets between the two cells is black
        if (isSkinFace(f, skin_colour)){
            graph.add_edge(u, v, -face_info[ind_f].area);
        }
        else{
            graph.add_edge(u, v, face_info[ind_f].area);
        }
    }
    
    // Then we perform the smoothing
    graph.smooth();

    // We retreive the classification
    for (size_t i = 0; i < cells.size(); i++) if(graph.classification[i]) setInternalCell(cells[i], internal_label);
}

void BSPcomplex::smoothing2(std::vector<int> inside, std::vector<int> outside, std::vector<Face_add>& face_info, uint32_t internal_label, uint32_t skin_colour){
    // We turn the complex in its dual graph

    Smoothing_graph graph(inside, outside);
    
    for (int ind_f = 0; ind_f < faces.size(); ind_f++)
    {
        BSPface f = faces[ind_f];

        if(!(IS_GHOST_CELL(f.conn_cells[0]) || IS_GHOST_CELL(f.conn_cells[1]))){
        if (isSkinFace(f, skin_colour)){
            graph.add_edge(f.conn_cells[0], f.conn_cells[1], -face_info[ind_f].area);
        }
        else{
            graph.add_edge(f.conn_cells[0], f.conn_cells[1], face_info[ind_f].area);
        }
        }
    }
    

    // ofstream f("number_of_ray_per_cell.hist");

    // for (int i = 0; i < cells.size(); i++)
    // {
    //     f << graph.inside[i] + graph.outside[i] << " ";
    // }
    // f.close();

    // We get back the potential additional paramater for the smoothing
    int number_of_iteration = 1000;
    if(global_parameters.size()>2) number_of_iteration = global_parameters[2];
    double proportion = 0.0;
    if(global_parameters.size()>3) proportion = global_parameters[3]/100.;

    // We do a certain number of iteration depending on if we want to reach convergence or no
    for (int i = 0; i < number_of_iteration; i++)
    {
        graph.iteration(proportion);
    }

    // Then with smoothed values of inside and outside, we get back the classification
    for (int i = 0; i < cells.size(); i++)
    {
        if (graph.inside[i] > graph.outside[i]){
            setInternalCell(cells[i], internal_label);
        }
    }

    // ofstream f2("number_of_ray_per_cell2.hist");

    // for (int i = 0; i < cells.size(); i++)
    // {
    //     f2 << graph.inside[i] + graph.outside[i] << " ";
    // }
    // f2.close();

}

void BSPcomplex::smoothing3(std::vector<int> inside, std::vector<int> outside, double area_sample, std::vector<Face_add>& face_info, uint32_t internal_label, uint32_t skin_colour){
    
    // We initialize the graph with two possible labels
    GCoptimizationGeneralGraph gc((GCoptimization::SiteID)cells.size() + 1, 2);

    // We add edges only between cells connected with a white facets
    for (size_t i = 0; i < faces.size(); i++)
    {
        BSPface& f = faces[i];
        const uint64_t cell1 = f.conn_cells[0];
        if (!isSkinFace(f, skin_colour))
        {
            const uint64_t cell2 = (f.conn_cells[1] == UINT64_MAX) ? (cells.size()) : (f.conn_cells[1]);
            gc.setNeighbors((GCoptimization::SiteID)cell1, (GCoptimization::SiteID)cell2, face_info[i].area);
        }
    }


    // for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, (int) (inside[i] < outside[i]), 2*(inside[i] < outside[i] ? 0.1*std::cbrt(outside[i]*outside[i]*outside[i] - inside[i]*inside[i]*inside[i]) : std::cbrt(inside[i]*inside[i]*inside[i] - outside[i]*outside[i]*outside[i]))*area_sample);
    // for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, (int) (inside[i] < outside[i]), 1.*(inside[i] < outside[i] ? std::sqrt(outside[i]*outside[i] - inside[i]*inside[i]) : std::sqrt(inside[i]*inside[i] - outside[i]*outside[i]))/3/number_of_rotation*area_sample);
    // for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, (int) (inside[i] < outside[i]), 0.1 *(inside[i] < outside[i] ? (outside[i] - inside[i]) : inside[i] - outside[i])/3/number_of_rotation*area_sample);

    // The label 0 mean Outside and 1 means Inside
    // We down weight the values of the ray to have a more proper output
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 0, 0.1*inside[i]*area_sample);
    for (size_t i = 0; i < cells.size(); i++) gc.setDataCost((GCoptimization::SiteID)i, 1, 0.1*outside[i]*area_sample);
    gc.setDataCost((GCoptimization::SiteID)cells.size(), 1, 10000000.0); // Ghost cell must be external

    // This mean that if the labels are the same, we multiply the weight in the mincut algorithm by 0, 
    // and by 1 if the label are different
    gc.setSmoothCost(0, 0, 0);
    gc.setSmoothCost(1, 1, 0);
    gc.setSmoothCost(0, 1, 1);
    gc.setSmoothCost(1, 0, 1);

    // Run graph cut algorithm
    // I.e., label all the cells so that the total data cost + smooth cost is minimized
    gc.swap();

    // We get back the classification
    for (size_t i = 0; i < cells.size(); i++)
        if (gc.whatLabel((GCoptimization::SiteID)i)) setInternalCell(cells[i], internal_label);
}


#pragma endregion

// Some geometric predicates
#pragma region 


bool segment_IntersectSegment(genericPoint* a, genericPoint* b, genericPoint* x, genericPoint*y){
    if (genericPoint::orient3D(*a, *b, *x, *y) == 0) return (genericPoint::segmentsCross(*a, *b, *x, *y));
    else return false;
}

#ifdef Time_info
long double intersection_time = 0;
#endif

// Return true if the closure [x,y] intersects the closorue of <a,b,c> (there can be issue if x, y, a, b, c are coplanar)
inline bool segment_IntersectTriangle(genericPoint *x, genericPoint *y, genericPoint *a, genericPoint *b, genericPoint *c){
    #ifdef Time_info
    auto start = std::chrono::high_resolution_clock::now();
    #endif
    bool res =genericPoint::innerSegmentCrossesInnerTriangle(*x, *y, *a, *b, *c) || segment_IntersectSegment(x, y, a, b) || segment_IntersectSegment (x, y, b, c) || segment_IntersectSegment(x, y, c, a);

    #ifdef Time_info
    auto end = std::chrono::high_resolution_clock::now();
    intersection_time += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    #endif
    return res;
}

// Return true if and only if the segment [a,b] intersect the face
bool BSPcomplex::segment_Intersectface(int face_ind, genericPoint* a, genericPoint* b){
    

    bool res = false;
    BSPface face = faces[face_ind];
    vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
    list_faceVertices(face, face_vrts);

    genericPoint *v0 = vertices[face_vrts[0]];
    for (int i = 0; i < face_vrts.size()-1; i++)
    {
        res = res || segment_IntersectTriangle(a, b, v0, vertices[face_vrts[i]], vertices[face_vrts[i+1]]);
    }

    return res;
}

inline bool segment_IntersectTriangle(points x, points y, points a, points b, points c){
    return (orient3d(x.data(), y.data(), a.data(), b.data()) == orient3d(x.data(), y.data(), b.data(), c.data()))
        && (orient3d(x.data(), y.data(), b.data(), c.data()) == orient3d(x.data(), y.data(), c.data(), a.data()))
        && (orient3d(x.data(), a.data(), b.data(), c.data()) != orient3d(y.data(), a.data(), b.data(), c.data()));
}
#pragma endregion

// First method, using AABB detector
#pragma region 

// Generate a point that can be considered at the infinite
points generate_point(){
    double x = random_float()-0.5, y = random_float()-0.5, z = random_float()-0.5;
    double norm = std::sqrt(x*x+y*y+z*z);
    return {1000*x/norm, 1000*y/norm, 1000*z/norm};
}

void BSPcomplex::markInternalCells_ray(uint32_t skin_colour, uint32_t internal_label, std::vector<points> barycenters, std::vector<std::vector<points>> triangles, std::vector<Face_add>& face_info, int smoothing){

    int N=20;
    if(global_parameters.size()>0) N = global_parameters[0];

    std::vector<points> outsides(N);

    for (int i = 0; i < N; i++)
    {
        outsides[i] = generate_point();
    }

    Octree detector(triangles);

    std::vector<int> threshold = {-1 , -1 , -1 , -1 , -1 , 0 , 0 , 0 , 1 , 1 , 1 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 5 , 5 , 5 , 6 , 6 , 7 , 7 , 7 , 8 , 8 , 9 , 9 , 10 , 10 , 10 , 11 , 11 , 12 , 12 , 13 , 13 , 13 , 14 , 14 , 15 , 15 , 16 , 16 , 16 , 17 , 17 , 18 , 18 , 19 , 19 , 20 , 20 , 20 , 21 , 21 , 22 , 22 , 23 , 23 , 24 , 24 , 24 , 25 , 25 , 26 , 26 , 27 , 27 , 28 , 28 , 28 , 29 , 29 , 30 , 30 , 31 , 31 , 32 , 32 , 33 , 33 , 33 , 34 , 34 , 35 , 35 , 36 , 36 , 37 , 37 , 38 , 38 , 38 , 39 , 39 , 40 , 40 , 41 , 41 , 42 , 42 , 43 , 43 , 44 , 44 , 44 , 45 , 45 , 46 , 46 , 47 , 47 , 48 , 48 , 49 , 49 , 50 , 50 , 50 , 51 , 51 , 52 , 52 , 53 , 53 , 54 , 54 , 55 , 55 , 56 , 56 , 56 , 57 , 57 , 58 , 58 , 59 , 59 , 60 , 60 , 61 , 61 , 62 , 62 , 63 , 63 , 63 , 64 , 64 , 65 , 65 , 66 , 66 , 67 , 67 , 68 , 68 , 69 , 69 , 70 , 70 , 70 , 71 , 71 , 72 , 72 , 73 , 73 , 74 , 74 , 75 , 75 , 76 , 76 , 77 , 77 , 78 , 78 , 78 , 79 , 79 , 80 , 80 , 81 , 81 , 82 , 82 , 83 , 83 , 84 , 84 , 85 , 85 , 85 , 86 , 86 , 87};

    int compteur_ray = 0;

    std::vector<int> inside(cells.size(), 0), outside(cells.size(), 0);

    // For each cells, we shoot ray until we go over the threshold or until we shoot to many rays
    for (int i = 0; i < cells.size(); i++){
        std::cerr << i << " ";
        int k=0;
        for (int j = 0; j < N; j++)
        {
            compteur_ray++;
            k += detector.number_intersection(barycenters[i], outsides[j])%2;
            inside[i] = k;
            outside[i] = j-k;
            if(k <= threshold[j]) {break;}
            if(k > j - threshold[j] ){break;}
        }
        
    }
    std::cout << "Number of shooted rays : " << compteur_ray / (double) cells.size() << "\n";

    if(smoothing == 1){
        smoothing1(inside, outside, face_info, internal_label, skin_colour);
    }
    else if(smoothing == 2){
        smoothing2(inside, outside, face_info, internal_label, skin_colour);
    }
    else{
        for (int i = 0; i < cells.size(); i++)
        {
            if(inside[i] > outside[i]) setInternalCell(cells[i], internal_label);
        }
    }

}

void BSPcomplex::constraintsSurface_complexPartition_ray_approx(std::vector<std::vector<points>> triangles, bool two_files, int smoothing){
    
    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);  

    std::vector<points> approxCoords(vertices.size(), points(3,0));
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i][0], approxCoords[i][1], approxCoords[i][2], true);


    // Precalculate approximate face informations for use in markInternalCells
    std::vector<Face_add> face_info(faces.size());
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_info[i].area = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_info[i].area;

        std::vector<uint32_t> vs(faces[i].edges.size(), 0);
        list_faceVertices(faces[i], vs);
        face_info[i].vertex_ind = vs;
    }

    // Let's precompute an approximation of the barycenter of each cells to know from where shoot ray
    std::vector<points> barycenters(cells.size());

    for (int i = 0; i < cells.size(); i++)
    {
        BSPcell cell = cells[i];
        uint64_t num_cellEdges = UINT64_MAX;
        uint32_t num_cellVrts = count_cellVertices(cell, &num_cellEdges);
        vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
        list_cellVertices(cell, num_cellEdges, cell_vrts);
        
        // Now that we have the vertices of the cells, we enumerate them, we get their coordinates
        std::map<int, int> correspondence; 
        std::vector<std::vector<double>> approx_cellVrts(cell_vrts.size(), {0,0,0});
        for (int j = 0; j < cell_vrts.size(); j++)
        {
            approx_cellVrts[j][0] = approxCoords[cell_vrts[j]][0];
            approx_cellVrts[j][1] = approxCoords[cell_vrts[j]][1];
            approx_cellVrts[j][2] = approxCoords[cell_vrts[j]][2];
            correspondence[cell_vrts[j]] = j;
        }


        // Then we are going to make the mean of the point using as weight the sum of the area of the face
        // touching the point
        std::vector<double> sum_area(num_cellVrts, 0);
        for (int f: cell.faces)
        {
            BSPface face = faces[f];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);
            
            for (int v : face_vrts)
            {
                sum_area[correspondence[v]] += face_info[f].area;
            }
        }

        auto centroid = barycenter_points(approx_cellVrts, sum_area);
        
        barycenters[i] = {centroid[0], centroid[1], centroid[2]};

    }

    std::cout << cells.size() << "\n\n\n\n" << std::endl;

    if (two_files)
    {
        markInternalCells_ray(BLACK_B, INTERNAL_B, barycenters, triangles, face_info, smoothing);
    }
    else markInternalCells_ray(BLACK_A, INTERNAL_A, barycenters, triangles,  face_info, smoothing);

}

#pragma endregion

// Second method using the topology of the mesh, but still shoot random ray
#pragma region

// Recursive function that count the number of intersections between the segment [a,b] and the complex
// once [a,b] as reach the cell cell_ind through the face face_ind
int BSPcomplex::count_intersection(genericPoint* a, genericPoint* b, int face_ind, int cell_ind, uint32_t skin_colour, Min_find &number_ray, std::vector<int>& inside_ray){
    int res = 0;
    if IS_GHOST_CELL(cell_ind) return 0;
    number_ray.increase(cell_ind);
    for( int f : cells[cell_ind].faces){
        if (f != face_ind){
            if (segment_Intersectface(f, a, b)){
                res += count_intersection(a, b, f, OTHER_ONE(faces[f].conn_cells, cell_ind), skin_colour, number_ray, inside_ray);
                if (isSkinFace(faces[f], skin_colour)) res += 1;
            }
        }
    }
    if (res%2){
        inside_ray[cell_ind] ++;
    }
    return res;
}

// Generate a point which is supposed to represent the infinite
genericPoint* generatePoint(){
    double x = random_float()-0.5, y = random_float()-0.5, z = random_float()-0.5;
    double norm = std::sqrt(x*x+y*y+z*z);
    return new explicitPoint3D(1000*x/norm, 1000*y/norm, 1000*z/norm);
}

void BSPcomplex::markInternalCells_ray(uint32_t skin_colour, uint32_t internal_label, std::vector<genericPoint*> barycenters, std::vector<Face_add>& face_info, int smoothing){

    // Get back the value of the minimum number of ray per cells,
    // passed in argument
    int N=10;
    if (global_parameters.size()>0) N = global_parameters[0];

    // We initialize the structure storing the number of ray in each cell
    // and gives in O(1) the number the cell with the less ray
    Min_find number_ray(cells.size());

    // Will store the number of ray saying inside for each cell
    std::vector<int> inside_ray(cells.size(), 0);

    int compteur = 0; // This value is not necessary and is just use to print more information

    // We send ray until each cell has been cross enough time
    while (number_ray.min < N){
        std::cout << compteur << " : " << number_ray.min << std::endl;
        compteur ++;
        genericPoint* outside = generatePoint();
        int i = number_ray.index_min();
        count_intersection(barycenters[i], outside, UINT32_MAX, i, skin_colour, number_ray, inside_ray);
        delete outside;
    }

    // We convert the data to the format required for the smoothing
    std::vector<int> inside(cells.size(), 0), outside(cells.size(), 0);

    for (int i = 0; i < cells.size(); i++)
    {
        inside[i] = inside_ray[i];
        outside[i] = number_ray.number_of_call(i) - inside_ray[i];
    }
    
    if(smoothing == 1){
        smoothing1(inside, outside, face_info, internal_label, skin_colour);
    }
    else if(smoothing == 2){
        smoothing2(inside, outside, face_info, internal_label, skin_colour);
    }
    else{
        for (int i = 0; i < cells.size(); i++)
        {
            if(inside[i] > outside[i]) setInternalCell(cells[i], internal_label);
        }
    }

    // ofstream f("number_of_ray_per_cell.hist");

    // for (int i = 0; i < cells.size(); i++)
    // {
    //     f << number_ray.number_of_call(i) << " ";
    // }
    // f.close();
    // int sum = 0;
    // for (int i = 0; i < cells.size(); i++)
    // {
    //     sum += number_ray.number_of_call(i);
    // }

    // std::cout << "Mean of number of ray : " << sum / cells.size() << std::endl;


}


void BSPcomplex::constraintsSurface_complexPartition_ray_exact(bool two_files, int smoothing){

    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);  

    std::vector<points> approxCoords(vertices.size(), points(3,0));
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i][0], approxCoords[i][1], approxCoords[i][2], true);


    // Precalculate approximate faces informations for use in markInternalCells
    std::vector<Face_add> face_info(faces.size());
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_info[i].area = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_info[i].area;

        std::vector<uint32_t> vs(faces[i].edges.size(), 0);
        list_faceVertices(faces[i], vs);
        face_info[i].vertex_ind = vs;
    }

    // Let's precompute an approximation of the barycenter of each cells to know from where shoot ray
    std::vector<genericPoint*> barycenters(cells.size());

    for (int i = 0; i < cells.size(); i++)
    {
        BSPcell cell = cells[i];
        uint64_t num_cellEdges = UINT64_MAX;
        uint32_t num_cellVrts = count_cellVertices(cell, &num_cellEdges);
        vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
        list_cellVertices(cell, num_cellEdges, cell_vrts);
        
        // Now that we have the vertices of the cells, we enumerate them, we get their coordinates
        std::map<int, int> correspondence; 
        std::vector<std::vector<double>> approx_cellVrts(cell_vrts.size(), {0,0,0});
        for (int j = 0; j < cell_vrts.size(); j++)
        {
            approx_cellVrts[j][0] = approxCoords[cell_vrts[j]][0];
            approx_cellVrts[j][1] = approxCoords[cell_vrts[j]][1];
            approx_cellVrts[j][2] = approxCoords[cell_vrts[j]][2];
            correspondence[cell_vrts[j]] = j;
        }

        // Then we are going to make the mean of the point using as weight the sum of the area of the face
        // touching the point
        std::vector<double> sum_area(num_cellVrts, 0);
        for (int f: cell.faces)
        {
            BSPface face = faces[f];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);
            
            for (int v : face_vrts)
            {
                sum_area[correspondence[v]] += face_info[f].area;
            }
        }

        auto centroid = barycenter_points(approx_cellVrts, sum_area);
        
        barycenters[i] = new explicitPoint3D(centroid[0], centroid[1], centroid[2]);

    }


    std::cout << cells.size() << "\n\n" << std::endl;

    if (two_files)
    {
        markInternalCells_ray(BLACK_A, INTERNAL_A, barycenters, face_info, smoothing);
        markInternalCells_ray(BLACK_B, INTERNAL_B, barycenters, face_info, smoothing);
    }
    else markInternalCells_ray(BLACK_A, INTERNAL_A, barycenters, face_info, smoothing);

    for (int i = 0; i < barycenters.size(); i++)
    {
        delete barycenters[i];
    }

    #ifdef Time_info
    std::cout << "Time to get intersections : " << intersection_time / 1000000000 << std::endl;
    #endif
}

#pragma endregion

// SORTING
#pragma region


std::vector<std::vector<std::vector<std::pair<int, genericPoint*>>>> BSPcomplex::classify_facets_exact(int axis, int sample, std::vector<Face_add> &face_info, std::vector<genericPoint*>& toBeDeleted){
    
    int axis_x = axis, axis_y = (axis+1)%3, axis_z = (axis+2)%3;

    // We are first going to compute a bounding box of the mesh
    std::vector<double> mini = {face_info[0].min[axis_x], face_info[0].min[axis_y], face_info[0].min[axis_z]}, maxi = {face_info[0].max[axis_x], face_info[0].max[axis_y], face_info[0].max[axis_z]};


    for (Face_add face : face_info)
    {
        mini[0] = std::min(face.min[axis_x], mini[0]);
        mini[1] = std::min(face.min[axis_y], mini[1]);
        mini[2] = std::min(face.min[axis_z], mini[2]);
        maxi[0] = std::max(face.max[axis_x], maxi[0]);
        maxi[1] = std::max(face.max[axis_y], maxi[1]);
        maxi[2] = std::max(face.max[axis_z], maxi[2]);
    }

    // We compute the space between two point of the grid
    double step_y = (maxi[1] - mini[1])/ ((float) sample-1), step_z = (maxi[2] - mini[2])/((float) sample-1);


    // auto start = explicitPoint3D(ROTATE_ARG2(mini[0]-1, mini[1]+1*step_y, mini[2]+2*step_z));
    // auto end = explicitPoint3D(ROTATE_ARG2(maxi[0]+1, mini[1]+1*step_y, mini[2]+2*step_z));
    // for (int i = 1; i < face_info[159].vertex_ind.size()-1; i++){
    //     std::cerr << segment_IntersectTriangle(&start, 
    //                                            &end, 
    //                                            vertices[face_info[159].vertex_ind[0]], 
    //                                            vertices[face_info[159].vertex_ind[i]], 
    //                                            vertices[face_info[159].vertex_ind[i+1]]) << std::endl;
    // }

    // std::cerr << segment_Intersectface(159, &start, &end);

    // exit(0);

    std::cout << step_y << " " << step_z << " " << std ::endl;
    print_vector(mini);
    print_vector(maxi);
    std::cout << std::endl;

    // Here we initialize the result, being the array representing the grid
    // and containg the list of the crossed faces and their intersection point
    std::vector<std::vector<std::vector<std::pair<int, genericPoint*>>>> res(sample, std::vector<std::vector<std::pair<int, genericPoint*>>>(sample));

    // Then we are going to add each face to all the tile of the grid the projection is in
    for (int f = 0; f < face_info.size(); f++)
    {        
        
        std::set<std::pair<int, int>> point_of_the_grid;
        auto vertex_coord = face_info[f].vertex_coord;
        double y_basis = vertex_coord[0][axis_y];
        double z_basis = vertex_coord[0][axis_z];

        for (int i = 1; i < face_info[f].vertex_ind.size()-1; i++)
        {
            // We project each triangle composing the face
            // We find a bounding rectangle containing this triangle
            int ind_y_min =  (std::min(y_basis, std::min(vertex_coord[i][axis_y], vertex_coord[i+1][axis_y])) -mini[1])/step_y;
            int ind_y_max =  (std::max(y_basis, std::max(vertex_coord[i][axis_y], vertex_coord[i+1][axis_y])) -mini[1])/step_y;
            int ind_z_min =  (std::min(z_basis, std::min(vertex_coord[i][axis_z], vertex_coord[i+1][axis_z])) -mini[2])/step_z;
            int ind_z_max =  (std::max(z_basis, std::max(vertex_coord[i][axis_z], vertex_coord[i+1][axis_z])) -mini[2])/step_z;
            
            // Then for each element of the grid also present in the bounding rectangle
            // we are going to test the intersection with the triangle, and if it succeeds
            // we will add it to the set of element of the grid where we should add the face
            for (int y = ind_y_min; y <= ind_y_max; y++)
            {
                for (int z = ind_z_min; z <= ind_z_max; z++)
                {
                    explicitPoint3D start, end;
                    switch (axis){
                        case 0:
                            start = explicitPoint3D(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z);
                            end = explicitPoint3D(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z);
                            break;
                        case 1:
                            start = explicitPoint3D(ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z));
                            end = explicitPoint3D(ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z));
                            break;
                        case 2:
                            start = explicitPoint3D(ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z));
                            end = explicitPoint3D(ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z));
                            break;
                    }
                    if(segment_IntersectTriangle(&start, 
                                                 &end, 
                                                 vertices[face_info[f].vertex_ind[0]], 
                                                 vertices[face_info[f].vertex_ind[i]], 
                                                 vertices[face_info[f].vertex_ind[i+1]])){
                        
                        point_of_the_grid.insert({y, z});

                    }
                    
                }
            }
        }

        // Then we actually add the ray to the grid
        for (std::pair<int, int> pair : point_of_the_grid){
            int y = pair.first, z = pair.second;
            explicitPoint3D *start, *end;
            switch (axis){
                case 0:
                    start = new explicitPoint3D(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z);
                    end = new explicitPoint3D(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z);
                    break;
                case 1:
                    start = new explicitPoint3D(ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z));
                    end = new explicitPoint3D(ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z));
                    break;
                case 2:
                    start = new explicitPoint3D(ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z));
                    end = new explicitPoint3D(ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z));
                    break;
            }
            BSPface& face = faces[f];
            implicitPoint3D_LPI *new_point = new implicitPoint3D_LPI(*start, *end, vertices[face.meshVertices[0]]->toExplicit3D(),  vertices[face.meshVertices[1]]->toExplicit3D(), vertices[face.meshVertices[2]]->toExplicit3D());
            res[y][z].push_back(std::pair<int, genericPoint*>(f, new_point));
            
            toBeDeleted.push_back(start); toBeDeleted.push_back(end); toBeDeleted.push_back(new_point);
        }
/*       
        int ind_y_min = (face_info[f].min[axis_y]-mini[1])/step_y;
        int ind_y_max = (face_info[f].max[axis_y]-mini[1])/step_y;
        int ind_z_min = (face_info[f].min[axis_z]-mini[2])/step_z;
        int ind_z_max = (face_info[f].max[axis_z]-mini[2])/step_z;
        
        for (int y = ind_y_min; y <= ind_y_max; y++)
        {
            for (int z = ind_z_min; z <= ind_z_max; z++)
            {
                explicitPoint3D *start, *end;
                switch (axis){
                    case 0:
                        start = new explicitPoint3D(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z);
                        end = new explicitPoint3D(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z);
                        break;
                    case 1:
                        start = new explicitPoint3D(ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z));
                        end = new explicitPoint3D(ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z));
                        break;
                    case 2:
                        start = new explicitPoint3D(ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z));
                        end = new explicitPoint3D(ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z));
                        break;
                }
                if(segment_Intersectface(f, start, end)){
                    BSPface& face = faces[f];
                    implicitPoint3D_LPI *new_point = new implicitPoint3D_LPI(*start, *end, vertices[face.meshVertices[0]]->toExplicit3D(),  vertices[face.meshVertices[1]]->toExplicit3D(), vertices[face.meshVertices[2]]->toExplicit3D());
                    res[y][z].push_back(std::pair<int, genericPoint*>(f, new_point));
                    
                    // double _x, _y, _z;
                    // std::get<1>(res[y][z].at(res[y][z].size()-1))->getApproxXYZCoordinates(_x, _y, _z);
                    // std::cout << _x << " " << _y << " " << _z << std::endl;
                    toBeDeleted.push_back(new_point);
                toBeDeleted.push_back(start); toBeDeleted.push_back(end);
                }
            }
        }
*/
    }

    return res;
}

void BSPcomplex::shoot_axis_ray_exact(int axis, int sample, std::vector<int>& inside, std::vector<int>& outside, std::vector<Face_add> &face_info, uint32_t skin_colour){

    std::vector<genericPoint*> toBeDeleted(0);
    auto classification = classify_facets_exact(axis, sample, face_info, toBeDeleted);

    // for(int a = 0; a < classification.size(); a++){
    //     for(int b = 0; b < classification[a].size(); b++){
    //         std::cout << a << " " << b << " : ";
    //         for(auto p: classification[a][b]){
    //             std::cout << p.first << " ";
    //         }
    //         std::cout << "\n";
    //     }
    // }

    // We are going to treat separetely each cell of the grid
    for (int i = 0; i < classification.size(); i++)
    {
        for (int j = 0; j < classification[i].size(); j++)
        {
            // We are going to sort the intersection points


            // std::cout << "stack : "<< i << " " << j << " of size " << classification[i][j].size() << " on axis " << (int) axis << std::endl;
            switch (axis)
            {
            case 0:
                std::sort(classification[i][j].begin(), classification[i][j].end(), 
                    [](const std::pair<int, genericPoint*>& a, const std::pair<int, genericPoint*>& b){
                        return (genericPoint::lessThanOnX(*std::get<1>(a), *std::get<1>(b)) == 1);
                    });
                break;
            case 1:
                std::sort(classification[i][j].begin(), classification[i][j].end(), 
                    [](const std::pair<int, genericPoint*>& a, const std::pair<int, genericPoint*>& b){
                        return (genericPoint::lessThanOnY(*std::get<1>(a), *std::get<1>(b)) == 1);
                    });
                break;
            case 2:
                std::sort(classification[i][j].begin(), classification[i][j].end(), 
                    [](const std::pair<int, genericPoint*>& a, const std::pair<int, genericPoint*>& b){
                        return (genericPoint::lessThanOnZ(*std::get<1>(a), *std::get<1>(b)) == 1);
                    });
                break;
            }


            // std::cerr << i << " " << j << " : " << std::endl;
            // for (int k = 0; k < classification[i][j].size(); k++)
            // {
            //     points intersec = {0,0,0}; classification[i][j][k].second->getApproxXYZCoordinates(intersec[0], intersec[1], intersec[2]);
            //     std::cerr << i << " " << j << " : " << classification[i][j][k].first << " " << intersec[axis] << std::endl;
            //     //std::cerr << classification[i][j][k].first << " " << faces[classification[i][j][k].first].conn_cells[0] << " " << faces[classification[i][j][k].first].conn_cells[1] << std::endl;
            // }
            // std::cerr << std::endl;


            // Then we are going to use the list of the crossed face to determine 
            // what the ray is saying about the crossed cells
            bool is_outside = true;

            int actual_cell = GHOST_CELL;

            std::vector<std::pair<int, bool>> to_be_change(0);

            for (int k = 0; k < classification[i][j].size(); k++)
            {
                
                BSPface face = faces[std::get<0>(classification[i][j][k])];
                if (!IS_IN(face.conn_cells, actual_cell)){
                    // We deal with the case where, beacause of rounding errors 
                    // or because we shoot a ray through a point, 
                    // we don't have a valid path along the complex
                    // (detected because two consecutive faces don't share a cell)
                    std::cout << "A ray failed " << i << " " << j << " at the face " << k << std::endl;
                    break;
                }
                if (isSkinFace(face, skin_colour)) is_outside = !is_outside;
                actual_cell = OTHER_ONE(face.conn_cells, actual_cell);

               if(!IS_GHOST_CELL(actual_cell)){
                    to_be_change.push_back(std::pair<int, bool>(actual_cell, is_outside));
                }
            }

            // We are going to update the information of the ray, only if at the end, we end up on outside, (which is supposed to be the case
            // as we are supposed to end up on the ghost cell which is outside), because otherwise it means that we cross an odd number of hole
            // (thus not 0)
            if (is_outside){
                for (auto &&p : to_be_change)
                {
                    if (p.second) outside[p.first]++;
                    else inside[p.first]++;
                }
            }
        }
    }

    // for (auto points : classification[9][8]){
    //     std::cerr << points.first << " " << faces[points.first].conn_cells[0] << " " << faces[points.first].conn_cells[1] << std::endl;
    // }

    for (auto p : toBeDeleted){
        delete p;
    }

}

void BSPcomplex::markInternalCells_exact(uint32_t skin_colour, uint32_t internal_label, std::vector<Face_add>& face_info, int smoothing){
    std::vector<int> outside(cells.size()), inside(cells.size());

    // We get the density of the grid we are going to use through the parameter
    int sample = std::cbrt(cells.size())+1;
    if (global_parameters.size()>0) sample = global_parameters[0]*sample;

    // We get the number of directions that we are going to use 
    // (knowing that for eahc rotation we have three directions)
    int number_of_rotation = 5;
    if (global_parameters.size()>1) number_of_rotation = global_parameters[1];


    for(int number_rotation = 0; number_rotation < number_of_rotation; number_rotation++){

        //Do a rotation
        double teta1 = random_float()*5, teta2 = random_float()*5, teta3 = random_float()*5;
        std::vector<std::vector<double>> matrix = rotation_matrix(teta1, teta2, teta3);
        for( auto point : explicit_vertices){
            point->rotation(matrix);
        }

        // Precalculate approximate vertex coordinates for use in approxFaceArea()
        std::vector<std::vector<double>> approxCoords(vertices.size(), {0,0,0});
        for (size_t i = 0; i < vertices.size(); i++)
            vertices[i]->getApproxXYZCoordinates(approxCoords[i][0], approxCoords[i][1], approxCoords[i][2]);


        // We compute some faces informations to use them in shoot_axis_ray_exact
        for (int f = 0; f < faces.size(); f++)
        {
            std::vector<std::vector<double>> face_coord(face_info[f].vertex_ind.size(), {0,0,0});

            for (int i = 0; i < face_coord.size(); i++)
            {
                face_coord[i] = approxCoords[face_info[f].vertex_ind[i]];
            }
            face_info[f].vertex_coord = face_coord;
            face_info[f].vertex_ind = face_info[f].vertex_ind;

            face_info[f].min = {face_coord[0][0], face_coord[0][1], face_coord[0][2]};
            face_info[f].max = face_info[f].min;

            for (int i = 0; i < face_coord.size(); i++)
            {
                face_info[f].min[0] = std::min(face_info[f].min[0], face_coord[i][0]);
                face_info[f].min[1] = std::min(face_info[f].min[1], face_coord[i][1]);
                face_info[f].min[2] = std::min(face_info[f].min[2], face_coord[i][2]);
                face_info[f].max[0] = std::max(face_info[f].max[0], face_coord[i][0]);
                face_info[f].max[1] = std::max(face_info[f].max[1], face_coord[i][1]);
                face_info[f].max[2] = std::max(face_info[f].max[2], face_coord[i][2]);
            }
            
        }

        shoot_axis_ray_exact(0, sample, inside, outside, face_info, skin_colour);
        shoot_axis_ray_exact(1, sample, inside, outside, face_info, skin_colour);
        shoot_axis_ray_exact(2, sample, inside, outside, face_info, skin_colour);

        // We undo the rotation
        matrix = inverse_rotation_matrix(teta1, teta2, teta3);
        for( auto point : explicit_vertices){
            point->rotation(matrix);
        }

    }


    std::cout << "Operation on the grid done.\nStart complementary actions"<< std::endl;
    
    if(smoothing == 1){
        smoothing1(inside, outside, face_info, internal_label, skin_colour);
    }
    else if(smoothing == 2){
        smoothing2(inside, outside, face_info, internal_label, skin_colour);
    }
    else{
        for (int i = 0; i < cells.size(); i++)
        {
            if(inside[i] > outside[i]) setInternalCell(cells[i], internal_label);
        }
    }
}

void BSPcomplex::constraintsSurface_complexPartition_grid_exact(bool two_files, int smoothing){

    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);  

    std::vector<points> approxCoords(vertices.size(), points(3,0));
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i][0], approxCoords[i][1], approxCoords[i][2]);


    // Precalculate approximate face informations for use in markInternalCells
    std::vector<Face_add> face_info(faces.size());
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_info[i].area = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_info[i].area;

        std::vector<uint32_t> vs(faces[i].edges.size(), 0);
        list_faceVertices(faces[i], vs);
        face_info[i].vertex_ind = vs;
    }

    if (two_files)
    {
        markInternalCells_exact(BLACK_A, INTERNAL_A, face_info, smoothing);
        markInternalCells_exact(BLACK_B, INTERNAL_B, face_info, smoothing);
    }
    else {markInternalCells_exact(BLACK_A, INTERNAL_A, face_info, smoothing);}

}


std::vector<std::vector<std::vector<std::pair<int, double>>>> BSPcomplex::classify_facets_approx(int axis, int sample, std::vector<Face_add>& face_info, std::vector<points>& vertex_coord){
    
    int axis_x = axis, axis_y = (axis+1)%3, axis_z = (axis+2)%3;

    // We are first going to compute a bounding box of the mesh
    std::vector<double> mini = {vertex_coord[0][axis_x], vertex_coord[0][axis_y], vertex_coord[0][axis_z]}, maxi = {vertex_coord[0][axis_x], vertex_coord[0][axis_y], vertex_coord[0][axis_z]};


    for (points p : vertex_coord)
    {
        mini[0] = std::min(p[axis_x], mini[0]);
        mini[1] = std::min(p[axis_y], mini[1]);
        mini[2] = std::min(p[axis_z], mini[2]);
        maxi[0] = std::max(p[axis_x], maxi[0]);
        maxi[1] = std::max(p[axis_y], maxi[1]);
        maxi[2] = std::max(p[axis_z], maxi[2]);
    }

    // We compute the space between two point of the grid
    double step_y = (maxi[1] - mini[1])/ ((float) sample-1), step_z = (maxi[2] - mini[2])/((float) sample-1);

    std::cout << step_y << " " << step_z << " " << std ::endl;
    print_vector(mini);
    print_vector(maxi);
    std::cout << std::endl;

    // Here we initialize the result, being the array representing the grid
    // and containg the list of the crossed faces and their intersection point
    std::vector<std::vector<std::vector<std::pair<int, double>>>> res(sample, std::vector<std::vector<std::pair<int, double>>>(sample));


    // Then we are going to add each face to all the tile of the grid the projection is in
    for (int f = 0; f < faces.size(); f++)
    {        
        
        std::set<std::pair<int, int>> point_of_the_grid;
        auto vertex_ind_f = face_info[f].vertex_ind;
        double y_basis = vertex_coord[vertex_ind_f[0]][axis_y];
        double z_basis = vertex_coord[vertex_ind_f[0]][axis_z];

        for (int i = 1; i < vertex_ind_f.size()-1; i++)
        {
            // We project each triangle composing the face
            // We find a bounding rectangle containing this triangle
            int ind_y_min =  (std::min(y_basis, std::min(vertex_coord[vertex_ind_f[i]][axis_y], vertex_coord[vertex_ind_f[i+1]][axis_y])) -mini[1])/step_y;
            int ind_y_max =  (std::max(y_basis, std::max(vertex_coord[vertex_ind_f[i]][axis_y], vertex_coord[vertex_ind_f[i+1]][axis_y])) -mini[1])/step_y;
            int ind_z_min =  (std::min(z_basis, std::min(vertex_coord[vertex_ind_f[i]][axis_z], vertex_coord[vertex_ind_f[i+1]][axis_z])) -mini[2])/step_z;
            int ind_z_max =  (std::max(z_basis, std::max(vertex_coord[vertex_ind_f[i]][axis_z], vertex_coord[vertex_ind_f[i+1]][axis_z])) -mini[2])/step_z;
            

            // Then for each element of the grid also present in the bounding rectangle
            // we are going to test the intersection with the triangle, and if it succeeds
            // we will add it to the set of element of the grid where we should add the face
            for (int y = ind_y_min; y <= ind_y_max; y++)
            {
                for (int z = ind_z_min; z <= ind_z_max; z++)
                {
                    points start, end;
                    switch (axis){
                        case 0:
                            start = {mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z};
                            end = {maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z};
                            break;
                        case 1:
                            start = {ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
                        case 2:
                            start = {ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
                    }
                    if(segment_IntersectTriangle(start, 
                                                 end, 
                                                 vertex_coord[vertex_ind_f[0]], 
                                                 vertex_coord[vertex_ind_f[i]], 
                                                 vertex_coord[vertex_ind_f[i+1]])){
                        point_of_the_grid.insert({y, z});
                    }
                    
                }
            }
        }

        // Then we actually add the ray to the grid
        for (std::pair<int, int> pair : point_of_the_grid){
            int y = pair.first, z = pair.second;
            points start, end;
                    switch (axis){
                        case 0:
                            start = {mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z};
                            end = {maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z};
                            break;
                        case 1:
                            start = {ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
                        case 2:
                            start = {ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
            }
            BSPface& face = faces[f];
            points new_point =  intersection_line_plane(start, end, vertex_coord[face.meshVertices[0]],  vertex_coord[face.meshVertices[1]], vertex_coord[face.meshVertices[2]]);
            res[y][z].push_back(std::pair<int, double>(f, new_point[axis]));
            
        }
    }

    return res;

}

void BSPcomplex::shoot_axis_ray_approx(int axis, int sample, std::vector<int>& inside, std::vector<int>& outside, std::vector<Face_add> &face_info, std::vector<points>& vertex_coord, uint32_t skin_colour){
    
    auto classification = classify_facets_approx(axis, sample, face_info, vertex_coord);


    for (int i = 0; i < classification.size(); i++)
    {
        for (int j = 0; j < classification[i].size(); j++)
        {
            // std::cout << "stack : "<< i << " " << j << " of size " << classification[i][j].size() << " on axis " << (int) axis << std::endl;
            std::sort(classification[i][j].begin(), classification[i][j].end(), 
                [](const std::pair<int, double>& a, const std::pair<int, double>& b){
                    return (a.second < b.second);
                });


            // std::cerr << i << " " << j << " : " << std::endl;
            // for (int k = 0; k < classification[i][j].size(); k++)
            // {
            //     std::cerr << i << " " << j << " : " << classification[i][j][k].first << " " << classification[i][j][k].second << std::endl;
            //     //std::cerr << classification[i][j][k].first << " " << faces[classification[i][j][k].first].conn_cells[0] << " " << faces[classification[i][j][k].first].conn_cells[1] << std::endl;
            // }
            // std::cerr << std::endl;


            bool is_outside = true;

            int actual_cell = GHOST_CELL;

            std::vector<std::pair<int, bool>> to_be_change(0);

            for (int k = 0; k < classification[i][j].size(); k++)
            {
                
                BSPface face = faces[std::get<0>(classification[i][j][k])];
                if (!IS_IN(face.conn_cells, actual_cell)){
                    // We deal with the case where, beacause of rounding errors 
                    // or because we shoot a ray through a point, 
                    // we don't have a valid path along the complex
                    // (detected because two consecutive faces don't share a cell)
                    std::cout << "A ray failed " << i << " " << j << " at the face " << k << std::endl;
                    break;
                }
                if (isSkinFace(face, skin_colour)) is_outside = !is_outside;
                actual_cell = OTHER_ONE(face.conn_cells, actual_cell);

                if(!IS_GHOST_CELL(actual_cell)){
                    to_be_change.push_back(std::pair<int, bool>(actual_cell, is_outside));
                }
            }

            // We are going to update the information of the ray, only if at the end, we end up on outside, (which is supposed to be the case
            // as we are supposed to end up on the ghost cell which is outside), because otherwise it means that we cross an odd number of hole
            // (thus not 0)
            if (is_outside){
                for (auto &&p : to_be_change)
                {
                    if (p.second) outside[p.first]++;
                    else inside[p.first]++;
                }
                
            }
        }
    }

}

void BSPcomplex::markInternalCells_approx(uint32_t skin_colour, uint32_t internal_label, std::vector<Face_add>& face_info, std::vector<points> vertex_coord_init, int smoothing){
    
    std::vector<int> outside(cells.size()), inside(cells.size());

    int sample = std::cbrt(cells.size())+1;
    if (global_parameters.size()>0) sample = global_parameters[0]*sample;


    int number_of_rotation = 5;
    if (global_parameters.size()>1) number_of_rotation = global_parameters[1];

    for(int number_rotation = 0; number_rotation < number_of_rotation; number_rotation++){
        
        std::vector<points> vertex_coord = vertex_coord_init;

        //Do a rotation
        double teta1 = random_float()*5, teta2 = random_float()*5, teta3 = random_float()*5;
        std::vector<std::vector<double>> matrix = rotation_matrix(teta1, teta2, teta3);
        for( points& point : vertex_coord){
            point = matrix_product(matrix, point);
        }

        shoot_axis_ray_approx(0, sample, inside, outside, face_info, vertex_coord, skin_colour);
        shoot_axis_ray_approx(1, sample, inside, outside, face_info, vertex_coord, skin_colour);
        shoot_axis_ray_approx(2, sample, inside, outside, face_info, vertex_coord, skin_colour);


    }


    std::cout << "Operation on the grid done.\nStart complementary actions"<< std::endl;


    if(smoothing == 1){
        smoothing1(inside, outside, face_info, internal_label, skin_colour);
    }
    else if(smoothing == 2){
        smoothing2(inside, outside, face_info, internal_label, skin_colour);
    }
    else{
        for (int i = 0; i < cells.size(); i++)
        {
            if(inside[i] > outside[i]) setInternalCell(cells[i], internal_label);
        }
    }
}

void BSPcomplex::constraintsSurface_complexPartition_grid_approx(bool two_files, int smoothing){

    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);  

    std::vector<points> approxCoords(vertices.size(), points(3,0));
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i][0], approxCoords[i][1], approxCoords[i][2], true);


    // Precalculate approximate face informations for use in markInternalCells
    std::vector<Face_add> face_info(faces.size());
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_info[i].area = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_info[i].area;

        std::vector<uint32_t> vs(faces[i].edges.size(), 0);
        list_faceVertices(faces[i], vs);
        face_info[i].vertex_ind = vs;
    }

    if (two_files)
    {
        markInternalCells_approx(BLACK_A, INTERNAL_A, face_info, approxCoords, smoothing);
        markInternalCells_approx(BLACK_B, INTERNAL_B, face_info, approxCoords, smoothing);
    }
    else {markInternalCells_approx(BLACK_A, INTERNAL_A, face_info, approxCoords, smoothing);}

}



std::vector<std::vector<std::vector<std::pair<int, double>>>> BSPcomplex::classify_facets_new(int axis, int sample, std::vector<Face_add>& face_info, std::vector<points>& vertex_coord, double& area_sample){
    
    int axis_x = axis, axis_y = (axis+1)%3, axis_z = (axis+2)%3;

    // We are first going to compute a bounding box of the mesh
    std::vector<double> mini = {vertex_coord[0][axis_x], vertex_coord[0][axis_y], vertex_coord[0][axis_z]}, maxi = {vertex_coord[0][axis_x], vertex_coord[0][axis_y], vertex_coord[0][axis_z]};



    for (points p : vertex_coord)
    {
        mini[0] = std::min(p[axis_x], mini[0]);
        mini[1] = std::min(p[axis_y], mini[1]);
        mini[2] = std::min(p[axis_z], mini[2]);
        maxi[0] = std::max(p[axis_x], maxi[0]);
        maxi[1] = std::max(p[axis_y], maxi[1]);
        maxi[2] = std::max(p[axis_z], maxi[2]);
    }

    // We compute the space between two point of the grid
    double step_y = (maxi[1] - mini[1])/ ((float) sample-1), step_z = (maxi[2] - mini[2])/((float) sample-1);

    area_sample = step_y*step_z;


    // Here we initialize the result, being the array representing the grid
    // and containg the list of the crossed faces and their intersection point
    std::vector<std::vector<std::vector<std::pair<int, double>>>> res(sample, std::vector<std::vector<std::pair<int, double>>>(sample));


    // Then we are going to add each face to all the tile of the grid the projection is in
    for (int f = 0; f < faces.size(); f++)
    {        
        
        std::set<std::pair<int, int>> point_of_the_grid;
        auto vertex_ind_f = face_info[f].vertex_ind;
        double y_basis = vertex_coord[vertex_ind_f[0]][axis_y];
        double z_basis = vertex_coord[vertex_ind_f[0]][axis_z];

        for (int i = 1; i < vertex_ind_f.size()-1; i++)
        {
            // We project each triangle composing the face
            // We find a bounding rectangle containing this triangle
            int ind_y_min =  (std::min(y_basis, std::min(vertex_coord[vertex_ind_f[i]][axis_y], vertex_coord[vertex_ind_f[i+1]][axis_y])) -mini[1])/step_y;
            int ind_y_max =  (std::max(y_basis, std::max(vertex_coord[vertex_ind_f[i]][axis_y], vertex_coord[vertex_ind_f[i+1]][axis_y])) -mini[1])/step_y;
            int ind_z_min =  (std::min(z_basis, std::min(vertex_coord[vertex_ind_f[i]][axis_z], vertex_coord[vertex_ind_f[i+1]][axis_z])) -mini[2])/step_z;
            int ind_z_max =  (std::max(z_basis, std::max(vertex_coord[vertex_ind_f[i]][axis_z], vertex_coord[vertex_ind_f[i+1]][axis_z])) -mini[2])/step_z;
            

            // Then for each element of the grid also present in the bounding rectangle
            // we are going to test the intersection with the triangle, and if it succeeds
            // we will add it to the set of element of the grid where we should add the face
            for (int y = ind_y_min; y <= ind_y_max; y++)
            {
                for (int z = ind_z_min; z <= ind_z_max; z++)
                {
                    points start, end;
                    switch (axis){
                        case 0:
                            start = {mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z};
                            end = {maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z};
                            break;
                        case 1:
                            start = {ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
                        case 2:
                            start = {ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
                    }
                    if(segment_IntersectTriangle(start, 
                                                 end, 
                                                 vertex_coord[vertex_ind_f[0]], 
                                                 vertex_coord[vertex_ind_f[i]], 
                                                 vertex_coord[vertex_ind_f[i+1]])){
                        point_of_the_grid.insert({y, z});
                    }
                    
                }
            }
        }

        // Then we actually add the ray to the grid
        for (std::pair<int, int> pair : point_of_the_grid){
            int y = pair.first, z = pair.second;
            points start, end;
                    switch (axis){
                        case 0:
                            start = {mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z};
                            end = {maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z};
                            break;
                        case 1:
                            start = {ROTATE_ARG1(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG1(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
                        case 2:
                            start = {ROTATE_ARG2(mini[0]-1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            end = {ROTATE_ARG2(maxi[0]+1, mini[1]+y*step_y, mini[2]+z*step_z)};
                            break;
            }
            BSPface& face = faces[f];
            points new_point =  intersection_line_plane(start, end, vertex_coord[face.meshVertices[0]],  vertex_coord[face.meshVertices[1]], vertex_coord[face.meshVertices[2]]);
            res[y][z].push_back(std::pair<int, double>(f, new_point[axis]));
            
        }
    }

    return res;

}

double BSPcomplex::shoot_axis_ray_new(int axis, int sample, std::vector<int>& inside, std::vector<int>& outside, std::vector<int>& sum_ray, std::vector<Face_add> &face_info, std::vector<points>& vertex_coord, uint32_t skin_colour){
    
    double area_sample;
    auto classification = classify_facets_new(axis, sample, face_info, vertex_coord, area_sample);


    for (int i = 0; i < classification.size(); i++)
    {
        for (int j = 0; j < classification[i].size(); j++)
        {
            // std::cout << "stack : "<< i << " " << j << " of size " << classification[i][j].size() << " on axis " << (int) axis << std::endl;
            std::sort(classification[i][j].begin(), classification[i][j].end(), 
                [](const std::pair<int, double>& a, const std::pair<int, double>& b){
                    return (a.second < b.second);
                });


            // std::cout << i << " " << j << " : " << std::endl;
            // for (int k = 0; k < classification[i][j].size(); k++)
            // {
            //     // std::cout << i << " " << j << " : " << classification[i][j][k].first << " " << classification[i][j][k].second << std::endl;
            //     std::cout << " : " << classification[i][j][k].first << " " << classification[i][j][k].second << " " << faces[classification[i][j][k].first].conn_cells[0] << " " << faces[classification[i][j][k].first].conn_cells[1] << std::endl;
            // }
            // std::cout << std::endl;


            bool is_outside = true;

            int actual_cell = GHOST_CELL;

            std::vector<std::pair<int, bool>> to_be_change(0);

            for (int k = 0; k < classification[i][j].size(); k++)
            {
                
                BSPface face = faces[std::get<0>(classification[i][j][k])];
                if (!IS_IN(face.conn_cells, actual_cell)){
                    // We deal with the case where, beacause of rounding errors 
                    // or because we shoot a ray through a point, 
                    // we don't have a valid path along the complex
                    // (detected because two consecutive faces don't share a cell)
                    std::cout << "A ray failed " << i << " " << j << " at the " << k << "-th face"<< std::endl;
                    break;
                }
                if (isSkinFace(face, skin_colour)) is_outside = !is_outside;
                actual_cell = OTHER_ONE(face.conn_cells, actual_cell);

                if(!IS_GHOST_CELL(actual_cell)){
                    to_be_change.push_back(std::pair<int, bool>(actual_cell, is_outside));
                    sum_ray[actual_cell]++;
                }
            }

            // We are going to update the information of the ray, only if at the end, we end up on outside, (which is supposed to be the case
            // as we are supposed to end up on the ghost cell which is outside), because otherwise it means that we cross an odd number of hole
            // (thus not 0)
            if (is_outside){
                for (auto &&p : to_be_change)
                {
                    if (p.second) outside[p.first]++;
                    else inside[p.first]++;
                }
                
            }

        }
    }

    return area_sample;
}

void BSPcomplex::markInternalCells_new(uint32_t skin_colour, uint32_t internal_label, std::vector<Face_add>& face_info, std::vector<points> vertex_coord_init, int smoothing){
    
    int sample = std::cbrt(cells.size())+1;
    if (global_parameters.size()>0) sample = global_parameters[0]*sample;
    std::cout << sample << "\n" << std::endl;;

    int number_of_rotation = 5;
    if (global_parameters.size()>1) number_of_rotation = global_parameters[1];

    std::vector<std::vector<int>> outside_vect(0), inside_vect(0);
    std::vector<int> sum_ray(cells.size(), 0);

    std::vector<double> area_samples(0);

    for(int number_rotation = 0; number_rotation < number_of_rotation; number_rotation++){
        
        std::vector<points> vertex_coord = vertex_coord_init;

        //Do a rotation
        double teta1 = random_float()*5, teta2 = random_float()*5, teta3 = random_float()*5;
        std::vector<std::vector<double>> matrix = rotation_matrix(teta1, teta2, teta3);

        print_vector(matrix); std::cout <<"\n";

        for( points& point : vertex_coord){
            point = matrix_product(matrix, point);
        }

        outside_vect.push_back(std::vector<int>(cells.size(),0));
        inside_vect.push_back(std::vector<int>(cells.size(),0));
        area_samples.push_back(shoot_axis_ray_new(0, sample, inside_vect[inside_vect.size()-1], outside_vect[outside_vect.size()-1], sum_ray, face_info, vertex_coord, skin_colour));
        
        outside_vect.push_back(std::vector<int>(cells.size(),0));
        inside_vect.push_back(std::vector<int>(cells.size(),0));
        area_samples.push_back(shoot_axis_ray_new(1, sample, inside_vect[inside_vect.size()-1], outside_vect[outside_vect.size()-1], sum_ray, face_info, vertex_coord, skin_colour));
        
        outside_vect.push_back(std::vector<int>(cells.size(),0));
        inside_vect.push_back(std::vector<int>(cells.size(),0));
        area_samples.push_back(shoot_axis_ray_new(2, sample, inside_vect[inside_vect.size()-1], outside_vect[outside_vect.size()-1], sum_ray, face_info, vertex_coord, skin_colour));

    }

    double area_sample = 0;
    for (double a : area_samples) area_sample += a;
    area_sample /= area_samples.size();

    std::cout << "Operation on the grid done.\nStart complementary actions"<< std::endl;

    std::vector<int> sum_inside(cells.size()), sum_outside(cells.size());

    for (int i = 0; i < inside_vect.size(); i++)
    {
        for (int j = 0; j < cells.size(); j++)
        {
            sum_inside[j] += inside_vect[i][j];
            sum_outside[j] += outside_vect[i][j];
        }
    }

    std::vector<int> inside(cells.size(), 0), outside(cells.size(), 0);

    for (int i = 0; i < inside_vect.size(); i++)
    {
        for (int j = 0; j < cells.size(); j++)
        {
            inside[j] += inside_vect[i][j] / ((double)(inside_vect[i][j] + outside_vect[i][j] + 1)) * (sum_ray[j]);
            outside[j] += outside_vect[i][j] / ((double)(outside_vect[i][j] + inside_vect[i][j] + 1)) * (sum_ray[j]);
        }
    }
    
    if(smoothing == 1){
        smoothing1(inside, outside, face_info, internal_label, skin_colour);
    }
    else if(smoothing == 2){
        smoothing2(inside, outside, face_info, internal_label, skin_colour);
    }
    else if (smoothing == 3){
        smoothing3(inside, outside, area_sample/3/number_of_rotation, face_info, internal_label, skin_colour);
    }
    else{
        for (int i = 0; i < cells.size(); i++)
        {
            if (inside[i]>outside[i]){
                setInternalCell(cells[i], internal_label);
            }
        }
    }
}

void BSPcomplex::constraintsSurface_complexPartition_grid_new(bool two_files, int smoothing){

    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);  

    std::vector<points> approxCoords(vertices.size(), points(3,0));
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(approxCoords[i][0], approxCoords[i][1], approxCoords[i][2], true);




    // Precalculate approximate face informations for use in markInternalCells
    std::vector<Face_add> face_info(faces.size());
    double tot_face_area = 0.0;
    for (size_t i = 0; i < faces.size(); i++)
    {
        face_info[i].area = approxFaceArea(faces[i], this, approxCoords);
        tot_face_area += face_info[i].area;

        std::vector<uint32_t> vs(faces[i].edges.size(), 0);
        list_faceVertices(faces[i], vs);
        face_info[i].vertex_ind = vs;
    }


    if (two_files)
    {
        markInternalCells_new(BLACK_A, INTERNAL_A, face_info, approxCoords, smoothing);
        markInternalCells_new(BLACK_B, INTERNAL_B, face_info, approxCoords, smoothing);
    }
    else {markInternalCells_new(BLACK_A, INTERNAL_A, face_info, approxCoords, smoothing);}

}



#pragma endregion
