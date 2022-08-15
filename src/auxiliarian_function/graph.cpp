#include "graph.h"

// Shortest_path_graph
#pragma region 
Shortest_path_graph::Shortest_path_graph(std::vector<std::vector<double>> barycenters){
    for (int i = 0; i < barycenters.size(); i++)
    {
        old_to_new[i] = i;
        new_to_old.push_back(i);
        position.push_back(barycenters[i]);
    }
    
}

void Shortest_path_graph::add_node(int v, std::vector<double> barycenter){
    if (!old_to_new.count(v)){
        old_to_new[v] = new_to_old.size();
        new_to_old.push_back(v);
        position.push_back(barycenter);
    }
}

void Shortest_path_graph::add_edge(int u, int v){

    int a = old_to_new[u], b = old_to_new[v];
    edges.insert({a,b});
    edges.insert({b,a});
}


void Shortest_path_graph::print(){
    std::cout << "Nodes : "; print_vector(new_to_old);
    std::cout << "\n Edges :\n"; print_set(edges);
    std::cout << "\n Positions : \n"; print_vector(position);
}


inline double Shortest_path_graph::weight(int a, int b){
    return euclidean_distance(position[a], position[b]);
}

int find_minimum(std::set<int> Q, std::vector<double>& key){
    if (Q.empty()) return -1;
    double minimum = 1./0.;
    int res = -1;
    for (int x: Q){
        if (key[x] < minimum){
            res = x;
            minimum = key[x];
        }
    }
    return res;
}

std::map<int, double> Shortest_path_graph::djikstra(int source){

    source = old_to_new[source];
    
    int N = new_to_old.size();
    std::vector<std::set<int>> adjacence(N);

    for (auto edge : edges){
        adjacence[edge[0]].insert(edge[1]);
    }
    
    
    double infinity = 1./0.;
    std::vector<double> dist(N, infinity);

    dist[source] = 0;

    std::set<int> Q;
    for (int i = 0; i < N; i++) Q.insert(i);
    
    int compteur = 0;
    while(!Q.empty()){
        std::cout << compteur << std::endl;
        compteur ++;
        int u = find_minimum(Q, dist);
        if (u == -1) break; // It happens when the graph is not connected
        Q.erase(u);
        for (int v: adjacence[u])
        {
            double alt = dist[u] + weight(u,v);
            if(alt < dist[v]) dist[v] = alt;
        }
    }
    
    // Now we will get back the old label to have a real correspondance
    

    std::map<int, double> res;
    for (int i = 0; i < N; i++)
    {
        res[new_to_old[i]] = dist[i];
    }

    return res;
}



// Save the graph at a .obj file
void Shortest_path_graph::save(const char* fileName){
    std::ofstream f(fileName);
    
    if (!f) { std::cout<< "Can't save the vector because the file cannot be opened"; exit(0);}

    for (int i = 0; i < position.size(); i++)
    {
        f << "v ";
        for (int j = 0; j < position[i].size(); j++)
        {
            f << position[i][j] << " ";
        }
        f << "\n";
    }

    for (auto pair: edges){
        if (pair[0] != -1 && pair[1] != -1)
            f << "f " << pair[0] << " " << pair[1] << "\n";
    }
    
    f.close();
}

#pragma endregion


double sum(std::vector<double>& a){
    double res = 0;
    for (double x: a){
        res += x;
    }
    return res;
}

std::vector<double> mean_sum(std::vector<double>& a, std::vector<double>& b, double weight = 0.5){
    std::vector<double> res(a.size(), 0);
    for (int i = 0; i < a.size(); i++)
    {
        res[i] = weight*a[i] + (1-weight)*b[i];
    }
    return res;    
}

void Smoothing_graph::iteration(double proportion){
    auto inside_res = inside;
    auto outside_res = outside;
    
    for(int v: to_be_changed){
        inside_res[v] = 0;
        outside_res[v] = 0;
        for (auto p : weight[v]){
            if(p.second > 0){
                inside_res[v] += p.second*inside[p.first];
                outside_res[v] += p.second*outside[p.first];
            }
            else{ // If the weight is negative, it means that the neighbour negatively influence the vertex
                inside_res[v] -= p.second*outside[p.first];
                outside_res[v] -= p.second*inside[p.first];
            }
        }
        inside_res[v] /= sum_weight[v];
        outside_res[v] /= sum_weight[v];
    }

    inside = mean_sum(inside, inside_res, proportion);
    outside = mean_sum(outside, outside_res, proportion);
}


void Minimize_area::smooth(){

    std::vector<int> new_classification(classification);

    for (int compteur = 0; compteur < 302; compteur++)
    {
        
        classification = new_classification;
        for (int i = 0; i < n_vertex; i++)
        {
            if (!unchangeable.count(i)){
                double sum = 0;
                for (auto p : weight[i]){
                    sum += p.second * (2*classification[p.first] - 1);
                }
                new_classification[i] = (sum >= 0);
            }
        }
    }
    
    classification = new_classification;

    // do
    // {
    //     classification = new_classification;
    //     for (int i = 0; i < n_vertex; i++)
    //     {
    //         if (!unchangeable.count(i)){
    //             double sum = 0;
    //             for (auto p : weight[i]){
    //                 sum += p.second * (2*classification[p.first] - 1);
    //             }
    //             new_classification[i] = (sum >= 0);
    //         }
    //     }
    //     for(int j: {4176,4210,9966,9991,10802,10803,11109,11919,11920,16656,20574,21702,22433,22434,22435,24726})
    //     std::cerr << classification[j] << " ";
    //     std::cerr <<"\n";
    // } while (new_classification != classification);
    
}