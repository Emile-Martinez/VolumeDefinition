#ifndef _GRAPH_SHORTEST_PATH
#define _GRAPH_SHORTEST_PATH

class BSPcomplex;

#include <vector>
#include <map>
#include "tools.h"
#include <set>
#include <unordered_set>
#include "debug.h"
#include <fstream>
#include <unordered_map>
#include <utility>


class Shortest_path_graph{

    protected:
        std::map<int, int> old_to_new;
        std::vector<int> new_to_old;
        std::vector<std::vector<double>> position;
        std::set<std::vector<int>> edges;



    public:
        Shortest_path_graph(){}
        Shortest_path_graph(std::vector<std::vector<double>>);
        void add_node(int v, std::vector<double> barycenter);
        void add_edge(int u, int v);
        void print();
        void save(const char* fileName);
        std::map<int, double> djikstra(int source);
        double weight(int a, int b);
};


Shortest_path_graph get_graph(BSPcomplex* cplx, uint32_t skin_colour);

struct pair_hash{
    std::size_t operator() (const std::pair<int,int> &p){
        return p.first*p.first + p.second;
    }
};


class Smoothing_graph{
    int n_vertex;
    std::vector<std::vector<std::pair<int, double>>> weight;
    std::vector<int> to_be_changed;
    std::vector<double> sum_weight;

    public:

        std::vector<double> inside, outside;
    
        Smoothing_graph(std::vector<double> inside_in, std::vector<double> outside_out)
            : inside(inside_in), outside(outside_out), n_vertex(inside.size()), sum_weight(std::vector<double>(inside_in.size(), 0)), weight(std::vector<std::vector<std::pair<int, double>>>(inside_in.size()))
        {
            for (int i = 0; i < inside_in.size(); i++)
            {
                if (1){
                    to_be_changed.push_back(i);
                }
            }
        }

        Smoothing_graph(std::vector<int> inside_in, std::vector<int> outside_out): n_vertex(inside.size()), sum_weight(std::vector<double>(inside_in.size(), 0)), weight(std::vector<std::vector<std::pair<int, double>>>(inside_in.size()))
        {
            outside = {};
            inside = {};

            for (int i = 0; i < inside_in.size(); i++)
            {
                if (1){//inside_in[i] + outside_out[i] < 5){
                    to_be_changed.push_back(i);
                }
                outside.push_back((double) outside_out[i]);
                inside.push_back((double) inside_in[i]);
            }
            std::cout << "Number of not-touched-enough cells : " << to_be_changed.size() << std::endl;
        }

        void add_edge(int u, int v, double weight_edge){
            weight[u].push_back(std::pair<int, double>(v, weight_edge)); 
            weight[v].push_back(std::pair<int, double>(u, weight_edge));
            sum_weight[u] += std::fabs(weight_edge);
            sum_weight[v] += std::fabs(weight_edge);
            }

        void iteration(double proportion = 0.0);

};


class Minimize_area{
    int n_vertex;
    std::vector<std::vector<std::pair<int, double>>> weight;
    std::unordered_set<int> unchangeable;

    public:
    std::vector<int> classification;
    
    Minimize_area(std::vector<int> is_inside): classification(is_inside), n_vertex(is_inside.size()), weight(std::vector<std::vector<std::pair<int, double>>>(is_inside.size())){};

    void add_edge(int u, int v, double weight_edge){
        weight[u].push_back(std::pair<int, double>(v, weight_edge)); 
        weight[v].push_back(std::pair<int, double>(u, weight_edge));
    }

    void set_unchangeable(int u){
        unchangeable.insert(u);
    }

    void smooth();
};

#endif