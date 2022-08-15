#include "debug.h"

void print_vector(std::vector<double> a){
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << a[i] << "  ";
    }
    std::cout << std::endl;
}

void print_vector(std::vector<int> a){
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << a[i] << "  ";
    }
    std::cout << std::endl;
}


void print_vector(std::vector<uint32_t> a){
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << a[i] << "  ";
    }
    std::cout << std::endl;
}


void print_vector(std::vector<uint64_t> a){
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << a[i] << "  ";
    }
    std::cout << std::endl;
}


void print_vector(std::vector<std::vector<double>> a){
    for (int i = 0; i < a.size(); i++)
    {
        print_vector(a[i]);
    }
    
}

void print_set(std::set<std::vector<int>> a){
    for(auto element : a){
        print_vector(element);
    }
}


void save_vector(const char* fileName, std::vector<std::vector<double>> a){
    std::ofstream f(fileName);

    if (!f) { std::cout<< "Can't save the vector because the file cannot be opened"; exit(0);}

    for (int i = 0; i < a.size(); i++)
    {
        f << "v ";
        for (int j = 0; j < a[i].size(); j++)
        {
            f << a[i][j] << " ";
        }
        f << "\n";
    }

    for (int i = 2; i < a.size(); i++)
    {
        f << "f 0 1 " << i << "\n";
    }
    

    f.close();

}

void print_set(std::set<int> a){
    for (int elem :a) std::cout << elem << " ";
    std::cout << "\n";
}

void print_vector(std::vector<std::vector<int>> a){
    for (int i = 0; i < a.size(); i++)
    {
        print_vector(a[i]);
    }
}


int number_of_difference(std::vector<int> a, std::vector<int> b){
    int res = 0;
    for (int i = 0; i < a.size(); i++)
    {
        if(a[i] != b[i]) std::cerr << i << " ";
        res += (a[i] != b[i]);
    }
    return res;
}