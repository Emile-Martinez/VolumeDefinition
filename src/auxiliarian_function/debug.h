#ifndef _DEBUG_ORIENTATION
#define _DEBUG_ORIENTATION

#include <iostream>
#include <vector>
#include <set>
#include <fstream>

class NotImplemented : public std::logic_error
{
public:
    NotImplemented() : std::logic_error("Function not yet implemented") { };
};

void print_vector(std::vector<double> a);

void print_vector(std::vector<int> a);

void print_vector(std::vector<uint32_t> a);

void print_vector(std::vector<uint64_t> a);

void print_vector(std::vector<std::vector<int>> a);

void print_vector(std::vector<std::vector<double>> a);

void print_set(std::set<std::vector<int>> a);

void print_set(std::set<int> a);

void save_vector(const char* fileName, std::vector<std::vector<double>> a);

int number_of_difference(std::vector<int> a, std::vector<int> b);
#endif