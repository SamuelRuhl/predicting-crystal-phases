#ifndef MAIN_H
#define MAIN_H

#include<iostream>
#include<time.h>
#include<random>
#include<fstream>
#include<vector>
#include<math.h>
#include<iterator>
#include<stdlib.h>
#include<string.h>
#include "individual.h"
#include "lattice.h"

using namespace std;
vector<char> to_binary(int n, int lenght_dna = 12);
int to_int(vector<char> c);
float Rosenbrock(vector<float> x);
Individual selection(vector<Individual> generation);
vector<Individual> sort_by_fitness(vector<Individual> I);
lattice lattice_selection(vector<lattice> generation);
vector<lattice> sort_lattice_by_fitness(vector<lattice> I);
vector<float> cross_product(vector<float> a, vector<float> b);
float skalar_product(vector<float> a);
float calc_surface(vector<vector<float>> x);
vector<float> sum_vectors(vector<float> a, vector<float> b, char o);
void print_vec(vector<float> a);

#endif // MAIN_H

