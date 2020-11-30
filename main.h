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
#include "individual.h"
#include "lattice.h"

using namespace std;
vector<char> to_binary(int n);
int to_int(vector<char> c);
float Rosenbrock(vector<float> x);
Individual selection(vector<Individual> generation);
vector<Individual> sort_by_fitness(vector<Individual> I);
lattice lattice_selection(vector<lattice> generation);
vector<lattice> sort_lattice_by_fitness(vector<lattice> I);

#endif // MAIN_H

