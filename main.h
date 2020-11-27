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

using namespace std;
Individual selection(vector<Individual> generation);
float random_num( float a, float b);
vector<Individual> sort_by_fitness(vector<Individual> I);
void write_data_to_file(vector<vector<Individual>> I);

#endif // MAIN_H

