#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include<iostream>
#include<time.h>
#include<fstream>
#include<vector>
#include<math.h>
#include<iterator>
#include<stdlib.h>
#include<random>
#include<chrono>

using namespace std;

const int dna_length = 12;
//const float mutation_prob = 1;
//const double lowest_fitness = 10e-20;
//const int n=20; //dim of the generation
//const int n_generations=100;
const float a=2;
const float b=0;


class Individual
{
private:
    float bmax=pow(2,dna_length)-1; //maximal binary number
    int numb_para=2;

public:
    //init
    Individual(vector<float> x, float a, float b); //x is element of [a,b]

    //attributes
    vector<float> genom;
    float a;
    float b;
    vector<vector<char>> dna = {};
    vector<char> long_dna = {};
    float phaenom;
    double fitness;
    double rel_fit;
    int generation;

    //methods
    void set_dna();
    void print_dna();
    void set_fitness();
    void set_rel_fit(float summit);
    bool mutation(float probability);
    Individual pairing(Individual &other);

    //overload bigger then operator for siort function
    bool operator< (const Individual &other) const{
        return fitness < other.fitness;
    }
};

#endif // INDIVIDUAL_H
