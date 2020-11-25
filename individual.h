#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include<iostream>
#include<time.h>
#include<fstream>
#include<vector>
#include<math.h>
#include<iterator>
#include<stdlib.h>

using namespace std;

const int dna_length = 12;

class Individual
{
private:
    float bmax=pow(2,dna_length)-1; //maximal binary number

public:
    //init
    Individual(vector<float> x, float a, float b); //x is element of [a,b]

    //attributes
    vector<float> genom;
    float a;
    float b;
    vector<vector<char>> dna = {};
    float phaenom;
    double fitness;
    double rel_fit;
    int generation;

    //methods
    void set_dna();
    void print_dna();
    void set_fitness();
    void set_rel_fit(float summit);
    void mutation(float probability);
    vector<Individual> pairing(Individual &other);

    //overload bigger then operator for siort function
    bool operator< (const Individual &other) const{
        return fitness < other.fitness;
    }
};

#endif // INDIVIDUAL_H
