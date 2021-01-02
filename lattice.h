#ifndef LATTICE_H
#define LATTICE_H
#define _USE_MATH_DEFINES
#include <cmath>
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

const int dna_length_x = 12; //strain lenght for x,y
const int dna_length_phi = 6; //strain lenght for phi,psi and theta


//general routine settings
const float mutation_prob = 0.01;
const double lowest_fitness = 10e-100;

const int n=1000; //Number of individuals per generation
const int n_generations=100; //Number of generations

// Interval boarders for x and y
const float a_x=1;
const float b_x=0;

// Interval boarders for theta and phi
const float a_phi=M_PI/2;
const float b_phi=0;

// Interval boarders for psi
const float a_psi=M_PI;
const float b_psi=0;



class lattice
{
private:
    float b_max_12=pow(2,dna_length_x);
    float b_max_6=pow(2,dna_length_phi);
    int numb_para=5;
    float ar;
public:
    //init: gets vector of lenght equal to number_para
    lattice(vector<float> z,double density, long double fcc);

    //attributes
    vector<float> genom; //vector of parameter x,y,theta,psi,phi in that order
    vector<vector<float>> x; //primitiv lattice vectors eq. (1) in ref.1
    vector<vector<char>> para_dna = {};
    vector<char> long_dna = {}; //all para DNA's attached to one
    long double lattice_sum;
    double surface;
    long double fitness;
    long double rel_fit; //relative fitness in ref. to the generation
    int generation;
    long double fcc_lattice_sum;

    //methods
    void set_dna();
    void print_dna();
    void print_long_dna();
    void set_primitiv_lattice();
    void set_surface();
    void minimize_surface(); //to ensure uniqueness
    void set_lattice_sum();
    void set_fitness();
    void set_rel_fit(float summit);
    bool mutation();
    lattice pairing(lattice &other);
    bool mutation(float probability);
    void calc_genom_from_x();
    void print_x();
    void print_genom();

    //overload bigger then operator for siort function
    bool operator< (const lattice &other) const{
        return fitness < other.fitness;
    }

};

#endif // LATTICE_H
