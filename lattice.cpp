#include "lattice.h"
#include "main.h"

lattice::lattice(vector<float> z):genom(z){
//    set_dna();
//    set_lattice_sum();
//    set_fitness();
//    set_prmitiv_lattice();
}

void lattice::set_dna(){
    for(int i=0; i < int(genom.size()); i++){
        if(i<2){
            //equation (5)  from paper Gottwald, Kahl, Likos
            int b_x = genom[i]*(b_max_12 + 1) - 1; //use equation (5)
            para_dna.push_back(to_binary(b_x));
        }
        else{
            //equation (6)  from paper Gottwald, Kahl, Likos
            int b_phi = 2/M_PI * genom[i] * (b_max_6 + 1) - 1;
            para_dna.push_back(to_binary(b_phi));
        }
    }
    for(int i=0; i < int(para_dna.size()); i++){
        for(int j=0; j < int(para_dna[i].size()); j++){
            long_dna.push_back(para_dna[i][j]);
        }
    }
}

void lattice::print_dna(){
    for(int i=0;i<int(para_dna.size());i++){
        for(int j=0;j<int(para_dna[i].size());j++){
            cout<<para_dna[i][j]<<endl;
        }
    }
}

void lattice::print_long_dna(){
    for(int i=0; i < int(long_dna.size());i++){
        cout<<long_dna[i];
    }
    cout<<endl;
}

void lattice::set_primitiv_lattice(){
    x={{1,0,0},                                                 // x1 = (1,0,0)
       {genom[0] * cos(genom[2]), genom[0] * sin(genom[2]), 1}, // x2 = (xcos(theta),xsin(theta),0)
       {genom[0] * genom[1] * cos(genom[3]) * cos(genom[4]),    // x3 = (xycos(psi)sin(phi),
        genom[0] * genom[1] * cos(genom[3]) * cos(genom[4]),    //       xycos(psi)sin(phi),
       genom[0] * genom[1] * sin(genom[4])}};                   //       xysin(phi))
}

//void lattice::minimize_surface(){} //eq.(7) from ref.(1)

//void lattice::set_lattice_sum(){} // eq.(11) from ref.(1)

//void lattice::set_fitness(){}



lattice lattice::pairing(lattice &other){

    //generate random numbers used as cuttign position for the crossover
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 mt(seed);
    uniform_int_distribution<int> dist(0, long_dna.size());
    int cut_pos= dist(mt);
    vector<vector<char>> tmp_dna={{},{}};
    for(int i=0;i<int(other.long_dna.size());i++){
        if(i<cut_pos){
            tmp_dna[0].push_back(long_dna[i]);
            tmp_dna[1].push_back(other.long_dna[i]);
        }
        if(i>=cut_pos){
            tmp_dna[0].push_back(other.long_dna[i]);
            tmp_dna[1].push_back(long_dna[i]);
        }
    }

    //chose a random kid
    uniform_int_distribution<int> ran(0, 1);
    int kid_number = ran(mt);
    //check every bit of the chosen kid for mutation
    for(int i = 0; i < int(long_dna.size());i++){
        if(mutation(mutation_prob)==0){
            if(tmp_dna[kid_number][i]=='0'){
                tmp_dna[kid_number][i]='1';
            }else{
                tmp_dna[kid_number][i]='0';
            }
        }
    }

    //calculate from the dna the genoms to build up the choosen kid
    vector<float> genom_kid;
    int jump=0;
    for(int i = 0; i < int(para_dna.size());i++){
        vector<char> tmp_dna;
        for(int j = 0; j < int(para_dna[i].size());i++){
            tmp_dna.push_back(long_dna[j+jump]);
        }
        //using here eq. (5) and (6) from Gottwald, Kahl and Likos
        if(i<2){
            genom_kid.push_back((to_int(tmp_dna) + 1)/(b_max_12 + 1));
            jump=jump + dna_length_x;
        }else{
            genom_kid.push_back(M_PI/2  * (to_int(tmp_dna) + 1)/(b_max_6 + 1));
            jump=jump + dna_length_phi;
        }
    }

    lattice kid = lattice(genom_kid);

    return kid;

}

bool lattice::mutation(float probability){
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 mt(seed);
    int p=100/probability;
    uniform_int_distribution<int> prob(0,p);
    return prob(mt);
}
