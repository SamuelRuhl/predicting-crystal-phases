#include "lattice.h"
#include "main.h"

lattice::lattice(vector<float> z, double density, long double fcc):genom(z){
    //it is importaned to set the variables in this order !
    ar=density;
    set_dna();
    set_primitiv_lattice();
    set_surface();
    //minimize_surface();
    set_lattice_sum();
    set_fitness();
    //set fcc lattice sum for the given density
    fcc_lattice_sum = fcc;
}

void lattice::set_dna(){
    para_dna={};
    for(int i=0; i < int(genom.size()); i++){
        if(i<2){
            //equation (5)  from paper Gottwald, Kahl, Likos
            int b_x = genom[i]*(b_max_12 + 1) - 1; //use equation (5)
            para_dna.push_back(to_binary(b_x, dna_length_x));
        }
        if(i == 2){
            //equation (6)  from paper Gottwald, Kahl, Likos
            int b_phi = 2/M_PI * genom[i] * (b_max_6 + 1) - 1;
            para_dna.push_back(to_binary(b_phi, dna_length_phi));
        }
        if(i == 3){
            int b_psi = 1/M_PI * genom[i] * (b_max_6 + 1) - 1;
            para_dna.push_back(to_binary(b_psi, dna_length_phi));
        }
        if(i == 4){
            int b_psi = 2/M_PI * genom[i] * (b_max_6 + 1) - 1;
            para_dna.push_back(to_binary(b_psi, dna_length_phi));
        }
    }
    long_dna={};
    for(int i=0; i < int(para_dna.size()); i++){
        for(int j=0; j < int(para_dna[i].size()); j++){
            long_dna.push_back(para_dna[i][j]);
        }
    }
}

void lattice::print_dna(){
    for(int i=0;i<int(para_dna.size());i++){
        for(int j=0;j<int(para_dna[i].size());j++){
            cout<<para_dna[i][j];
        }
        cout<<endl;
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
       {genom[0] * sin(genom[2]), genom[0] * cos(genom[2]), 0}, // x2 = (xsin(theta),xcos(theta),0)
       {genom[0] * genom[1] * sin(genom[3]) * sin(genom[4]),    // x3 = (xysin(psi)sin(phi),
        genom[0] * genom[1] * cos(genom[3]) * sin(genom[4]),    //       xycos(psi)sin(phi),
        genom[0] * genom[1] * cos(genom[4])}};                   //       xycos(phi))

}

void lattice::set_surface(){
    surface=calc_surface(x);
}

void lattice::minimize_surface(){               //eq.(7) from ref.(1)
    vector<vector<float>> best_x = x;
    vector<vector<float>> tmp_x;
    vector<vector<float>> last_best_x;
    bool check = true;
    //  could be more sophisticated
    // try all 12 cases
    char c[] = {'+','-'}; //should also run over {'-'}
    while(check){
        last_best_x = best_x;
        for(int i = 0; i < 2; i++){
            //1.
            tmp_x = {sum_vectors(x[0],x[1],c[i]),x[1],x[2]};
            if(calc_surface(tmp_x) < calc_surface(best_x)){
               best_x = tmp_x;
            }
            //2.
            tmp_x = {x[0],sum_vectors(x[1],x[0],c[i]),x[2]};
            if(calc_surface(tmp_x) < calc_surface(best_x)){
                best_x = tmp_x;
            }
            //3.
            tmp_x = {x[0],x[1],sum_vectors(x[2],x[0],c[i])};
            if(calc_surface(tmp_x) < calc_surface(best_x)){
                best_x = tmp_x;
            }
            //4.
            tmp_x = {sum_vectors(x[0],x[2],c[i]),x[1],x[2]};
            if(calc_surface(tmp_x) < calc_surface(best_x)){
               best_x = tmp_x;
            }
            //5.
            tmp_x = {x[0],sum_vectors(x[1],x[2],c[i]),x[2]};
            if(calc_surface(tmp_x) < calc_surface(best_x)){
                best_x = tmp_x;
            }
            //6.
            tmp_x = {x[0],x[1],sum_vectors(x[2],x[1],c[i])};
            if(calc_surface(tmp_x) < calc_surface(best_x)){
                best_x = tmp_x;
            }
        }
        if(calc_surface(last_best_x)<=calc_surface(best_x)){
            check=false;
        }
    }
    if(surface > calc_surface(best_x)){
        x=last_best_x;
        //calc_genom_from_x();
        //set_dna();
        set_surface();
        set_lattice_sum();
        set_fitness();
        //cout<<"find better x"<<endl;
        //have to add something
    }
}

void lattice::set_lattice_sum(){    // eq.(11) from ref.(1)
    long double L = 0;
    int interaction_range = 5; //could also define as constant //8 gives the same value as 200
    for(int l = 1; l < interaction_range; l++){
        for(int k = 0; k < interaction_range; k++){
            for(int m = 0; m < interaction_range; m++){
                L = L + exp((-1/(ar*ar))*(l * l * skalar_product(x[0]) * skalar_product(x[0])
                        + k * k * skalar_product(x[1]) * skalar_product(x[1])
                        + m * m * skalar_product(x[2]) * skalar_product(x[2])));
            }
        }
    }
    for(int k=1; k < interaction_range; k++ ){
        for(int m=0; m < interaction_range; m++){
            L = L + exp((-1/(ar*ar))*( k *k * skalar_product(x[1]) * skalar_product(x[1])   //(0,k,m)
                   + m * m * skalar_product(x[2]) * skalar_product(x[2])));
        }
    }
    for(int m=1; m < interaction_range; m++){
        L = L + exp((-1/(ar*ar)) * m * m * skalar_product(x[2]) * skalar_product(x[2])); //(0,0,m)
    }
    lattice_sum = 2*L;
}

void lattice::set_fitness(){
    //use f=(exp(1-L/L(fcc))
    fitness = exp(1-(lattice_sum/fcc_lattice_sum));
}



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
        vector<char> tmp;
        for(int j = 0; j < int(para_dna[i].size());j++){
            tmp.push_back(tmp_dna[kid_number][j+jump]);
        }
        //using here eq. (5) and (6) from Gottwald, Kahl and Likos
        if(i<2){
            genom_kid.push_back((to_int(tmp) + 1)/(b_max_12 + 1));
            jump=jump + dna_length_x;
        }
        if(i==2){
            genom_kid.push_back(M_PI/2  * (to_int(tmp) + 1)/(b_max_6 + 1));
            jump=jump + dna_length_phi;
        }
        if(i==3){
            genom_kid.push_back(M_PI  * (to_int(tmp) + 1)/(b_max_6 + 1));
            jump=jump + dna_length_phi;
        }
        if(i==4){
            genom_kid.push_back(M_PI/2  * (to_int(tmp) + 1)/(b_max_6 + 1));
            jump=jump + dna_length_phi;
        }
    }

    lattice kid = lattice(genom_kid,ar,fcc_lattice_sum);

    return kid;

}

bool lattice::mutation(float probability){
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    mt19937 mt(seed);
    int p=100/probability;
    uniform_int_distribution<int> prob(0,p);
    return prob(mt);
}

void lattice::print_x(){
    cout<<"x1 =("<<x[0][0]<<","<<x[0][1]<<","<<x[0][2]<<")"<<endl;
    cout<<"x2 =("<<x[1][0]<<","<<x[1][1]<<","<<x[1][2]<<")"<<endl;
    cout<<"x3 =("<<x[2][0]<<","<<x[2][1]<<","<<x[2][2]<<")"<<endl;
}

void lattice::print_genom(){
    cout<<"x="<<genom[0]<<" "<<"y="<<genom[1]<<" "<<"phi="<<genom[2]<<" "<<"psi="<<genom[3]<<" "<<"theta="<<genom[4]<<endl;
}
