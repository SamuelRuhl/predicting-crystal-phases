#include "individual.h"
#include "main.h"


auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mt(seed);

//----------------------class starts here-------------------------------

Individual::Individual(vector<float> x,float a, float b):genom(x),a(a),b(b){
    set_dna();
    set_fitness();
};

//for 2d functions (x,y) take
void Individual::set_dna(){
    for(int i=0;i<int(genom.size());i++){
        float tmp_x=(genom[i]-b);   //shift x by -b
        float resolution=abs(a-b)/bmax;  //used here n=(genom[i]-b)*bmax/|a-b|
        int n=tmp_x/resolution;
        dna.push_back(to_binary(n, dna_length));
    }
    for(int i=0; i < int(dna.size()); i++){
        for(int j=0; j < int(dna[i].size()); j++){
            long_dna.push_back(dna[i][j]);
        }
    }
}

void Individual::print_dna(){
    for(int i=0;i<int(dna.size());i++){
        for(int j=0;j<int(dna[i].size());j++){
            cout<<long_dna[j + i*dna_length];
        }
    }
    cout<<endl;
}

//-----------------fitness function-------------------------
//optimization proplem is set up here
void Individual::set_fitness(){
    phaenom=Rosenbrock(genom);   //optimization for Rosenbrock function (nullstell)
    fitness=exp(-abs(phaenom)); //fitness function exp(-|phanom|)
    if(fitness==0){
        fitness=lowest_fitness;
    }
}

Individual Individual::pairing(Individual &other){
    //generate random numbers used as cuttign position for the crossover
    vector<float> genom_kid;
    //still generic, works also for more then two variables and for different DNA lenghts
    uniform_int_distribution<int> dist(0, numb_para*dna_length);
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
    //calculate from the dna the genoms to build up the kids
    //chose a random kid
    uniform_int_distribution<int> ran(0, 1);
    int random_number = ran(mt);
    for(int i=0;i<int(numb_para);i++){
        vector<char> geno;
        for(int j=0;j<int(dna_length);j++){
            //mutation
            if(mutation(mutation_prob)==0){
                if(tmp_dna[random_number][j+i*dna_length]=='0'){
                    tmp_dna[random_number][j+i*dna_length]='1';
                }else{
                    tmp_dna[random_number][j+i*dna_length]='0';
                }
            }
            geno.push_back(tmp_dna[random_number][j+i*dna_length]);
        }
        genom_kid.push_back(to_int(geno) * abs(a-b) / bmax + b); //used transofmation genom=n*|a-b|/bmax + b
    }
//    this->print_dna();
//    other.print_dna();
    Individual kid = Individual(genom_kid,a,b);
//    cout<<cut_pos<<endl;
//    kid.print_dna();
//    cout<<endl;
    return kid;
}

void Individual::set_rel_fit(float summit){
    rel_fit=fitness/summit;
}

bool Individual::mutation(float probability){
    int p=100/probability;
    uniform_int_distribution<int> prob(0,p);
    return prob(mt);
}
