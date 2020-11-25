#include "individual.h"



//transforms int to a binary number stoared in a charakter vector
vector<char> to_binary(int n){
    vector<char> c;
    for(int i=dna_length-1;i>=0;i--){
        if((n-pow(2,i)) >= 0){
            c.push_back('1');
            n=n-pow(2,i);
        }else{
            c.push_back('0');
        }
    }
    return c;
}

int to_int(vector<char> c){
    int n=0;
    for(int i=dna_length-1;i>=0;i--){
        if(c[i]=='1'){
            n=n+pow(2,dna_length-i-1);
        }
    }
    return n;
}

//Rosenbrock function f(x,y) = (1-x)^2 + 100(y-x^2)^2 only 2d vector possible x[0] = x; x[1]=y
float Rosenbrock(vector<float> x){
    return (1-x[0])*(1-x[0]) + 100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0]);
}




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
        dna.push_back(to_binary(n));
    }
}

void Individual::print_dna(){
    for(int i=0;i<int(dna.size());i++){
        for(int j=0;j<int(dna[i].size());j++){
            cout<<dna[i][j];
        }
        cout<<" int:"<<to_int(dna[i])<<endl;
    }
}

//-----------------fitness function-------------------------
//optimization proplem is set up here
void Individual::set_fitness(){
    phaenom=Rosenbrock(genom);   //optimization for Rosenbrock function (nullstell)
    fitness=exp(-abs(phaenom)); //fitness function exp(-|phanom|)
    if(fitness==0){
        fitness=1e-26;
    }
}

vector<Individual> Individual::pairing(Individual &other){
    //generate random numbers used as cuttign position for the crossover
    vector<float> genom_kid1;
    vector<float> genom_kid2;
    //still generic, works also for more then two variables and for different DNA lenghts
    for(int i=0;i<int(other.dna.size());i++){
        int cut_pos= (rand() % 12 - 1) + 1; //twister (andere zufallstahlengenerator)
        vector<char> tmp_dna1;
        vector<char> tmp_dna2;
        for(int j=0;j<int(other.dna[i].size());j++){
            //generate crossover dna
            if(j<=cut_pos){
                tmp_dna1.push_back(dna[i][j]);          //nur ein Kind
                tmp_dna2.push_back(other.dna[i][j]);
            }
            if(j>cut_pos){
                tmp_dna1.push_back(other.dna[i][j]);
                tmp_dna2.push_back(dna[i][j]);
            }
        }
        //calculate from the dna the genoms to build up the kids
        genom_kid1.push_back(to_int(tmp_dna1) * abs(a-b) / bmax + b); //used transofmation genom=n*|a-b|/bmax + b
        genom_kid2.push_back(to_int(tmp_dna2) * abs(a-b) / bmax + b);
    }
    return {Individual(genom_kid1,a,b),Individual(genom_kid2,a,b)};
}

void Individual::set_rel_fit(float summit){
    rel_fit=fitness/summit;
}
