#include "main.h"

using namespace std;

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


vector<Individual> sort_by_fitness(vector<Individual> I){
    sort(I.begin(),I.end());
    //reverse(I.begin(),I.end());
    return I;
}

Individual selection(vector<Individual> generation){
    //calculate summit fitness and find the minima fitness by the way
    double fit_sum = 0;

    for(int i = 0; i<int(generation.size());i++){
        fit_sum=fit_sum+generation[i].fitness;
        }
   //--------------------------wheel of fortune-------------------------
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    auto real_rand = std::bind(std::uniform_real_distribution<float>(0,1),
                               mt19937(seed));
    double r1 = real_rand();
    generation = sort_by_fitness(generation);
    //wenn es zu problemen kommt dann hier.
        //partner1
        double sum=0;
        int iter=-1;
        float r1_fit=fit_sum * r1;
        while(r1_fit>sum){
            iter++;
            sum=sum + generation[iter].fitness;
        }
        Individual partner1 = generation[iter];
        //partner2
        fit_sum=fit_sum - generation[iter].fitness;
        generation.erase(generation.begin()+iter);
        sum=0;
        iter=-1;
        double r2 = real_rand();
        float r2_fit=fit_sum * r2;
        while(r2_fit>=sum){
            iter++;
            sum=sum + generation[iter].fitness;
        }
        Individual partner2 = generation[iter];
    return partner1.pairing(partner2);
}


//-------------------------functions for lattice class---------------------

vector <lattice> sort_lattice_by_fitness(vector<lattice> I){
    sort(I.begin(),I.end());
    //reverse(I.begin(),I.end());
    return I;
}

lattice lattice_selection(vector<lattice> generation){
    //calculate summit fitness and find the minima fitness by the way
    double fit_sum = 0;

    for(int i = 0; i<int(generation.size());i++){
        fit_sum=fit_sum+generation[i].fitness;
        }
   //--------------------------wheel of fortune-------------------------
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    auto real_rand = std::bind(std::uniform_real_distribution<float>(0,1),
                               mt19937(seed));
    double r1 = real_rand();
    generation = sort_lattice_by_fitness(generation);
    //wenn es zu problemen kommt dann hier.
        //partner1
        double sum=0;
        int iter=-1;
        float r1_fit=fit_sum * r1;
        while(r1_fit>sum){
            iter++;
            sum=sum + generation[iter].fitness;
        }
        lattice partner1 = generation[iter];
        //partner2
        fit_sum=fit_sum - generation[iter].fitness;
        generation.erase(generation.begin()+iter);
        sum=0;
        iter=-1;
        double r2 = real_rand();
        float r2_fit=fit_sum * r2;
        while(r2_fit>=sum){
            iter++;
            sum=sum + generation[iter].fitness;
        }
        lattice partner2 = generation[iter];
    return partner1.pairing(partner2);
}

vector<float> cross_product(vector<float> a, vector<float> b){
    vector<float> c;
    float x1 = a[1] * b[2] - a[2] * b[1];
    float x2 = a[2] * b[0] - a[0] * a[2];
    float x3 = a[0] * b[1] - a[1] * b[0];
    c={x1,x2,x3};
    return c;
}

float vec_product(vector<float> a){
    float abs=0;
    for(int i=0; i < int(a.size()); i++){
        abs = abs + a[i] * a[i];
    }
    return sqrt(abs);
}

float calc_surface(vector<vector<float>> x){
    return  (vec_product(cross_product(x[0],x[1]))
            + vec_product(cross_product(x[0],x[2]))
            + vec_product(cross_product(x[1],x[2])));
}

vector<float> sum_vectors(vector<float> a, vector<float> b, char o){
    vector<float> c;
    for(int i = 0 ; i < int(a.size()); i++){
        if(o=='-'){
            c.push_back(a[i]-b[i]);
        }
        else{
            c.push_back(a[i]+b[i]);
        }
    }
    return c;
}
