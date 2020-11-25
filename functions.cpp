#include "main.h"

using namespace std;


vector<Individual> sort_by_fitness(vector<Individual> I){
    sort(I.begin(),I.end());
    //reverse(I.begin(),I.end());
    return I;
}

float random_num( float a, float b){
    return a + (rand() / ( float(RAND_MAX) / (b-a) ) );
}

vector<Individual> selection(vector<Individual> generation){
    //calculate summit fitness and find the minima fitness by the way
    double fit_sum = 0;

    for(int i = 0; i<int(generation.size());i++){
        fit_sum=fit_sum+generation[i].fitness;
        }
   //--------------------------wheel of fortune-------------------------
    double r1 = rand() / double(RAND_MAX);
    generation = sort_by_fitness(generation);
    //wenn es zu problemen kommt dann hier.
        //partner1
        double sum=0;
        int iter=-1;
        float r1_fit=fit_sum * r1;
        while(r1_fit>=sum){
            iter++;
            sum=sum + generation[iter].fitness;
        }
        Individual partner1 = generation[iter];
        //partner2
        fit_sum=fit_sum - generation[iter].fitness;
        generation.erase(generation.begin()+iter);
        sum=0;
        iter=-1;
        double r2 = rand() / double(RAND_MAX);
        float r2_fit=fit_sum * r2;
        while(r2_fit>=sum){
            iter++;
            generation[iter].set_rel_fit(fit_sum);
            sum=sum + generation[iter].fitness;
        }
        Individual partner2 = generation[iter];
    return partner1.pairing(partner2);
}

void write_data_to_file(vector<vector<Individual>> I){
    ofstream f;
    f.open("/Users/samuelruhl/Documents/Uni/Bachlorarbeit/Code/predicting-crystal-phases/Data/genetic_data.txt");
    f<<"Generation"<<" "<<"x"<<" "<<"y"<<" "<<"f(x,y)"<<" "<<"fitness"<<"\n";
    for(int i=0;i<int(I.size());i++){
       for(int j=0;j<int(I[i].size());j++){
           I[i][j].generation=i;
           f<<i<<" "<<I[i][j].genom[0]<<" "<<I[i][j].genom[1]<<" "<<I[i][j].phaenom<<" "<<I[i][j].fitness<<'\n';
       }
    }
    f.close();
}
