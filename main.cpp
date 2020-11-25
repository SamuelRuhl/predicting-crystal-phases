#include "main.h"

using namespace std;

int main()
{
    srand(time(NULL));

    int n=40; //dim of the generation
    float a=1.5;
    float b=-1.5;
    int n_generations=100; //anzahl generations


    //create generation with random values.
    vector<Individual> gen;
    for(int i=0;i<n;i++){
        float x = random_num(a,b);
        float y = random_num(a,b);
        gen.push_back(Individual({x,y},a,b));
    }

    ofstream f;
    f.open("/Users/samuelruhl/Documents/Uni/Bachlorarbeit/Code/predicting-crystal-phases/Data/genetic_data.txt");
    f<<"Generation"<<" "<<"x"<<" "<<"y"<<" "<<"f(x,y)"<<" "<<"fitness"<<"\n";

    //perform evolution
    vector<vector<Individual>> Generations = {gen};
    for(int i = 0;i<n_generations;i++){
        vector<Individual> tmp_gen;
        for(int j = 0;j<n/2;j++){
            vector<Individual> kids = selection(Generations[i]);
            tmp_gen.push_back(kids[0]);
            tmp_gen.push_back(kids[1]);
            tmp_gen[j].generation=i;
            tmp_gen[j+1].generation=i;
            f<<i<<" "<<tmp_gen[j].genom[0]<<" "<<tmp_gen[j].genom[1]<<" "<<tmp_gen[j].phaenom<<" "<<tmp_gen[j].fitness<<'\n';
            f<<i<<" "<<tmp_gen[j+1].genom[0]<<" "<<tmp_gen[j+1].genom[1]<<" "<<tmp_gen[j+1].phaenom<<" "<<tmp_gen[j+1].fitness<<'\n';
        }
        Generations.push_back(tmp_gen);
    }
    f.close();
    Generations[n_generations] = sort_by_fitness(Generations[n_generations]);
    Individual winner=Generations[n_generations][n-1];
    cout<<"x="<<winner.genom[0]<<" y="<<winner.genom[0]<<" bei "<<winner.phaenom<<" mit  DNA:";
    winner.print_dna();
    cout<<winner.fitness<<endl;
    //write_data_to_file(Generations);
}
