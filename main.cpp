#include "main.h"

using namespace std;


#if 1
// ------------------ predicting crystal structures --------------------
int main(){
    //check if fcc is the fittest structure
    lattice test=lattice({0.3,0.7,1.1,1.9,0.5},1);
    lattice fcc=lattice({1,1,M_PI/4,0.0001,M_PI/4},1);
    fcc.print_x();
    cout<<fcc.fitness<<endl;

    cout<<"1)"<<endl;

    test.print_genom();
    //test.print_x();
    test.print_dna();

    cout<<"2)"<<endl;

    test=test.pairing(test);
    test.print_genom();
    //test.print_x();
    test.print_dna();

    cout<<"3)"<<endl;
    test.print_x();

    cout<<"4) Fläsche:"<<endl;
    cout<<test.surface<<endl;

    cout<<"5) Minimierung:"<<endl;
    test.minimize_surface();
    test.print_x();
    cout<<test.surface<<endl;

//    float rho = 1;
//#else

    //run for rho [1,10]
    for(float rho = 1; rho <= 20; rho=rho+5){
    //generate random Population
    ofstream f;
    string s = to_string(rho);
    f.open("/Users/samuelruhl/predicting_crystal_structures/Data/crystal_data_"+s+".txt");
    f<<"mutation probability per bit="<<mutation_prob<<endl;
    f<<"bias fitness="<<lowest_fitness<<endl;
    f<<"number of Individuals per Generation="<<n<<endl;
    f<<"number of Generations="<<n_generations<<endl;

    f<<"Generation"<<" "<<"x"<<" "<<"y"<<" "<<"theta"<<" "<<"psi"<<" "
     <<"phi"<<" "<<"lattice_sum"<<" "<<"fitness"<<"\n";



    //init random number generator
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    //random number between 0 < x <= 1
    auto x_y_rand = std::bind(std::uniform_real_distribution<float>(0.001,1),
                               mt19937(seed));
    //random number between 0 < phi <= pi/2
    auto theta_phi_rand = std::bind(std::uniform_real_distribution<float>(a_phi,b_phi),
                                   mt19937(seed));
    //random number between 0 < theta <= pi
    auto psi_rand = std::bind(std::uniform_real_distribution<float>(a_psi,b_psi),
                                   mt19937(seed));

    //set up random generation
    vector<lattice> gen = {};
    for(int i = 0; i < n; i++){
        vector<float> random_para = {x_y_rand(), x_y_rand(),
                                     theta_phi_rand(), psi_rand(),
                                     theta_phi_rand()};
        gen.push_back(lattice(random_para,rho));
        f<<0<<" "<<gen[i].genom[0]<<" "<<gen[i].genom[1]<<" "<<gen[i].genom[2]<<
           " "<<gen[i].genom[3]<<" "<<gen[i].genom[4]<<" "<<gen[i].lattice_sum<<
           " "<<gen[i].fitness<<"\n";
    }

    //perform evolution
    vector<vector<lattice>> Generations = {};
    Generations.push_back(gen);
    for(int i = 0; i < n_generations; i++){
        vector<lattice> tmp_gen = {};
        //elithism: store the fittest in the new generation
        Generations[i] = sort_lattice_by_fitness(Generations[i]);
        tmp_gen.push_back(Generations[i][n-1]);
        for(int j = 0; j < n - 1; j++){
            lattice kid = lattice_selection(Generations[i]);
            tmp_gen.push_back(kid);
            tmp_gen[j].generation = i + 1;
            f<<i+1<<" "<<tmp_gen[j].genom[0]<<" "<<tmp_gen[j].genom[1]<<
               " "<<tmp_gen[j].genom[2]<<" "<<tmp_gen[j].genom[3]<<
               " "<<tmp_gen[j].genom[4]<<" "<<tmp_gen[j].lattice_sum<<
               " "<<tmp_gen[j].fitness<<"\n";
        }
        Generations.push_back(tmp_gen);
    }
    f.close();
    Generations[n_generations] = sort_lattice_by_fitness(Generations[n_generations]);
    lattice winner = Generations[n_generations][n - 1];
    cout<<"x="<<winner.genom[0]<<" y="<<winner.genom[1]<<
          " theta="<<winner.genom[2]<<" psi="<<winner.genom[3]<<
          " phi="<<winner.genom[4]<<
          " mit  DNA:"<<endl;
    winner.print_long_dna();
    winner.set_primitiv_lattice();
    cout<<"x1 =("<<winner.x[0][0]<<","<<winner.x[0][1]<<","<<winner.x[0][2]<<")"<<endl;
    cout<<"x2 =("<<winner.x[1][0]<<","<<winner.x[1][1]<<","<<winner.x[1][2]<<")"<<endl;
    cout<<"x3 =("<<winner.x[2][0]<<","<<winner.x[2][1]<<","<<winner.x[2][2]<<")"<<endl;
    cout<<"fitness:"<<winner.fitness<<endl;
    cout<<"lattice sum:"<<winner.lattice_sum<<endl;

}

#else
// ----------- find zero-points of the Rosenbrockfunktion ----------------
int main()
{
    ofstream f;
    f.open("/Users/samuelruhl/predicting_crystal_structures/Data/genetic_data.txt");
    f<<"genom dna lenght="<<dna_length<<endl;
    f<<"mutation probability per bit="<<mutation_prob<<endl;
    f<<"bias fitness="<<lowest_fitness<<endl;
    f<<"number of Individuals per Generation="<<n<<endl;
    f<<"number of Generations="<<n_generations<<endl;
    f<<"Search space [a,b]=["<<a<<','<<b<<']'<<endl;

    f<<"Generation"<<" "<<"x"<<" "<<"y"<<" "<<"f(x,y)"<<" "<<"fitness"<<"\n";

    //init random number generator
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    auto real_rand = std::bind(std::uniform_real_distribution<float>(a,b),
                               mt19937(seed));

    //set up random generation
    vector<Individual> gen;
    for(int i=0;i<n;i++){
        float x = real_rand();
        float y = real_rand();
        gen.push_back(Individual({x,y},a,b));
        f<<0<<" "<<gen[i].genom[0]<<" "<<gen[i].genom[1]<<" "<<gen[i].phaenom<<" "<<gen[i].fitness<<'\n';
    }
    //perform evolution
    vector<vector<Individual>> Generations = {gen};
    for(int i = 0;i<n_generations;i++){
        vector<Individual> tmp_gen={};
        //elithismus: take the fittest into the new generation
        Generations[i] = sort_by_fitness(Generations[i]);
        tmp_gen.push_back(Generations[i][n-1]);
        for(int j = 0;j<n-1;j++){
            Individual kid = selection(Generations[i]);
            tmp_gen.push_back(kid);
            tmp_gen[j].generation=i+1;
            f<<i+1<<" "<<tmp_gen[j].genom[0]<<" "<<tmp_gen[j].genom[1]<<" "<<tmp_gen[j].phaenom<<" "<<tmp_gen[j].fitness<<'\n';
        }
        Generations.push_back(tmp_gen);
    }
    f.close();
    Generations[n_generations] = sort_by_fitness(Generations[n_generations]);
    Individual winner=Generations[n_generations][n-1];
    cout<<"x="<<winner.genom[0]<<" y="<<winner.genom[1]<<" bei "<<winner.phaenom<<" mit  DNA:";
    winner.print_dna();
    cout<<winner.fitness<<endl;
#endif
}

