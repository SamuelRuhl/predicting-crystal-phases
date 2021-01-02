#include "main.h"

using namespace std;


#if 1
// ------------------ predicting crystal structures --------------------
int main(){
lattice fcc=lattice({1,1,M_PI/4,M_PI/2,M_PI/4},0.4,1);
lattice bcc = lattice({1,1,M_PI/2,3*M_PI/10,M_PI/6},0.4,1);
fcc.print_x();
bcc.print_x();
cout<<setprecision(8);
cout<<pow(0.5,3)<<endl<<endl;
cout<<fcc.lattice_sum<<endl;
cout<<bcc.lattice_sum<<endl<<endl;
cout<<fcc.lattice_sum - bcc.lattice_sum<<endl<<endl;







    //run for rho [1,10]
    for(float ar = 500; ar <= 600 ; ar=ar + 20){
    //generate random Population
    lattice fcc=lattice({1,1,M_PI/4,M_PI/2,M_PI/4},ar/1000,1);
    cout<<"fcc lattice sum = "<<fcc.lattice_sum<<endl;
    cout<<"rho = "<< pow(ar/1000,3)<<endl;
    ofstream f;
    string s = to_string(float(ar));
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
    auto x_y_rand = std::bind(std::uniform_real_distribution<float>(0,1),
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
        gen.push_back(lattice(random_para,ar/1000,fcc.lattice_sum));
        f<<0<<" "<<gen[i].genom[0]<<" "<<gen[i].genom[1]<<" "<<gen[i].genom[2]<<
           " "<<gen[i].genom[3]<<" "<<gen[i].genom[4]<<" "<<gen[i].lattice_sum<<
           " "<<gen[i].fitness<<"\n";
    }

    //perform evolution
    vector<vector<lattice>> Generations = {};
    Generations.push_back(gen);
    bool check = true;
    int i=0;
    while(check){
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
        Generations[i] = sort_lattice_by_fitness(Generations[i]);
        if(i >= n_generations){
            check=false;
        }else{
            i++;
        }
    }
    f.close();
    lattice winner = Generations[i][n - 1];
    cout<<endl<<"x="<<winner.genom[0]<<" y="<<winner.genom[1]<<
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
    float angle_x1_x2 = acos((winner.x[0][0] * winner.x[1][0] + winner.x[0][1] * winner.x[1][1]
                              + winner.x[0][2] * winner.x[1][2])/(sqrt(skalar_product(winner.x[0]))
                                * sqrt(skalar_product(winner.x[1]))));
    float angle_x1_x3 = acos((winner.x[0][0] * winner.x[2][0] + winner.x[0][1] * winner.x[2][1]
                              + winner.x[0][2] * winner.x[2][2])/(sqrt(skalar_product(winner.x[0]))
                                * sqrt(skalar_product(winner.x[2]))));
    float angle_x2_x3 = acos((winner.x[2][0] * winner.x[1][0] + winner.x[2][1] * winner.x[1][1]
                              + winner.x[2][2] * winner.x[1][2])/(sqrt(skalar_product(winner.x[2]))
                                * sqrt(skalar_product(winner.x[1]))));
    cout<<"|x1|="<<sqrt(skalar_product(winner.x[0]))<<endl;
    cout<<"|x2|="<<sqrt(skalar_product(winner.x[1]))<<endl;
    cout<<"|x3|="<<sqrt(skalar_product(winner.x[2]))<<endl;
    cout<<"angle between x1 and x2 :"<<(angle_x1_x2/M_PI) * 180<<endl;
    cout<<"angle between x1 and x3 :"<<(angle_x1_x3/M_PI) * 180<<endl;
    cout<<"angle between x2 and x3 :"<<(angle_x2_x3/M_PI) * 180<<endl;
    cout<<"----------ende Generation a="<<ar<<"------------------"<<endl<<endl;
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

