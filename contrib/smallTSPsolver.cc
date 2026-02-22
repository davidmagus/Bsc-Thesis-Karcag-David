#pragma region includok_namespacek
#include <iostream>
#include <fstream>
#include <lemon/maps.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/time_measure.h>
#include <lemon/adaptors.h>
#include "subTask_Bounds.h"
#include <limits>
#include <fstream>
#include <string>
using namespace std;
using namespace lemon;
#pragma endregion

struct task{
    int mask;
    int node;
    double Lbound;

    task( const int _mask, const int _node, const double _Lbound): mask(_mask), node(_node), Lbound(_Lbound){}
    task():mask(0),node(0),Lbound(0.0){}
};


struct Logging {
    std::ofstream logFile;
    int counter;
    Logging() : logFile("debug_log.log") {counter = 1500;}

    template<typename... Args>
    void log(Args... args) {
        if(counter){
            counter--;
            (logFile << ... << args);
        }
    }
};

struct Silent {
    Silent(){}
    static void log(const std::string& msg) {}
};



template <typename DEBUG = Silent, BoundFunc LB = SST>
class STSPsolver

{
    /*This class is used for optimally solving small instances of Assymetric TSP. Use:
        The constructor has the following parameters:
            const ListDigraph& _G,                              The target graph
            const ListDigraph::NodeMap<int>& _Label,            Names of the Nodes (index from 0 to n-1)
            const ListDigraph::ArcMap<double>& _weight,         The weights on the arcs
            const double _Upper_bound)                          An Upper bound on the solution

        Methods:
            .solve(); Returns the value of an optimal TSP solution
    */
private:
#pragma region adattagok
    //Kapott adatok
    int n;
    const ListDigraph& G;
    const ListDigraph::NodeMap<int>& Label;
    vector<ListDigraph::Node> V;
    const ListDigraph::ArcMap<double>& weight;
    double Upper_bound;

    //belso Eszkozok
    double Current_bestval = 0;
    vector<int> Best_route;

    vector<int> Current_route;
    SparseMap<int, double> Lbounds{0.0};
    vector<task>  tasklist; // Listája a A részfeladatoknak amit még ki kell vizsgálni. Ha tasklist[x][0] == mask, tasklist[x][1] == i tasklist[x][2] == k: az x-edik feladat az hogy az i. csúcsból indulva be kell járni azokat a csúcsokat ahol a maskban 1 van úgy hogy idiág k költséggel jutottunk el.
    DEBUG logger;

    #pragma endregion

#pragma region alsokorlat
    void constructBound(const int mask){
        ListDigraph::NodeMap<bool> node_filter(G);

        for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
        {
            node_filter[i] = (mask & (1 << Label[i]));
        }

        FilterNodes<const ListDigraph> subGraph(G, node_filter);
        Lbounds.set(mask, LB(G, weight));
        
    }
#pragma endregion

public:

#pragma region konstruktor
    STSPsolver(
                const ListDigraph& _G,
                const ListDigraph::NodeMap<int>& _Label,
                const ListDigraph::ArcMap<double>& _weight,
                const double _Upper_bound = std::numeric_limits<double>::max()):
                G(_G), weight(_weight),Label(_Label) , Upper_bound(_Upper_bound){
                    n = countNodes(_G);
                    V.resize(n);
                    for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
                    {
                        V[_Label[i]] = i;
                    }
                    logger.log( "Nodes: ", n, " Upperbound: ", Upper_bound, "\n"); 

                };
#pragma endregion


    double solve(){
        //init
        int firstmask = (1 << n) - 2; //összes csúcs kivéve 0
        tasklist.emplace_back( firstmask,0,0.0);
        
        //Futás
        task current_task;
        while (!tasklist.empty())
        {
            current_task = tasklist.back();
            tasklist.pop_back();
            logger.log("Task processed: mask: ", current_task.mask, " Entry Node: ", current_task.node, " Lbound: ", current_task.Lbound, "\n"); 
            for (ListDigraph::OutArcIt a(G, V[current_task.node]); a != INVALID; ++a)
            {
                logger.log("Looping on arcs", "\n");
                int next_mask = current_task.mask & ~(1 << Label[(G).target(a)]);
                if (!next_mask){
                    if ( current_task.Lbound + (weight)[a] < Upper_bound)
                    {
                        Current_bestval = Upper_bound = current_task.Lbound + (weight)[a];
                        Best_route = Current_route;
                        logger.log("foundroute!", "\n");
                    }

                }else if(current_task.Lbound + (weight)[a] + Lbounds[next_mask] <= Upper_bound){
                    if(!Lbounds[next_mask]){
                        constructBound(next_mask);
                        logger.log("Bound built for ", next_mask, " val: ", Lbounds[next_mask],"\n");
                    }
                    if(current_task.Lbound + (weight)[a] + Lbounds[current_task.mask] <= Upper_bound){
                        tasklist.emplace_back(current_task.mask & ~(1 << Label[(G).target(a)]),Label[(G).target(a)],current_task.Lbound + (weight)[a]);
                        logger.log("Task added: mask: ", next_mask, " Entry Node: ", (Label[(G).target(a)]), " Lbound: ", (current_task.Lbound + (weight)[a]), "\n"); 

                    }
                }
            }
            
        }
        
        return Current_bestval;
    }


};


#pragma region main_beolvasas
int main() {
    ListDigraph G;
    ListDigraph::NodeMap<int> Label(G);
    ListDigraph::ArcMap<double> weight(G);
    string filename = "digraph_tsp.lgf";

    // LEMON Timer példányosítása
    // Azonnal elindul a mérés a létrehozás pillanatában
    Timer t;

    cout << "Beolvasas kezdete..." << "\n";

    ifstream f(filename);
    if (!f) {
        cerr << "Hiba: A " << filename << " nem talalhato!" << "\n";
        return 1;
    }

    try {
        DigraphReader<ListDigraph> reader(G, f);
        reader
            .nodeMap("label", Label)
            .arcMap("weight", weight)
            .run();
    } catch (const Exception& e) { // LEMON saját kivétel osztálya
        cerr << "LEMON Hiba: " << e.what() << "\n";
        return 1;
    }

    // Időmérés megállítása (opcionális, a t.userTime() enélkül is az aktuális állást adja)
    t.halt();
#pragma endregion

#pragma region futatas
    std::cout << "Futás kezdete..." << "\n";
    STSPsolver<Logging, SST> ALG(G,Label,weight);
    std::cout << "Opt value: " << ALG.solve() << "\n";
#pragma endregion

}