#pragma region includok_namespacek
#include <iostream>
#include <fstream>
#include <lemon/maps.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/time_measure.h>
#include <lemon/adaptors.h>
#include <lemon/kruskal.h>
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


void logDebug(const std::string& message) {
    // std::ios::app = append mód (a fájl végéhez fűzi az új sorokat)
    std::ofstream logFile("debug_log.txt", std::ios::app);
    if (logFile.is_open()) {
        logFile << message << std::endl;
        logFile.close();
    }
}




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
    const ListDigraph*const G;
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
    std::ofstream debugFile;

    #pragma endregion

#pragma region alsokorlat
    void constructBound(const int mask){
        ListDigraph::NodeMap<bool> node_filter(*G);

        for (ListDigraph::NodeIt i(*G); i != INVALID; ++i)
        {
            node_filter[i] = (mask & (1 << Label[i]));
        }

        FilterNodes<const ListDigraph> subGraph(*G, node_filter);
        std::vector<ListDigraph::Arc> tree;
        kruskal(subGraph, weight, std::back_inserter(tree));
        double _lbound = 0;
        for (ListDigraph::Arc a : tree)
        {
            _lbound += (weight)[a]; 
        }
        Lbounds.set(mask, _lbound);
        
    }
#pragma endregion

public:

#pragma region konstruktor
    STSPsolver(
                const ListDigraph& _G,
                const ListDigraph::NodeMap<int>& _Label,
                const ListDigraph::ArcMap<double>& _weight,
                const double _Upper_bound = std::numeric_limits<double>::max()):
                G(&_G), weight(_weight),Label(_Label) , Upper_bound(_Upper_bound){
                    debugFile.open("debug.log");
                    n = countNodes(_G);
                    V.resize(n);
                    for (ListDigraph::NodeIt i(*G); i != INVALID; ++i)
                    {
                        V[_Label[i]] = i;
                    }
                    
                    debugFile << "Nodes: " <<  n << " Upperbound: " << Upper_bound << endl; 

                };
#pragma endregion


    double solve(){
        //init
        tasklist.emplace_back((1 << n) - 1,0,0.0);
        
        //Futás
        task current_task;
        while (!tasklist.empty())
        {
            current_task = tasklist.back();
            tasklist.pop_back();
            debugFile << "Task processed: mask: " << current_task.mask << " Entry Node: " << current_task.node << " Lbound: " << current_task.Lbound << endl; 
            for (ListDigraph::OutArcIt a(*G, V[current_task.node]); a != INVALID; ++a)
            {
                debugFile << "Looping on arcs" << endl;
                if (!(current_task.mask & ~(1 << Label[(*G).target(a)]))){
                        Current_bestval = Upper_bound = current_task.Lbound + (weight)[a];
                        Best_route = Current_route;
                        debugFile << "foundroute!" << endl;
                }else if( (current_task.mask & (1 << Label[(*G).target(a)])) != 0  && current_task.Lbound + (weight)[a] + Lbounds[current_task.mask] <= Upper_bound){
                    if(!Lbounds[current_task.mask]){
                        constructBound(current_task.mask);
                        debugFile << "Bound built for" << current_task.mask << "val: " << Lbounds[current_task.mask] << endl;
                    }
                    if(current_task.Lbound + (weight)[a] + Lbounds[current_task.mask] <= Upper_bound){
                        tasklist.emplace_back(current_task.mask & ~(1 << Label[(*G).target(a)]),Label[(*G).target(a)],current_task.Lbound + (weight)[a]);
                        debugFile << "Task added: mask: " << (current_task.mask & ~(1 << Label[(*G).target(a)])) << " Entry Node: " << (Label[(*G).target(a)]) << " Lbound: " << (current_task.Lbound + (weight)[a]) << endl; 

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

    cout << "Beolvasas kezdete..." << endl;

    ifstream f(filename);
    if (!f) {
        cerr << "Hiba: A " << filename << " nem talalhato!" << endl;
        return 1;
    }

    try {
        DigraphReader<ListDigraph> reader(G, f);
        reader
            .nodeMap("label", Label)
            .arcMap("weight", weight)
            .run();
    } catch (const Exception& e) { // LEMON saját kivétel osztálya
        cerr << "LEMON Hiba: " << e.what() << endl;
        return 1;
    }

    // Időmérés megállítása (opcionális, a t.userTime() enélkül is az aktuális állást adja)
    t.halt();
#pragma endregion

#pragma region futatas
    STSPsolver ALG(G,Label,weight);
    std::cout << "Opt value: " << ALG.solve() << endl;
#pragma endregion

}