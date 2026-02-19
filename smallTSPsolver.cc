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

using namespace std;
using namespace lemon;
#pragma endregion

struct task{
    int mask;
    int node;
    double Lbound;
};


class STSPsolver
{
private:
#pragma region adattagok
    //Kapott adatok
    const int n;
    const ListDigraph*const G;
    const ListDigraph::NodeMap<int>& Label;
    vector<ListDigraph::Node> V;
    const ListDigraph::ArcMap<double>* weight;
    double Upper_bound;

    //belso Eszkozok
    double current_best = 0;
    vector<int> Best_route;

    vector<int> Current_route;
    SparseMap<int, double> Lbounds{0.0};
    vector<task>  tasklist; // Listája a A részfeladatoknak amit még ki kell vizsgálni. Ha tasklist[x][0] == mask, tasklist[x][1] == i tasklist[x][2] == k: az x-edik feladat az hogy az i. csúcsból indulva be kell járni azokat a csúcsokat ahol a maskban 1 van úgy hogy idiág k költséggel jutottunk el.
#pragma endregion

#pragma region alsokorlat
    double constructBound(const int mask){
        ListDigraph::NodeMap<bool> node_filter(*G);

        for (ListDigraph::NodeIt i(*G); i != INVALID; ++i)
        {
            node_filter[i] = (mask & (1 << Label[i]));
        }

        FilterNodes<const ListDigraph> subGraph(*G, node_filter);
        std::vector<ListDigraph::Arc> tree;
        kruskal(subGraph, *weight, std::back_inserter(tree));
        double _lbound = 0;
        for (ListDigraph::Arc a : tree)
        {
            _lbound += (*weight)[a]; 
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
                const double _Upper_bound): G(&_G), weight(&_weight),Label(_Label) , Upper_bound(_Upper_bound), n(countNodes(G)){
                    V.resize(n);
                    for (ListDigraph::NodeIt i(*G); i != INVALID; ++i)
                    {
                        V[_Label[i]] = i;
                    }
                    
                };
#pragma endregion


    double solve(){
        //init
        tasklist.push_back(task({(1 << n) - 1,0,0}));
        
        //Futás
        task current_task;
        while (tasklist.empty())
        {
            current_task = tasklist.back();
            tasklist.pop_back();
            for (ListDigraph::OutArcIt a(*G, V[current_task.node]); a != INVALID; ++a)
            {
                if (!(current_task.mask & ~(1 << Label[(*G).target(a)]))){
                        current_best = Upper_bound = current_task.Lbound + (*weight)[a];
                        Best_route = Current_route;
                }else if( (current_task.mask & (1 << Label[(*G).target(a)])) != 0  && current_task.Lbound + (*weight)[a] + Lbounds[current_task.mask] <= Upper_bound){
                    if(!Lbounds[current_task.mask]){
                        constructBound(current_task.mask);
                    }
                    if(current_task.Lbound + (*weight)[a] + Lbounds[current_task.mask] <= Upper_bound){
                        int next_mask = current_task.mask & ~(1 << Label[(*G).target(a)]);
                        double next_bound = current_task.Lbound + (*weight)[a];
                        int next_node = Label[(*G).target(a)];
                        tasklist.emplace_back(next_mask,next_node,next_bound); // Ez vmi olyat csinál hogy egyből a tömb végén építí fel az adatszerekezetet
                    }
                }
            }
            
        }
        
    }


};


#pragma region main_beolvasas
int main() {
    ListDigraph G;
    ListDigraph::NodeMap<int> Label(G);
    ListDigraph::ArcMap<int> weight(G);
    string filename = "digraph.lgf";

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
}