#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include "Heuristic.h"
using namespace std;
using namespace lemon;

int main()
{
#pragma region main_beolvasas
    ListDigraph G;
    ListDigraph::NodeMap<int> Label(G);
    ListDigraph::ArcMap<double> weight(G);
    string filename = "digraph_tsp.lgf";

    cout << "Beolvasas kezdete..." << "\n";

    ifstream f(filename);
    if (!f)
    {
        cerr << "Hiba: A " << filename << " nem talalhato!" << "\n";
        return 1;
    }

    try
    {
        DigraphReader<ListDigraph> reader(G, f);
        reader
            .nodeMap("label", Label)
            .arcMap("weight", weight)
            .run();
    }
    catch (const Exception &e)
    { // LEMON saját kivétel osztálya
        cerr << "LEMON Hiba: " << e.what() << "\n";
        return 1;
    }
#pragma endregion

    Heuristic::ClosestNeighbour CL{G, Label, weight};
    cout << Label[(G.source(CL.Route[0]))];
    for (size_t i = 0; i < CL.Route.size(); i++)
    {
        cout << Label[(G.target(CL.Route[i]))];
    }
    cout << " Tour Length: " << CL.Length << endl;
}