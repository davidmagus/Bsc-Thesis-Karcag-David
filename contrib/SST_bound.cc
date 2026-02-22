#include "subTask_Bounds.h"
#include <lemon/list_graph.h>
#include <lemon/kruskal.h>

double SST(const ListDigraph& G, const ListDigraph::ArcMap<double>& weight){
        std::vector<ListDigraph::Arc> tree;
        kruskal(G, weight, std::back_inserter(tree));
        double _lbound = 0;
        double minarc = 0;
        for (ListDigraph::Arc a : tree)
        {
            if(weight[a] < minarc){ minarc = weight[a];}
            _lbound += (weight)[a]; 
        }
        return _lbound + minarc;
}