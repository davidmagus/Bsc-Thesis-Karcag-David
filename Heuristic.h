#ifndef Heuristic_H
#define Heuristic_H
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
using namespace std;
using namespace lemon;

namespace Heuristic{
    struct ClosestNeighbour
    {
        int n;
        const ListDigraph& G;
        const ListDigraph::NodeMap<int>& Label;
        const ListDigraph::ArcMap<double>& weight;
        vector<ListDigraph::Arc> Route;
        double Length;

        #pragma region Constructor, run
        ClosestNeighbour(
                    const ListDigraph& _G,
                    const ListDigraph::NodeMap<int>& _Label,
                    const ListDigraph::ArcMap<double>& _weight
                ): G(_G), weight(_weight), Label(_Label){
                        n = countNodes(_G);
                        Length = 0;
                        ListDigraph::NodeMap<int> seen(G, 0);
                        ListDigraph:: NodeIt v(G);
                        ListDigraph::Node firstnode = v;
                        ListDigraph::Node currentnode = v;
                        seen[currentnode] = 1;
                        while (Route.size() < n-1)
                        {
                            ListDigraph::OutArcIt a(G, currentnode);
                            while(seen[(G).target(a)])
                            {
                                ++a;
                            }
                            ListDigraph::Arc Bestout = a;
                            ++a;
                            while (a !=  INVALID)
                            {
                                if (0 == seen[(G).target(a)] && weight[Bestout] > weight[a])
                                {
                                    Bestout = a;
                                }
                                ++a;    
                            }  
                            Route.push_back(Bestout);
                            currentnode = (G).target(Bestout);
                            Length += weight[Bestout];
                            seen[(G).target(Bestout)] = 1;                   
                        }
                        
                        ListDigraph::Arc lastarc = INVALID;

                        for (ListDigraph::OutArcIt a(G, G.target(Route[n-1])); a != INVALID; ++a) {
                            if (G.target(a) == firstnode) {
                                lastarc = a;
                                break;
                            }
                        }
                        Route.push_back(lastarc);
                        Length += weight[lastarc];
                    };
    #pragma endregion
    };
    
}
#endif