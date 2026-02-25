#ifndef subTask_bounds_H
#define subTask_bounds_H
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#include <lemon/matching.h>
#include <iostream>
using namespace std;
using namespace lemon;

namespace bound{
#pragma region Mincost_tree_bound

    template<typename DG, typename maptype>
    double SSTbound(const DG& G, const maptype& weight, const double Upper_bound){
        std::vector<typename DG::Arc> tree;
        kruskal(G, weight, std::back_inserter(tree));
        double _lbound = 0;
        double minarc = 0;
        for(typename DG::Arc a : tree)
        {
            if(weight[a] < minarc){ minarc = weight[a];}
            _lbound += (weight)[a]; 
        }
        return _lbound;
    };

    struct SST
    {
        template <typename DG, typename maptype>
        static double calculate(const DG& G, const maptype& weight, const double Upper_bound){
            return SSTbound<DG>(G, weight, Upper_bound);
        }
    };
#pragma endregion

#pragma region Assigment_matching_problem_bound

        template<typename DG, typename maptype>
    double APbound(const DG& G, const maptype& weight, const double Upper_bound){
        typename DG::NodeMap<int> G_index(G); 
        
        ListGraph g;
        ListGraph::NodeMap<int> index(g);
        ListGraph::EdgeMap<double> cost(g);
        ListGraph::Node START = g.addNode();
        ListGraph::Node END = g.addNode();

        vector<ListGraph::Node> V;
        V.resize(2 * countNodes(G));


        int i = 0;
        for (typename DG::NodeIt v(G); v != INVALID; ++v)
        {
            G_index[v] = i;
            ListGraph::Node v_source = g.addNode();
            index[v_source] = 2*i;
            V[2*i] = v_source;
            ListGraph::Edge e = g.addEdge(END, v_source);
            cost[e] = 0;

            ListGraph::Node v_sink = g.addNode();
            index[v_sink] = 2*i+1;
            V[2*i+1] = v_sink;
            i++;
            ListGraph::Edge uw = g.addEdge(START, v_sink);
            cost[uw] = 0;
        }
        for (typename DG::ArcIt a(G); a != INVALID; ++a)
        {
            ListGraph::Edge uv = g.addEdge(V[2* G_index[G.source(a)]], V[1 + 2* G_index[G.target(a)]]);
            cost[uv] = -weight[a];
        }

        try{
        MaxWeightedPerfectMatching<ListGraph, ListGraph::EdgeMap<double>> M_ALG(g, cost);
        M_ALG.run();
        return -M_ALG.matchingWeight();
        } catch (const Exception& e) {
            return Upper_bound + 1;
        }

    };

    struct AP
    {
        template <typename DG, typename maptype>
        static double calculate(const DG& G, const maptype& weight, const double Upper_bound){
            return APbound<DG>(G, weight,Upper_bound);
        }
    };
#pragma endregion

#pragma region bruteforce
    struct BF
    {
    };
#pragma endregion
}

#endif
