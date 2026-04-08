#pragma once
#pragma region includes namespaces
#include <iostream>
#include <fstream>
#include <lemon/maps.h>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/time_measure.h>
#include <lemon/adaptors.h>
#include "subTask_Bounds.h"
#include "Heuristic.h"
#include <limits>
#include <fstream>
#include <stdexcept>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
#include <iostream>
#include <lemon/unionfind.h>
#include <type_traits>
using namespace std;
using namespace lemon;
#pragma endregion

namespace Heldkarp
{
#pragma region BnB tree node
    struct tree_node
    {
        bool finalize;
        boost::dynamic_bitset<> mask;
        ListDigraph::Node node;
        ListDigraph::Arc arc;
        double LB;
        tree_node(const bool _finalize,
                  boost::dynamic_bitset<> &_mask,
                  const ListDigraph::Node &_node,
                  const ListDigraph::Arc &_arc) : finalize(_finalize),
                                                  mask(_mask),
                                                  node(_node),
                                                  arc(_arc)
        {
        }

        tree_node(double N) : LB() {}
    };
#pragma endregion

#pragma region Logging tools
    struct Logging
    {
        std::ofstream logFile;
        int counter;
        Logging() : logFile("debug.log") { counter = 1500; }

        template <typename... Args>
        void log(Args... args)
        {
            if (counter)
            {
                counter--;
                (logFile << ... << args);
            }
        }
    };
    struct Silent
    {
        Silent() {}
        template <typename... Args>
        void log(Args... args)
        {
        }
    };

#pragma endregion

    template <typename DEBUG = Silent> // Template parameters: <DEBUG = Heldkarp::Silent / Heldkarp::Logging, LB = bound::BF / bound::AP / bound::SST, alpha_approx = Heldkarp::approx(double Alpha)> parameters: const ListDigraph& _G, const ListDigraph::NodeMap<int>& _Label, const ListDigraph::ArcMap<double>& _weight, const double _Upper_bound
    class Heldkarp
    {
    private:
#pragma region Members
        // Kapott adatok
        // const double approx_ratio;
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        vector<ListDigraph::Node> V;
        const ListDigraph::ArcMap<double> &weight;
        double Upper_bound;

        // Belső változók
        vector<ListDigraph::Arc> Arcs; // Az élek egy rendezett listája

        // Eredmények
        double Current_bestval = 0;
        vector<ListDigraph::Arc> Route;
        vector<ListDigraph::Arc> Best_Route;
        int solved = 0; // A query függvények számára fentartott belső ellenörző 0 ha nem hívták még meg a .solve metódust, 1 egyébként.

        // DP
        vector<SparseMap<boost::dynamic_bitset<>, double>> LB;                // Alsóbecslések szótára
        vector<SparseMap<boost::dynamic_bitset<>, ListDigraph::Arc>> NextArc; // Következő él
        vector<SparseMap<boost::dynamic_bitset<>, double>> DP;                // Konkrét értékek szótára

        // Branch and Bound fa
        ListDigraph r;              // Held Karp féle gyökér
        vector<tree_node> tasklist; // Listája a A részfeladatoknak amit még ki kell vizsgálni.
        DEBUG logger;               // A logoláshoz használt objektum

#pragma endregion

#pragma region Initialization
        void init()
        {
            // Élek rendezése
            for (ListDigraph::ArcIt a(G); a != INVALID; a++)
            {
                Arcs.push_back(a);
            }

            std::sort(Arcs.begin(), Arcs.end(),
                      [&weight](const ListDigraph::Arc &a, const ListDigraph::Arc &b)
                      {
                          return weight[a] < weight[b];
                      });

            // Felső becslés beállítása
            Heuristic::Max_insert maxins{G, Label, weight};
            maxins.Run();
            if (Upper_bound > maxins.Length)
            {
                Upper_bound = maxins.Length;
                Route = maxins.Route;
            }

            Heuristic::Repetitive_Nearest_Neighbour RNN{G, Label, weight};
            if (Upper_bound > RNN.Length)
            {
                Upper_bound = RNN.Length;
                Route = RNN.Route;
            }

            // Tasklist Vector előkészítése
            tasklist.clear();
            boost::dynamic_bitset<> H(n - 1);
            H.set();
            for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
            {
                if (Label[v] = n - 1)
                {
                    r = v;
                }
                else
                {

                    H.reset(Label[v]);
                    tasklist.emplace_back(false, H, v);
                }
            }
        }
#pragma endregion

#pragma region Route found
        void Route_found(tree_node &X, tree_node &START)
        {
            if (G.source(Route.back()) != r)
            {
                Route.push_back(NextArc[X.node][X.mask]);
                boost::dynamic_bitset<> newmask = X.mask;
                newmask.reset(Label[G.source(NextArc[X.node][X.mask])]);
                ListDigraph::Node newnode = G.source(NextArc[X.node][X.mask]);
                tree_node new_X(false, newmask, newnode, NextArc[X.node][X.mask]);

                Route_found(new_X, START);
                return;
            }

            double val = 0.0;
            for (size_t i = 0; i < Route.size(); i++)
            {
                val += weight[Route[i]];
            }

            if (val < Upper_bound)
            {
                Upper_bound = val;
                Best_Route = Route;
            }

            while (START.arc != Route.back())
            {
                Route.pop_back();
            }
        }
#pragma endregion

#pragma region Lower Bound
        struct DSU
        {
            vector<int> parent;
            DSU(int n)
            {
                parent.resize(n);
                for (int i = 0; i < n; ++i)
                    parent[i] = i;
            }

            int find(int i)
            {
                if (parent[i] == i)
                    return i;
                return parent[i] = find(parent[i]); // Útvonal-tömörítés
            }

            void unite(int i, int j)
            {
                int root_i = find(i);
                int root_j = find(j);
                if (root_i != root_j)
                    parent[root_i] = root_j;
            }
        };

        double Bound(tree_node &X)
        {
            DSU dsu{n};
            int edge_count = 0;
            double rin{Upper_bound + 1}, rout{Upper_bound + 1}, tin{Upper_bound + 1}, FH{0}, FC{0};
            for (size_t i = 0; i < Arcs.size(); i++)
            {
                ListDigraph::Arc &a = Arcs[i];

                if (G.source(a) == r)
                {
                    if (X.mask[Label[G.target(a)]] && rout == Upper_bound + 1)
                    {
                        rout = weight[a];
                    }
                }
                else if (G.target(a) == r)
                {
                    if (G.source(a) != X.node && !X.mask[Label[G.source(a)]] && rin == Upper_bound + 1)
                    {
                        rin = weight[a];
                    }
                }
                else if (X.node == G.target(a))
                {
                    if (X.mask[Label[G.source(a)]] && tin == Upper_bound + 1)
                    {
                        tin = weight[a];
                    }
                }
                else if (X.node != G.source(a))
                {
                    if (!X.mask[Label[G.target(a)]] && !X.mask[Label[G.source(a)]])
                    {
                        if (find(Label[G.target(a)]) != find(Label[G.source(a)]))
                        {
                            edge_count++;
                            FC += weight[a];
                            dsu.unite(Label[G.target(a)], Label[G.source(a)]);
                        }
                    }
                    else if (X.mask[Label[G.target(a)]] && X.mask[Label[G.source(a)]])
                    {
                        if (find(Label[G.target(a)]) != find(Label[G.source(a)]))
                        {
                            edge_count++;
                            FH += weight[a];
                            dsu.unite(Label[G.target(a)], Label[G.source(a)]);
                        }
                    }
                }
            }
            if (edge_count < n - 4)
            {
                return Upper_bound + 1;
            }
            return LB[X.node][X.mask] = rin + rout + tin + FH + FC;
        }
#pragma endregion

#pragma region processnode
        void process(tree_node &X, ArcLookUp<ListDigraph> &lookUp)
        {
            if (X.finalize)
            {
                Route.pop_back();
                double minval = Upper_bound + 1;
                ListDigraph::Arc tempNextArc = a;
                ListDigraph::InArcIt a(G, X.node);
                for (ListDigraph::InArcIt a(G, X.node); a != INVALID; ++a)
                {
                    X.mask.flip(Label[G.source(a)]);
                    if (minval > DP[Label[G.source(a)]][X.mask])
                    {
                        minval = DP[Label[G.source(a)]][X.mask];
                        tempNextArc = NextArc[Label[G.source(a)]][X.mask];
                    }
                    X.mask.flip(Label[G.source(a)]);
                }
                DP[Label[X.node]][X.mask] = LB[Label[X.node]][X.mask] = minval;
                NextArc[Label[X.node]][X.mask] = tempNextArc;
            }
            else
            {
                Route.push_back(X.arc);
                if (!X.mask.any())
                {
                    if (lookUp(r, X.node) == INVALID)
                    {
                        DP[Label[X.node]][X.mask] = LB[Label[X.node]][X.mask] = Upper_bound + 1;
                        return;
                    }

                    Route.push_back(lookUp(r, X.node));
                    DP[Label[X.node]][X.mask] = LB[Label[X.node]][X.mask] = weight[lookUp(r, X.node)];
                    NextArc[Label[X.node]][X.mask] = lookUp(r, X.node);
                    Route_found(X, X);
                    return;
                }
                vector<tree_node> tasks;
                for (ListDigraph::InArcIt a(G, X.node); a != INVALID; ++a)
                {
                    boost::dynamic_bitset<> newmask = X.mask;
                    newmask.reset(Label[G.source(a)]);
                    ListDigraph::Node newnode = G.source(a);
                    tree_node new_X(false, newmask, newnode, a);
                    if (LB[newnode][newmask] == 0)
                    {
                        new_X.LB = bound(new_X) + weight[a];
                        if (new_X.LB < Upper_bound)
                        {
                            tasks.push_back(new_X);
                        }
                    }
                    else
                    {
                        Route_found(X, X);
                    }
                }
                std::sort(tasks.begin(), tasks.end(), [](const tree_node &a, const tree_node &b)
                          {
                              return a.LB > b.LB; // Csökkenő sorrend
                          });
                X.finalize = true;
                tasklist.push_back(X);
                for (size_t i = 0; i < tasks.size(); i++)
                {
                    if (tasks[i].LB < Upper_bound)
                    {
                        tasklist.push_back(tasks[i]);
                    }
                }
            }
        }

#pragma endregion

    public:
#pragma region Constructor
        Heldkarp(
            const ListDigraph &_G,
            const ListDigraph::NodeMap<int> &_Label,
            const ListDigraph::ArcMap<double> &_weight,
            const double _Upper_bound = std::numeric_limits<double>::max()) : G(_G), weight(_weight), Label(_Label), Upper_bound(_Upper_bound)
        {
            n = countNodes(_G);
            V.resize(n);
            for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
            {
                V[_Label[i]] = i;
            }
            logger.log("Nodes: ", n, " Upperbound: ", Upper_bound, "\n");

            // 1. LB inicializálása, ha még nem számoltunk ki jobbat legyen 0
            LB.assign(n - 1, SparseMap<boost::dynamic_bitset<>, double>(0.0));

            // 2. NextArc inicializálása INVALID alapból
            NextArc.assign(n - 1, SparseMap<boost::dynamic_bitset<>, ListDigraph::Arc>(INVALID));

            // 3. Val inicializálása legyen alapból egy nagy szám
            DP.assign(n - 1, SparseMap<boost::dynamic_bitset<>, double>(Upper_bound + 1));
        };
#pragma endregion

#pragma region Query functions
        double solve()
        {
            lemon::ArcLookUp<ListDigraph> lookUp(G);
            init();

            solved = 1;
            return Current_bestval;
        }

        double OPTval()
        {
            if (!solved)
            {
                throw std::invalid_argument("Most call method .solve() before using query functions");
            }
            return Current_bestval;
        }

        vector<ListDigraph::Arc> OPTroute()
        {
            if (!solved)
            {
                throw std::invalid_argument("Most call method .solve() before using query functions");
            }
            return Best_;
        }

        void printroute()
        {
            for (int i = 1; i < Best_route.size(); i++)
            {
                cout << " " << Best_route[i];
            }
            cout << endl;
        }
#pragma endregion
    };
}
