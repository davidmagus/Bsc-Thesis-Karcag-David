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
        uint64_t mask;
        ListDigraph::Node node;
        ListDigraph::Arc arc;
        double LB;
        tree_node(const bool _finalize,
                  uint64_t &_mask,
                  const ListDigraph::Node &_node,
                  const ListDigraph::Arc &_arc) : finalize(_finalize),
                                                  mask(_mask),
                                                  node(_node),
                                                  arc(_arc)
        {
        }
    };
#pragma endregion

#pragma region Logging tools
    struct Logging
    {
        std::ofstream logFile;
        int counter;
        Logging() : logFile("HK.log") { counter = 15000; }

        template <typename... Args>
        void log(Args... args)
        {
            if (counter)
            {
                counter--;
                (logFile << ... << args) << endl;
            }
        }
    };
struct Silent
{
    Silent() {}
    template <typename... Args>
    void log(Args...)
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
        double Current_bestval = 1000000000;
        vector<ListDigraph::Arc> Route;
        double R_val = 0;
        void add(const ListDigraph::Arc &arc)
        {
            Route.push_back(arc);
            R_val += weight[arc];
            logger.log("                               Added: ", Label[G.source(arc)], "->", Label[G.target(arc)]);
        }
        void pop()
        {
            if (!Route.empty())
            {
                ListDigraph::Arc a = Route.back();
                logger.log("                               Removed: ", Label[G.source(a)], "->", Label[G.target(a)]);
                R_val -= weight[a];
                Route.pop_back();
            }
            else
            {
                throw invalid_argument(" Popped when empty");
            }
        }
        void R_update()
        {
            double newval = 0.0;
            for (ListDigraph::Arc a : Route)
            {
                newval += weight[a];
            }
            R_val = newval;
        }

        vector<ListDigraph::Arc> Best_Route;
        int solved = 0; // A query függvények számára fentartott belső ellenörző 0 ha nem hívták még meg a .solve metódust, 1 egyébként.

        // DP
        vector<SparseMap<uint64_t, double>> LB;                // Alsó Becslések
        vector<SparseMap<uint64_t, double>> Done;              // Alsó Becslések
        vector<SparseMap<uint64_t, ListDigraph::Arc>> NextArc; // Következő élek
        vector<SparseMap<uint64_t, double>> DP;                // Egzakt értékek

        // LowerBoundhoz feszítőfák
        SparseMap<uint64_t, double> SST{0.0};

        // Branch and Bound fa
        ListDigraph::Node r;        // Held Karp féle gyökér
        vector<tree_node> tasklist; // Listája a A részfeladatoknak amit még ki kell vizsgálni.
        DEBUG logger;               // A logoláshoz használt objektum

#pragma endregion

#pragma region Initialization
        void init(ArcLookUp<ListDigraph> &lookUp)
        {
            // Élek rendezése
            for (ListDigraph::ArcIt a(G); a != INVALID; ++a)
            {
                Arcs.push_back(a);
            }

            std::sort(Arcs.begin(), Arcs.end(),
                      [this](const ListDigraph::Arc &a, const ListDigraph::Arc &b)
                      {
                          return weight[a] < weight[b];
                      });

            // Felső becslés beállítása
            Heuristic::Max_insert maxins{G, Label, weight};
            maxins.Run();
            if (Upper_bound > maxins.Length)
            {
                Upper_bound = maxins.Length;
                Best_Route = maxins.Route;
            }

            Heuristic::Repetitive_Nearest_Neighbour RNN{G, Label, weight};
            if (Upper_bound > RNN.Length)
            {
                Upper_bound = RNN.Length;
                Best_Route = RNN.Route;
            }

            // Tasklist Vector előkészítése
            tasklist.clear();
            uint64_t H = (1ULL << (n - 2)) - 1;
            for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
            {
                if (Label[v] == n - 1)
                {
                    r = v;
                    break;
                }
            }
            for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
            {
                if (Label[v] == n - 1)
                {
                }
                else
                {
                    if (lookUp(v, r) != INVALID)
                    {
                        H &= ~(1ULL << Label[v]);
                        tasklist.emplace_back(false, H, v, lookUp(v, r));
                        logger.log("Added: ", tasklist.back().finalize, " ", tasklist.back().mask, " ", Label[tasklist.back().node], " ", Label[G.source(tasklist.back().arc)], "->", Label[G.target(tasklist.back().arc)]);
                        H |= (1ULL << Label[v]);
                    }
                }
            }
        }
#pragma endregion

#pragma region Route found
        void Route_found_(tree_node &X, tree_node &START)
        {
            if (G.source(Route.back()) != r)
            {
                add(NextArc[Label[X.node]][X.mask]);
                logger.log("Iter on: ", X.mask, " ", Label[X.node], "\n", Label[G.source(NextArc[Label[X.node]][X.mask])], "->", Label[G.target(NextArc[Label[X.node]][X.mask])], "\n");
                uint64_t newmask = X.mask;
                newmask &= ~(1ULL << Label[G.source(NextArc[Label[X.node]][X.mask])]);
                ListDigraph::Node newnode = G.source(NextArc[Label[X.node]][X.mask]);
                tree_node new_X(false, newmask, newnode, NextArc[Label[X.node]][X.mask]);

                Route_found_(new_X, START);
                return;
            }
            if (int(Route.size()) < n)
            {
                std::map<int, bool> Check;
                logger.log("Too short route: ");
                for (size_t i = 0; i < Route.size(); i++)
                {
                    logger.log(Label[G.source(Route[i])], "->", Label[G.target(Route[i])]);
                    Check.insert({Label[G.source(Route[i])], true});
                }

                for (int i = 0; i < n; i++)
                {
                    if (Check.find(i) == Check.end())
                    {
                        logger.log(i, " is not inluded");
                    }
                }
                throw invalid_argument("Route too short");
            }

            std::map<int, bool> Check;
            bool error = false;
            double val = 0.0;
            for (size_t i = 0; i < Route.size(); i++)
            {
                logger.log(Label[G.source(Route[i])], "->", Label[G.target(Route[i])], " ", val, "+", weight[Route[i]]);
                if (Check.find(Label[G.source(Route[i])]) != Check.end())
                {
                    logger.log(Label[G.source(Route[i])], " is duplicate");
                }
                Check.insert({Label[G.source(Route[i])], true});
                val += weight[Route[i]];
            }
            if (error)
            {
                throw invalid_argument("Not a route");
            }

            if (val < Upper_bound)
            {
                Upper_bound = val;
                Best_Route = Route;
            }
        }
        void Route_found(tree_node &X, tree_node &START)
        {
            double tempv = R_val;
            vector<ListDigraph::Arc> temp = Route;
            Route_found_(X, START);
            Route = temp;
            R_val = tempv;
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
            if (!X.mask)
            {
                return 0;
            }
            double rout{Upper_bound + 1}, tin{Upper_bound + 1}, FH{0.0};
            for (ListDigraph::InArcIt a(G, X.node); a != INVALID; ++a)
            {
                if ((X.mask & (1ULL << Label[G.source(a)])) && weight[a] < tin)
                {
                    tin = weight[a];
                }
            }

            for (ListDigraph::OutArcIt a(G, r); a != INVALID; ++a)
            {
                if ((X.mask & (1ULL << Label[G.target(a)])) && weight[a] < rout)
                {
                    rout = weight[a];
                }
            }

            if (SST[X.mask] != 0)
            {
                return SST[X.mask] + rout + tin;
            }
            DSU dsu{n};
            int edge_count = 0;
            for (size_t i = 0; i < Arcs.size(); i++)
            {
                ListDigraph::Arc &a = Arcs[i];
                if ( (X.mask & (1ULL << Label[G.target(a)])) && (X.mask & (1ULL << Label[G.source(a)])))
                {
                    if (dsu.find(Label[G.target(a)]) != dsu.find(Label[G.source(a)]))
                    {
                        edge_count++;
                        FH += weight[a];
                        dsu.unite(Label[G.target(a)], Label[G.source(a)]);
                    }
                }
            }
            if (edge_count < int(__builtin_popcountll(X.mask)) - 1)
            {
                FH = Upper_bound + 1;
            }
            SST[X.mask] = FH;
            return LB[Label[X.node]][X.mask] = FH + rout + tin;
        }
#pragma endregion

#pragma region processnode
        void process(tree_node &X, ArcLookUp<ListDigraph> &lookUp)
        {
            logger.log("Proccesing: ", X.mask, " ", Label[X.node], " Finalize: ", X.finalize, " Depth: ", Route.size());
            // Ha befejezzük
            if (X.finalize)
            {
                logger.log("------------------------------------------");
                pop();
                double minval = Upper_bound + 1;
                ListDigraph::Node tempNextNode;
                ListDigraph::InArcIt a(G, X.node);
                for (ListDigraph::InArcIt a(G, X.node); a != INVALID; ++a)
                {
                    X.mask ^= (1ULL << Label[G.source(a)]);
                    if (minval > DP[Label[G.source(a)]][X.mask] + weight[a])
                    {
                        minval = DP[Label[G.source(a)]][X.mask] + weight[a];
                        tempNextNode = G.source(a);
                    }
                    X.mask ^= (1ULL << Label[G.source(a)]);
                }
                DP[Label[X.node]][X.mask] = LB[Label[X.node]][X.mask] = minval;
                if (minval < Upper_bound + 1)
                {
                    NextArc[Label[X.node]][X.mask] = lookUp(tempNextNode, X.node);
                }
                Done[Label[X.node]][X.mask] = 1;
            }
            // Ha most kezdjük
            else
            {
                // if (weight[X.arc] + R_val >= Upper_bound)
                // {
                //     return;
                // }
                add(X.arc);
                if (!X.mask)
                {
                    if (lookUp(r, X.node) == INVALID)
                    {
                        DP[Label[X.node]][X.mask] = LB[Label[X.node]][X.mask] = Upper_bound + 1;
                        return;
                    }
                    add(lookUp(r, X.node));
                    DP[Label[X.node]][X.mask] = LB[Label[X.node]][X.mask] = weight[lookUp(r, X.node)];
                    NextArc[Label[X.node]][X.mask] = lookUp(r, X.node);
                    Route_found(X, X); // Ha már minden csúcsot meglátogattunk
                    pop();
                    pop();
                    return;
                }
                vector<tree_node> tasks;
                for (ListDigraph::InArcIt a(G, X.node); a != INVALID; ++a)
                {
                    if (G.source(a) == r || !(X.mask & (1ULL << Label[G.source(a)])))
                    {
                        continue;
                    }

                    // Új mask = H - source(a)
                    uint64_t newmask = X.mask;
                    newmask &= ~(1ULL << Label[G.source(a)]);

                    // Új csúcs source(a)
                    ListDigraph::Node newnode = G.source(a);
                    tree_node new_X(false, newmask, newnode, a);

                    // Ha nincs még kész nézzük meg az alsó becslést, csináljuk meg
                    if (Done[Label[newnode]][newmask] == 0)
                    {
                        new_X.LB = Bound(new_X);
                        if (new_X.LB + weight[a] < Upper_bound)
                        {
                            tasks.push_back(new_X);
                        }
                    }

                    // Ha kész van: vagy az arc INVALID azaz nincs jó út, vagy van egy utunk
                    else
                    {

                        if (NextArc[Label[newnode]][newmask] != INVALID && DP[Label[newnode]][newmask] + weight[a] + R_val < Upper_bound)
                        {
                            logger.log("Route end: ", newmask, " ", Label[newnode]);
                            add(a);
                            Route_found(new_X, new_X); // Ha van egy utunk, és kész van előtte minden
                            pop();
                        }
                    }
                }
                std::sort(tasks.begin(), tasks.end(), [](const tree_node &a, const tree_node &b)
                          {
                              return a.LB > b.LB; // Csökkenő sorrend
                          });
                X.finalize = true;
                logger.log("Added: ", X.finalize, " ", X.mask, " ", Label[X.node], " ", Label[G.source(X.arc)], "->", Label[G.target(X.arc)]);
                tasklist.push_back(X);
                R_update();
                for (size_t i = 0; i < tasks.size(); i++)
                {
                    string rute = "";
                    for (size_t i = 0; i < Route.size(); i++)
                    {
                        auto &a = Route[i];
                        string arc = to_string(Label[G.source(a)]) + "->" + to_string(Label[G.target(a)]) + " ";
                        rute += arc;
                    }
                    if (tasks[i].LB + weight[tasks[i].arc] < Upper_bound)
                    {
                        logger.log("Added: ", tasks[i].finalize, " ", tasks[i].mask, " ", Label[tasks[i].node], " ", Label[G.source(tasks[i].arc)], "->", Label[G.target(tasks[i].arc)]);
                        tasklist.push_back(tasks[i]);
                    }
                    else
                    {
                        logger.log(rute, " ", tasks[i].mask, " ", Label[tasks[i].node], "val: ", R_val, "+", tasks[i].LB, "=", R_val + tasks[i].LB, ">=", Upper_bound);
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
            const double _Upper_bound = std::numeric_limits<double>::max()) : G(_G), Label(_Label), weight(_weight), Upper_bound(_Upper_bound)
        {
            n = countNodes(_G);
            V.resize(n);
            for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
            {
                V[_Label[i]] = i;
            }
            logger.log("Nodes: ", n, " Upperbound: ", Upper_bound, "\n");

            // 1. LB inicializálása, ha még nem számoltunk ki jobbat legyen 0
            for (int i = 0; i < n; i++)
            {
                SparseMap<uint64_t, double> lb{0.0};
                LB.push_back(lb);
            }

            // 1.2. Done inicializálása, ha még nem számoltuk ki 0
            for (int i = 0; i < n; i++)
            {
                SparseMap<uint64_t, double> dn{0.0};
                Done.push_back(dn);
            }

            // 2. NextArc inicializálása INVALID alapból
            for (int i = 0; i < n; i++)
            {
                SparseMap<uint64_t, ListDigraph::Arc> Arcs{INVALID};
                NextArc.push_back(Arcs);
            }

            // 3. Val inicializálása legyen alapból egy nagy szám
            for (int i = 0; i < n; i++)
            {
                SparseMap<uint64_t, double> vals{Upper_bound + 1};
                DP.push_back(vals);
            }
        };
#pragma endregion

#pragma region solve()
        double solve()
        {
            lemon::ArcLookUp<ListDigraph> lookUp(G);
            init(lookUp);
            logger.log("Root: ", Label[r]);
            while (!tasklist.empty())
            {
                tree_node X = tasklist.back();
                tasklist.pop_back();
                process(X, lookUp);
            }
            solved = 1;
            Current_bestval = Upper_bound;
            return Current_bestval;
        }

#pragma endregion
#pragma region Query functions
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
            return Best_Route;
        }

        void printroute()
        {
            vector<ListDigraph::Arc> Tour = this->OPTroute();
            for (size_t i = 0; i < Tour.size(); i++)
            {
                cout << "->" << Label[(G.target(Tour[i]))];
            }
        }
#pragma endregion
    };
}
