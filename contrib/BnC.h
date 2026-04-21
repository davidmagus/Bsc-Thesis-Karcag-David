#pragma once
#pragma region Includes, Namespaces
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/time_measure.h>
#include <stdexcept>
#include <lemon/lp.h>
#include <type_traits>
#include <list>
#include <lemon/preflow.h>
#include <lemon/adaptors.h>
#include "Heuristic.h"
using namespace std;
using namespace lemon;
#pragma endregion

#pragma region Bitmask extensions
namespace Bit
{
    inline bool is(const boost::dynamic_bitset<> &bits, size_t i)
    {
        return (i < bits.size()) && bits[i];
    }

    inline boost::dynamic_bitset<> AND(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b)
    {
        size_t max_size = std::max(a.size(), b.size());

        if (a.size() < max_size)
            a.resize(max_size);
        if (b.size() < max_size)
            b.resize(max_size);

        return a & b;
    }

    inline boost::dynamic_bitset<> OR(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b)
    {
        size_t max_size = std::max(a.size(), b.size());
        if (a.size() < max_size)
            a.resize(max_size);
        if (b.size() < max_size)
            b.resize(max_size);

        return a | b;
    }
}
#pragma endregion

namespace BnCnP
{

#pragma region Logging Tools
    struct Logging
    {
        std::ofstream logFile;
        int counter;
        Logging() : logFile("BnC.log") { counter = 1500; }

        template <typename... Args>
        void log(Args... args)
        {
            if (counter)
            {
                counter--;
                (logFile << ... << args);
                logFile << endl;
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

    template <typename DEBUG = Silent> // Template parameters: <DEBUG = STSP::Silent / STSP::Logging
    class Algorithm
    {
    private:
#pragma region Member variables
        // Kapott adatok
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        vector<ListDigraph::Node> V;
        const ListDigraph::ArcMap<double> &weight;
        ListDigraph::ArcMap<bool> IsCol;
        double Upper_bound;

        // belso Eszkozok
        bool OPTsolved = false;
        int solved = 0; // A query függvények számára fentartott belső ellenörző 0 ha nem hívták még meg a .solve metódust 1 egyébként.
        DEBUG logger;   // A logoláshoz használt objektum

        // LP / IP Eszközök
        Lp coreLP;
        vector<std::pair<ListDigraph::Arc, Lp::Col>> Columns; // Oszlopok
        list<ListDigraph::Arc> Arcs_notincluded;

        // Paraméterek
        const int col_add_number;
        const int sepequality_add_number;
        const int exploring_steps;
        const size_t initarcs;
        const double timeout;

        // Branch and cut Eszközök
        int idgiver = 0;
        double Best_Val = 0;
        vector<ListDigraph::Arc> Best_Tour;

#pragma endregion

#pragma region BnC Tree Nodes
        struct node_BnCnPtree
        {
            int Parent_ID;
            boost::dynamic_bitset<> Edges_setto_0;
            boost::dynamic_bitset<> Edges_setto_1;
            double LB;
            int ID;
            node_BnCnPtree(
                Algorithm &O,
                const int _Pid,                              // ID of the Parent Task
                const boost::dynamic_bitset<> edges_setto_0, // Edges for which the value is fixed to 1
                const boost::dynamic_bitset<> edges_setto_1, // Edges for which the value is fixed to 0
                const double LB,                             // Lower Bound
                const int _ID                                // ID
                ) : Parent_ID(_Pid),
                    Edges_setto_0(edges_setto_0),
                    Edges_setto_1(edges_setto_1),
                    LB(LB),
                    ID(_ID)
            {
                O.logger.log("A new task added! ", "ID: ", ID, " Parents ID: ", Parent_ID);
            }
        };
        list<node_BnCnPtree> nodes; // A fegoldolgozatlan nodeok listája
#pragma endregion

#pragma region Inequalities
        struct degree_eq
        {
            ListDigraph::Node v;
            int direction;
            Lp::Row eq;
            degree_eq(Algorithm &O, ListDigraph::Node _v, int _direction) : v(_v), direction(_direction)
            {
                eq = O.coreLP.addRow(1, 0, 1);
                O.logger.log("Degree equation added: ", O.Label[v], " Direction: ", direction);
                for (std::pair<ListDigraph::Arc, Lp::Col> e : O.Columns)
                {
                    addcol(O, e);
                }
            }
            void addcol(Algorithm &O, std::pair<ListDigraph::Arc, Lp::Col> e)
            {
                switch (direction)
                {
                case 1:
                    if (O.G.target(e.first) == v)
                    {
                        O.coreLP.coeff(eq, e.second, 1);
                    }
                    break;

                case 0:
                    if (O.G.source(e.first) == v)
                    {
                        O.coreLP.coeff(eq, e.second, 1);
                    }
                    break;
                }
            }
            int coeff(Algorithm &O, ListDigraph::Arc &a)
            {
                switch (direction)
                {
                case 1:
                    if (O.G.target(a) == v)
                    {
                        return 1;
                    }
                    else
                    {
                        return 0;
                    }
                    break;

                case 0:
                    if (O.G.source(a) == v)
                    {
                        return 1;
                    }
                    else
                    {
                        return 0;
                    }
                    break;
                }
                return 0;
            }
        };
        vector<degree_eq> deg_eqs; // A fokszám egyenlőségek
        struct separation_ineq
        {
            boost::dynamic_bitset<> inside;
            Lp::Row eq;
            separation_ineq(Algorithm &O, boost::dynamic_bitset<> &_inside) : inside(_inside)
            {
                eq = O.coreLP.addRow(0, 0, inside.count() - 1);
                for (std::pair<ListDigraph::Arc, Lp::Col> e : O.Columns)
                {
                    addcol(O, e);
                }

                O.logger.log("Separation equation added for: ", inside);
            }
            void addcol(Algorithm &O, std::pair<ListDigraph::Arc, Lp::Col> e)
            {
                if (inside[O.Label[O.G.target(e.first)]] && inside[O.Label[O.G.source(e.first)]])
                {
                    O.coreLP.coeff(eq, e.second, 1);
                }
            }
            int coeff(Algorithm &O, ListDigraph::Arc &a)
            {
                if (inside[O.Label[O.G.target(a)]] && inside[O.Label[O.G.source(a)]])
                {
                    return 1;
                }
                else
                {
                    return 0;
                }
            }
        };
        vector<separation_ineq> sep_ineqs; // A SEC egyenlőtlenségek listája

#pragma endregion

#pragma region Adding new column
        void addcol(ListDigraph::Arc a)
        {
            if (!IsCol[a])
            {
                Lp::Col Col = coreLP.addCol();
                std::pair<ListDigraph::Arc, Lp::Col> e{a, Col};
                Columns.push_back(e);
                for (degree_eq i : deg_eqs)
                {
                    i.addcol(*this, e);
                }
                for (separation_ineq i : sep_ineqs)
                {
                    i.addcol(*this, e);
                }
                IsCol[a] = true;
                coreLP.objCoeff(Col, weight[a]);
                coreLP.colBounds(Col, 0, 1);
                logger.log("Column added for arc: (", Label[G.source(a)], ",", Label[G.target(a)], ")");
            }
        }
#pragma endregion

#pragma region Initialization
        void init()
        {
            logger.log("\nIntialization...");
            logger.log("\nClosest Neighbour Arcs: ");
            Heuristic::Max_insert CL{G, Label, weight};
            CL.Run();
            for (ListDigraph::Arc a : CL.Route)
            {
                addcol(a);
            }
            Best_Val = CL.Length;
            cout << Best_Val;
            Best_Tour = CL.Route;
            logger.log("Closest Neighbour Value: ", Best_Val);
            logger.log("\nShortest Arcs: ");
            vector<ListDigraph::Arc> Arcs; // Adjuk hozzá a legolcsóbb kimenő éleket minden pontból.
            for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
            {
                std::vector<ListDigraph::Arc> arcsL;

                // 1. Összes kimenő él összegyűjtése
                for (ListDigraph::OutArcIt a(G, v); a != INVALID; ++a)
                {
                    if (!IsCol[a])
                    {
                        arcsL.push_back(a);
                    }
                }

                // 2. Ha több van mint 15, kiválasztjuk a legkisebbeket
                if (arcsL.size() > initarcs)
                {
                    std::nth_element(arcsL.begin(), arcsL.begin() + initarcs, arcsL.end(),
                                     [&](const ListDigraph::Arc &left, const ListDigraph::Arc &right)
                                     {
                                         return weight[left] < weight[right]; // Kisebb súly előre
                                     });
                    // Csak az első  maradjon
                    arcsL.resize(initarcs);
                }

                // 3. Hozzáadás a globális listához
                for (auto &a : arcsL)
                {
                    Arcs.push_back(a);
                }
            }
            for (ListDigraph::Arc a : Arcs)
            {
                addcol(a);
            }

            logger.log("\nFilling list of not included arcs");
            for (ListDigraph::ArcIt a(G); a != INVALID; ++a)
            {
                if (!IsCol[a])
                {
                    Arcs_notincluded.push_back(a);
                }
            }

            logger.log("\nAdding starting constraints: ");
            ListDigraph::NodeIt v(G);

            // Első csúcsnál csak az egyik kell
            deg_eqs.emplace_back(*this, v, 0);
            ++v;

            // A többinél mindkettő
            while (v != INVALID)
            {
                deg_eqs.emplace_back(*this, v, 0);
                deg_eqs.emplace_back(*this, v, 1);
                ++v;
            }
            boost::dynamic_bitset<> mask{Columns.size()};
            nodes.emplace_front(*this, -1, mask, mask, 0, ++idgiver);
            coreLP.min();
            logger.log("\nEnd of Intialization");
        }

#pragma endregion
#pragma region Colum generation
        void Pricing()
        {
            int counter = 0;
            for (ListDigraph::Arc a : Arcs_notincluded)
            {
                double Dualcost = -weight[a];
                for (degree_eq i : deg_eqs)
                {
                    Dualcost += coreLP.dual(i.eq) * i.coeff(*this, a);
                }
                for (separation_ineq i : sep_ineqs)
                {
                    Dualcost += coreLP.dual(i.eq) * i.coeff(*this, a);
                }
                if (Dualcost >= 0)
                {
                    counter++;
                    addcol(a);
                }
                if (col_add_number && counter >= col_add_number)
                {
                    break;
                }
            }
            Arcs_notincluded.remove_if([&](const ListDigraph::Arc &a)
                                       { return IsCol[a]; });
        }
#pragma endregion

#pragma region Separation
        bool Separation() // Returns true if separation equalities have been found
        {
            FilterArcs<const ListDigraph, ListDigraph::ArcMap<bool>> Flowgraph(G, IsCol);
            ListDigraph::ArcMap<double> primalsol(G);

            for (std::pair<ListDigraph::Arc, Lp::Col> a : Columns)
            {
                primalsol[a.first] = coreLP.primal(a.second);
            }

            vector<boost::dynamic_bitset<>> cuts;

            ListDigraph::NodeIt k(G);

            ListDigraph::Node s = k;

            ++k;

            while (k != INVALID)
            {

                typedef FilterArcs<const ListDigraph, ListDigraph::ArcMap<bool>> FilteredGraph;
                typedef ListDigraph::ArcMap<double> CapacityMap;
                Preflow<FilteredGraph, CapacityMap> st(Flowgraph, primalsol, s, static_cast<ListDigraph::Node>(k));
                st.runMinCut();

                if (st.flowValue() < 1)
                {
                    boost::dynamic_bitset<> temp(n);
                    for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
                    {
                        if (st.minCut(v))
                        {
                            temp.set(Label[v]);
                        }
                    }
                    bool added = false;
                    for (boost::dynamic_bitset<> i : cuts)
                    {
                        if (temp == i)
                        {
                            added = true;
                            break;
                        }
                    }
                    if (!added)
                    {
                        cuts.push_back(temp);
                    }
                }
                ++k;
            }

            for (boost::dynamic_bitset<> S : cuts)
            {
                sep_ineqs.emplace_back(*this, S);
            }

            return !cuts.empty();
        }
#pragma endregion

#pragma region Branching
        void Branch(size_t i, node_BnCnPtree &X)
        {
            ListDigraph::Node t = G.target(Columns[i].first);

            vector<size_t> Candidate_arcs;
            vector<size_t> Unused_tarcs;
            for (size_t j = 0; j < Columns.size(); j++)
            {
                if (G.target(Columns[j].first) == t)
                {
                    if (coreLP.primal(Columns[j].second) > 0.3)
                    {
                        Candidate_arcs.push_back(j);
                    }
                    else
                    {
                        Unused_tarcs.push_back(j);
                    }
                }
            }
            boost::dynamic_bitset<> new0s = X.Edges_setto_0;
            boost::dynamic_bitset<> new1s = X.Edges_setto_1;
            new0s.resize(Columns.size());
            new1s.resize(Columns.size());

            for (size_t k : Candidate_arcs)
            {
                new0s.set(k);
            }
            nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver); // Where all the arcs with high value are set to 0.
            for (size_t k : Unused_tarcs)
            {
                new0s.set(k);
            }

            for (size_t i : Candidate_arcs) // Where a Single high value arc is set to 1, and all others to 0.
            {
                new1s.set(i);
                new0s.reset(i);
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver);
                new1s.reset(i);
                new0s.set(i);
            }
        }
#pragma endregion

#pragma region Processing a Node
        void process_Node(node_BnCnPtree &X)
        {
            logger.log("\nProcessing Node: ", X.ID, " Lower Bound: ", X.LB, " Current Upper Bound: ", Best_Val);
            // Setting up Lp
            for (size_t i = 0; i < Columns.size(); i++)
            {
                bool set = false;
                if (Bit::is(X.Edges_setto_0, i))
                {
                    coreLP.colBounds(Columns[i].second, 0, 0);
                    set = true;
                }
                if (!set && Bit::is(X.Edges_setto_1, i))
                {
                    coreLP.colBounds(Columns[i].second, 1, 1);
                    set = true;
                }
                if (!set)
                {
                    coreLP.colBounds(Columns[i].second, 0, 1);
                }
            }

            // Pricing step
            coreLP.solve();
            Pricing();

            // Now we can Compute a Lowerbound
            coreLP.solve();
            X.LB = coreLP.primal();
            if (!(X.LB <= Best_Val)) // Cutting step
            {
                return;
            }

            // Separation step
            bool found_violated_eqs = Separation();

            if (found_violated_eqs)
            {
                nodes.emplace_front(*this, X.ID, X.Edges_setto_0, X.Edges_setto_1, X.LB, ++idgiver);
                return;
            }

            // If the solution is feasablie constraintwise we start fixing non intiger arcs as integers
            size_t i = 0;
            while (i < Columns.size())
            {
                if (std::abs(coreLP.primal(Columns[i].second) - std::round(coreLP.primal(Columns[i].second))) < 1e-9)
                {
                    i++;
                }
                else
                {
                    break;
                }
            }
            if (i < Columns.size())
            {
                Branch(i, X);
            }
            else
            {
                Tourfound();
            }
        }
#pragma endregion
    public:
#pragma region Constructor
        Algorithm(
            const ListDigraph &_G,
            const ListDigraph::NodeMap<int> &_Label,
            const ListDigraph::ArcMap<double> &_weight,
            const double _Upper_bound = std::numeric_limits<double>::max(),
            const int col = 0,
            const int sep = 5,
            const int exploring_steps = 100,
            const size_t _initarcs = 10,
            const double timelimt = 600) : G(_G),
                                           Label(_Label),
                                           weight(_weight),
                                           IsCol(_G, false),
                                           Upper_bound(_Upper_bound),
                                           col_add_number(col),
                                           sepequality_add_number(sep),
                                           exploring_steps(exploring_steps),
                                           initarcs(_initarcs),
                                           timeout(timelimt)
        {
            n = countNodes(_G);
            V.resize(n);
            for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
            {
                V[_Label[i]] = i;
            }
            logger.log("Nodes: ", n, " Upperbound: ", Upper_bound, "\n");
            Best_Val = 0;
        };
#pragma endregion

#pragma region Solve
        double solve()
        {
            Timer t;
            init();

            // Exploring Phase
            for (int i = 0; i < exploring_steps; i++)
            {
                if (nodes.empty() || t.realTime() > timeout)
                {
                    break;
                }

                auto it = std::min_element(nodes.begin(), nodes.end(),
                                           [](const node_BnCnPtree &a, const node_BnCnPtree &b)
                                           {
                                               return a.LB < b.LB;
                                           });

                node_BnCnPtree X = *it;
                nodes.erase(it);

                process_Node(X);
            }

            while (!nodes.empty() && t.realTime() < timeout)
            {

                node_BnCnPtree X = nodes.front();
                nodes.pop_front();
                process_Node(X);
            }
            if (t.realTime() < timeout)
            {
                OPTsolved = true;
            }
            solved = 1;
            return Best_Val;
        }

#pragma endregion

#pragma region Tour Found
        void Tourfound()
        {
            logger.log("A new Tour found with value: ", coreLP.primal());
            if (coreLP.primal() > Best_Val)
            {
                return;
            }
            vector<ListDigraph::Arc> Tour;
            double val = 0;
            for (std::pair<ListDigraph::Arc, Lp::Col> e : Columns)
            {
                if (std::round(coreLP.primal(e.second)) == 1)
                {
                    Tour.push_back(e.first);
                    val += weight[e.first];
                }
            }
            Best_Tour = Tour;
            Best_Val = val;

            nodes.remove_if([&](const node_BnCnPtree &v)
                            { return (v.LB >= Best_Val); });
        }
#pragma endregion

#pragma region Query Functions
        bool get_OPTsolved() { return OPTsolved; }
        double OPTval()
        {
            if (!solved)
            {
                throw std::invalid_argument("Most call method .solve() before using query functions");
            }
            return Best_Val;
        }

        vector<ListDigraph::Arc> OPTroute()
        {
            if (!solved)
            {
                throw std::invalid_argument("Most call method .solve() before using query functions");
            }
            vector<ListDigraph::Arc> Tour;
            ListDigraph::Node O;
            for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
            {
                if (Label[v] == 0)
                {
                    O = v;
                    break;
                }
            }
            ListDigraph::Node Nextnode;
            for (ListDigraph::Arc a : Best_Tour)
            {
                if (G.source(a) == O)
                {
                    Nextnode = G.target(a);
                    Tour.push_back(a);
                }
            }

            bool notDone = true;
            while (notDone)
            {
                for (ListDigraph::Arc a : Best_Tour)
                {

                    if (G.source(a) == Nextnode)
                    {
                        Nextnode = G.target(a);
                        Tour.push_back(a);
                        if (Nextnode == O)
                        {
                            notDone = false;
                        }
                    }
                }
            }
            Best_Tour = Tour;
            return Best_Tour;
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