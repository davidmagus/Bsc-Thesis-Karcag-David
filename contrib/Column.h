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

namespace Column_generation
{
    class Algorithm
    {
    private:
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

        // LP / IP Eszközök
        Lp coreLP;
        vector<std::pair<ListDigraph::Arc, Lp::Col>> Columns; // Oszlopok
        list<ListDigraph::Arc> Arcs_notincluded;

        const size_t initarcs;
        const double timeout;

        // Branch and cut Eszközök
        int idgiver = 0;
        double Best_Val = 0;
        vector<ListDigraph::Arc> Best_Tour;

        #pragma region Inequalities
        struct degree_eq
        {
            ListDigraph::Node v;
            int direction;
            Lp::Row eq;
            degree_eq(Algorithm &O, ListDigraph::Node _v, int _direction) : v(_v), direction(_direction)
            {
                eq = O.coreLP.addRow(1, 0, 1);
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
#pragma region Initialization

        void init()
        {
            Heuristic::Nearest_Neighbour CL{G, Label, weight};
            for (ListDigraph::Arc a : CL.Route)
            {
                addcol(a);
            }
            Best_Val = CL.Length;
            Best_Tour = CL.Route;
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

            for (ListDigraph::ArcIt a(G); a != INVALID; ++a)
            {
                if (!IsCol[a])
                {
                    Arcs_notincluded.push_back(a);
                }
            }

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
        }
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
            const size_t _initarcs = 10,
            const double timelimt = 600) : G(_G),
                                           Label(_Label),
                                           weight(_weight),
                                           IsCol(_G, false),
                                           Upper_bound(_Upper_bound),
                                           initarcs(_initarcs),
                                           timeout(timelimt)
        {
            n = countNodes(_G);
            V.resize(n);
            for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
            {
                V[_Label[i]] = i;
            }
            Best_Val = 0;
        };
#pragma endregion
    };
}