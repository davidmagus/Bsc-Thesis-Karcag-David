#pragma once
constexpr double EPS = 1e-9; // Vagy 1e-7, a feladat pontosságától függően
#pragma region Includes_Namespaces
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <lemon/list_graph.h>
#include <lemon/time_measure.h>
#include <stdexcept>
#include <lemon/lp.h>
#include <list>
#include <unordered_set>
#include <lemon/preflow.h>
#include <lemon/adaptors.h>
#include "Heuristic.h"
using namespace std;
using namespace lemon;
#pragma endregion

#pragma region Bitmask_extensions
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
#pragma region Logging
    struct Logging
    {
        std::ofstream logFile;
        int counter;
        Logging() : logFile("analysis/Logs/BnC.log") { counter = 5000; }

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
        static void log(Args...)
        {
        }
    };

#pragma endregion

    template <typename DEBUG = Silent, int BRANCHING = 1> // Template parameters: <DEBUG = STSP::Silent / STSP::Logging
    class Algorithm
    {
    private:
#pragma region Member_variables
        // Kapott adatok
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        vector<ListDigraph::Node> V;
        const ListDigraph::ArcMap<double> &weight;
        ListDigraph::ArcMap<bool> IsCol;
        ArcLookUp<ListDigraph> Lookup;
        double Upper_bound;

        // belso Eszkozok
        bool OPTsolved = false;
        int solved = 0; // A query függvények számára fentartott belső ellenörző 0 ha nem hívták még meg a .solve metódust 1 egyébként.
        DEBUG logger;   // A logoláshoz használt objektum

        // LP / IP Eszközök
        Lp coreLP;
        vector<std::pair<ListDigraph::Arc, Lp::Col>> Columns; // Oszlopok
        list<ListDigraph::Arc> Arcs_notincluded;
        ListDigraph::ArcMap<Lp::Col> Col;

        // Paraméterek
        const int sepequality_add_number;
        const int col_add_number;
        const int exploring_steps;
        size_t initarcs;
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
                const boost::dynamic_bitset<>& edges_setto_0, // Edges for which the value is fixed to 1
                const boost::dynamic_bitset<>& edges_setto_1, // Edges for which the value is fixed to 0
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
            degree_eq(Algorithm &O, ListDigraph::Node _v, int _direction, bool silent = true) : v(_v), direction(_direction)
            {
                eq = O.coreLP.addRow(1, 0, 1);
                if (!silent)
                {
                    O.logger.log("Degree equation added: ", O.Label[v], " Direction: ", direction);
                }
                for (std::pair<ListDigraph::Arc, Lp::Col> e : O.Columns)
                {
                    addcol(O, e);
                }
            }
            int coeff(const Algorithm &O, const ListDigraph::Arc &a)
            {
                switch (direction)
                {
                case 1:
                    if (O.G.target(a) == v)
                    {
                        return 1;
                    }
                    return 0;
                    break;

                case 0:
                    {
                        if (O.G.source(a) == v)
                        {
                            return 1;
                        }
                        return 0;
                    }
                    break;
                default:
                    throw std::runtime_error("Unknown direction");
                }
                return 0;
            }
            void addcol(Algorithm &O, const std::pair<ListDigraph::Arc,const Lp::Col>& e)
            {
                if (coeff(O, e.first))
                {
                    O.coreLP.coeff(eq, e.second, 1);
                }
            }
        };
        vector<degree_eq> deg_eqs; // A fokszám egyenlőségek
        struct separation_ineq
        {
            const boost::dynamic_bitset<> inside;
            Lp::Row eq;
            separation_ineq(Algorithm &O, const boost::dynamic_bitset<> &_inside) : inside(_inside)
            {
                eq = O.coreLP.addRow(0, 0, inside.count() - 1);
                for (std::pair<ListDigraph::Arc, Lp::Col> e : O.Columns)
                {
                    addcol(O, e);
                }

                O.logger.log("sec_separation equation added for: ", inside);
            }
            int coeff(const Algorithm &O, const ListDigraph::Arc &a)
            {
                if (inside[O.Label[O.G.target(a)]] && inside[O.Label[O.G.source(a)]])
                {
                    return 1;
                }
                return 0;
            }
            void addcol(Algorithm &O, const std::pair<ListDigraph::Arc, Lp::Col>& e)
            {
                if (coeff(O, e.first))
                {
                    O.coreLP.coeff(eq, e.second, 1);
                }
            }
        };
        vector<separation_ineq> sep_ineqs; // A SEC egyenlőtlenségek listája

        class Comb_ineqs // Az alternatív alakot használjuk
        {
            public:
            const boost::dynamic_bitset<> relevant;       // Ez a halmaz tartalmaz minden csúcsot ami benne van bármelyik halmazban
            const boost::dynamic_bitset<> H;              // H csúcsai
            const vector<boost::dynamic_bitset<>> Teeth;  // Fogak csúcsai
            Lp::Row eq;                                   // Hozzá tartozó LP sor
            Comb_ineqs(Algorithm &O, const boost::dynamic_bitset<>& relevant, const boost::dynamic_bitset<>& h,
                const vector<boost::dynamic_bitset<>>& teeth)
                : relevant(relevant),
                  H(h),
                  Teeth(teeth)
            {
                eq = O.coreLP.addRow(Teeth.size() * 3 + 1, 0, INFINITY); // nagyobb mint 6k + 4 ami 3t + 1
                O.logger.log("Comb inequality added for: ", H, " Number of teeth: ", Teeth.size());
                for (std::pair<ListDigraph::Arc, Lp::Col> e : O.Columns)
                {
                    addcol(O, e);
                }
            }
            int coeff(const Algorithm &O, const ListDigraph::Arc &a)
            {
                if (const int s = O.Label[O.G.source(a)], t = O.Label[O.G.target(a)]; relevant[s] || relevant[t])
                // ha legalább az egyik csúcs releváns átvizsgáljuk.
                {
                    int cnt = 0;
                    // Legfeljebb háromszor lehet a az egyenletben egyszer H miatt és mehet az él két fog között

                    auto tooth = Teeth.begin();
                    while(cnt < 2 && tooth != Teeth.end()) // Csak kettő lehet a fogakból
                    {
                        if ( (*tooth)[s] != (*tooth)[t]) // Azaz adott fogban csak az egyik van benne meaning eleme d(T)
                        {
                            ++cnt;
                        }
                        ++tooth;
                    }

                    if (H[s] != H[t]) //Nézzük meg H-t is
                    {
                        ++cnt;
                    }

                    return cnt;
                }
                return 0;
            }
            void addcol(Algorithm &O, const std::pair<ListDigraph::Arc, Lp::Col>& e)
            {
                if (int cf = coeff(O, e.first))
                {
                    O.coreLP.coeff(eq, e.second, cf);
                }
            }
        };

#pragma endregion

#pragma region column_adding
        void addcol(ListDigraph::Arc a, bool silent = true)
        {
            if (!IsCol[a])
            {
                Lp::Col Col_ = coreLP.addCol();
                std::pair<ListDigraph::Arc, Lp::Col> e{a, Col_};
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
                Col[a] = Col_;
                coreLP.objCoeff(Col_, weight[a]);
                coreLP.colBounds(Col_, 0, 1);
                if (!silent)
                {
                    logger.log("Column added for arc: (", Label[G.source(a)], ",", Label[G.target(a)], ")");
                }
            }
        }
#pragma endregion

#pragma region Initialization //Legyen tobb heu
        void init()
        {
            logger.log("\nIntialization...");
            logger.log("\nClosest Neighbour Arcs: ");
            Heuristic::Max_insert CL{G, Label, weight};
            CL.Run();
            for (ListDigraph::Arc a : CL.Route)
            {
                addcol(a, true);
            }
            Best_Val = CL.Length;
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
                addcol(a, true);
            }

            logger.log("\nFilling list of not included arcs");
            for (ListDigraph::ArcIt a(G); a != INVALID; ++a)
            {
                if (!IsCol[a])
                {
                    Arcs_notincluded.push_back(a);
                }
            }

            logger.log("\nAdding starting constraints ");
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

#pragma region Pricing // Ha félbehagyjuk onnan folytassuk
        bool Pricing() // A bool értéke azt jelzi kell-e a tovább folytatni az oszlopgenerálást, ha true akkor igen
        {
            struct ActiveConstraint {
                double dual_val;
                void* ptr;
                enum Type { DEG, SEP } type;
            };
            int counter{0};
            std::vector<ActiveConstraint> active_constraints;
            active_constraints.reserve(deg_eqs.size() + sep_ineqs.size());

            for (auto& i : deg_eqs) {
                double d = coreLP.dual(i.eq);
                if (d != 0) {
                    active_constraints.push_back({d, static_cast<void*>(&i), ActiveConstraint::DEG});
                }
            }
            for (auto& i : sep_ineqs) {
                double d = coreLP.dual(i.eq);
                if (d != 0) {
                    active_constraints.push_back({d, static_cast<void*>(&i), ActiveConstraint::SEP});
                }
            }

            for (ListDigraph::Arc a : Arcs_notincluded) {
                double Dualreducedcost = weight[a];

                for (const auto& constr : active_constraints) {
                    if (constr.type == ActiveConstraint::DEG) {
                        auto* p = static_cast<degree_eq*>(constr.ptr);
                        Dualreducedcost -= constr.dual_val * p->coeff(*this, a);
                    } else {
                        auto* p = static_cast<separation_ineq*>(constr.ptr);
                        Dualreducedcost -= constr.dual_val * p->coeff(*this, a);
                    }
                }
                if (Dualreducedcost < 0) {
                    counter++;
                    if (counter <= col_add_number)
                    {
                        addcol(a);
                    }
                }
            }
            Arcs_notincluded.remove_if([&](const ListDigraph::Arc &a)
                                       { return IsCol[a]; });
            if (counter == 0)
            {
            return false;  //Nem találtunk sértő oszlopot azaz a jelenlegi megoldás optimális.
            }
            logger.log("Generating columns");
            return true;   //Ha találtunk akár 1-et is akkor tovább kell iterálni.
        }
#pragma endregion

#pragma region SEC_Separation
        bool sec_separation() // Returns true if separation equalities have been found
        {
            FilterArcs<const ListDigraph, ListDigraph::ArcMap<bool>> Flowgraph(G, IsCol);
            ListDigraph::ArcMap<double> primalsol(G);

            for (std::pair<ListDigraph::Arc, Lp::Col> a : Columns)
            {
                primalsol[a.first] = coreLP.primal(a.second);
            }

            std::unordered_set<boost::dynamic_bitset<>> cuts;

            ListDigraph::NodeIt k(G);

            ListDigraph::Node s = k;

            ++k;
            typedef FilterArcs<const ListDigraph, ListDigraph::ArcMap<bool>> FilteredGraph;
            typedef ListDigraph::ArcMap<double> CapacityMap;
            Preflow<FilteredGraph, CapacityMap> st(Flowgraph, primalsol, s, static_cast<ListDigraph::Node>(k));
            while (k != INVALID)
            {
                st.target(static_cast<ListDigraph::Node>(k));
                st.runMinCut();
                if (st.flowValue() < 1.0 - EPS)
                {
                    boost::dynamic_bitset<> temp(n);
                    for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
                    {
                        if (st.minCut(v))
                        {
                            temp.set(Label[v]);
                        }
                    }
                    if (cuts.find(temp) == cuts.end())
                    {
                        cuts.insert(temp);
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

#pragma  region Comb_Separation
        class maxback_Handle_growing
        {
            public:
            Algorithm& O;
            vector<int> maxbackval; // A H-hoz tartozó max_back érték, az n. indexen mindig 0
            // Egy itarátorokból álló lista aminek az i. helyén egy i-re mutató iterátor áll.
            double coboundary = 0;      // H-ból kivezető összérték
            vector<int> pos;
            //Denotes where the current node is. 0: inside 1: subto1 2: infto1 3: extremof1

            vector<int> is_odd;         // Ha nem odd csúcs 0, egyébként hány csúcs teszi oddá;
            vector<int> pointer_to_outside; // Melyik csúcs a kinti 1-el látott szomszéd, ha is_odd > 0
            vector<int> pointer_to_odd; // Általában O.n, ha a csúcs extrém akkor a e helyén az van kit tesz odd-á
            int size = 0; // Size of the currently built Handle.

            maxback_Handle_growing(Algorithm& O, int i, int j) : O(), maxbackval(O.n + 1, 0), pos(O.n, 2),
                                                                 is_odd(O.n+1, 0), pointer_to_odd(O.n, O.n)
            {
                //Init
                pos[i] = 0;
                pos[j] = 3;

                for (ListDigraph::OutArcIt a(O.G, O.V[i]); a != INVALID; ++a)
                {
                    double x_a = 0;
                    if (!O.IsCol[a] || (x_a = O.coreLP.primal(O.Col[a])) == 0)
                    {
                        continue; //Ha $a$ nem oszlop vagy x_a = 0 marad 0
                    }
                    const int t = O.Label[O.G.target(static_cast<ListDigraphBase::Arc>(a))];
                    maxbackval[t] += x_a;
                    if (x_a == 1)
                    {
                        pointer_to_outside[i] = t;      // i-nek mostantól t a kinti 1-el látott szomszédja
                        ++is_odd[i];
                        --is_odd[pointer_to_odd[t]];
                        pointer_to_odd[t] = i;
                        pos[t] = 3;
                    }
                    coboundary += x_a;
                }

                for (ListDigraph::InArcIt a(O.G, O.V[i]); a != INVALID; ++a)
                {
                    double x_a = 0;
                    if (!O.IsCol[a] || (x_a = O.coreLP.primal(O.Col[a])) == 0)
                    {
                        continue; //Ha $a$ nem oszlop vagy x_a = 0 marad 0
                    }
                    const int s = O.Label[O.G.source(static_cast<ListDigraphBase::Arc>(a))];
                    maxbackval[s] += x_a;
                    if (x_a == 1)
                    {
                        pos[s] = 3;
                        pointer_to_outside[i] = s;      // i-nek mostantól s a kinti 1-el látott szomszédja
                        ++is_odd[i];
                        --is_odd[pointer_to_odd[s]];
                        pointer_to_odd[s] = i;
                    }
                    if (maxbackval[s] >= 1)
                    {
                        if (maxbackval[s] >= 1.66 || pos[s] != 3)
                        {
                            pos[s] = 1;
                        }
                    }
                    coboundary += x_a;
                }

            }

            void add_node_to_handle(const int v)
            /*Egy csúcs hozzáadásakor:
                - Frissíteni kell a maxback tömböt
                - Frissíteni Kell a listákat
                - Frissíteni kell a coboundary értéket (Kivonni ami bele ment (Ez pont a maxback érték, És hozzáadni amí kimegy))
            */
            {

                --is_odd[pointer_to_odd[v]];    // Ha v tett valakit odd csúccsá az ne legyen az
                pos[v] = 0; ++size;             // Mostantól benne van, eggyel nőtt a nyél
                coboundary -= maxbackval[v];    // Pont annyival járult hozzá a vágás értékéhez amennyivel látjuk
                for (ListDigraph::InArcIt a(O.G, O.V[v]); a != INVALID; ++a)    // Nézzük végig az éleit
                {
                    double x_a = 0;
                    if (!O.IsCol[a] || (x_a = O.coreLP.primal(O.Col[a])) == 0)
                    {
                        continue; //Ha $a$ nem oszlop vagy x_a = 0 menjünk tovább
                    }
                    const int s = O.Label[O.G.source(static_cast<ListDigraphBase::Arc>(a))];

                    maxbackval[s] += x_a;
                    if (pos[s] == 0){continue;} // Ha benne van biztos nem kerül ki
                    coboundary += x_a;          // Ha itt járunk kívül van hozzáadjuk az él értékét a vágáshoz
                    if (x_a == 1)               // Ha az él egy legyen extrém csúcs, de ezt később felülírhatjuk
                    {
                        pointer_to_outside[v] = s;      // v-nek mostantól s a kinti 1-el látott szomszédja
                        ++is_odd[v];
                        --is_odd[pointer_to_odd[s]];
                        pointer_to_odd[s] = v;
                        pos[s] = 3;
                    }
                    if ((maxbackval[s] >= 1 && pos[s] != 3) || maxbackval[s] >= 1.66)
                        // Ha nagyon látjuk, vagy eléggé látjuk és nem extrém legyen subto1
                    {
                        pos[s] = 1;
                    }
                }

                for (ListDigraph::OutArcIt a(O.G, O.V[v]); a != INVALID; ++a)    // Nézzük végig a másik irányú éleit
                {
                    double x_a = 0;
                    if (!O.IsCol[a] || (x_a = O.coreLP.primal(O.Col[a])) == 0)
                    {
                        continue; //Ha $a$ nem oszlop vagy x_a = 0 menjünk tovább
                    }
                    const int t = O.Label[O.G.target(static_cast<ListDigraphBase::Arc>(a))];

                    maxbackval[t] += x_a;
                    if (pos[t] == 0){continue;} // Ha benne van biztos nem kerül ki
                    coboundary += x_a;          // Ha itt járunk kívül van hozzáadjuk az él értékét a vágáshoz
                    if (x_a == 1)               // Ha az él egy legyen extrém csúcs, de ezt később felülírhatjuk
                    {
                        pos[t] = 3;
                        pointer_to_outside[v] = t;      // v-nek mostantól t a kinti 1-el látott szomszédja
                        ++is_odd[v];
                        --is_odd[pointer_to_odd[t]];  //Egy csúcs mindig max egy csúcsot tegyen oddá
                        pointer_to_odd[t] = v;
                    }
                    if ((maxbackval[t] >= 1 && pos[t] != 3) || maxbackval[t] >= 1.66)
                        // Ha nagyon látjuk, vagy eléggé látjuk és nem extrém, legyen subto1
                    {
                        pos[t] = 1;
                    }
                }
            };

            bool produce()
            {

                 do
                     {
                     int new_v{O.n};
                     vector<int> Candidates(4, O.n);
                     //Csinálunk egy vektort ami eltárolja minden tipusból a legjobb értékű indexét

                     for (int i = 0; i < O.n; ++i) // Nézzük végig az összes csúcsot
                     {
                         if (pos[i]) //Ha nem belső csúcs
                         {
                             if (maxbackval[Candidates[pos[i]]] <= maxbackval[i])
                                 // Nézzük meg, hogy az ilyen tipus legjobbnál jobb-e
                             {
                                 Candidates[pos[i]] = i;
                             }
                         }
                     }


                     int j = 1;
                     //Melyiket rakjuk be, Ha van subto1 biztos azt, ha van infto1 azt, egyébként egy extremof1
                     while (Candidates[j] == O.n)
                     {
                         ++j;
                     }
                     new_v = Candidates[j];
                     this->add_node_to_handle(new_v);
                     }
                while(std::abs(coboundary / 2 - std::round(coboundary)) <= 0.2 && size < O.n/2);
                //Legalább 0.2-re van a legközelebbi páros egésztől visszadjuk, ha nem növeljünk.
                //Mindig tudunk növelni mert a halmazok az összes csúcsot tartalmazzák

                if (size < O.n/2) // A nyél legfeljebb n/2 nagy legyen.
                {
                    return true;
                }
                return false;
            }
        };

        void lookforteeth(maxback_Handle_growing& H)
        {
            int num_of_found_teeth = 0;

            int target_num_of_teeth = std::floor(H.coboundary / 2.0) * 2.0 + 1.0;
            //A keresett fogak száma 2k+1 ha 2k < H.coboundary < 2k+2
            std::vector<boost::dynamic_bitset<>> teeth(target_num_of_teeth, boost::dynamic_bitset<>(n));
            //A Lendő fogak vektora.

            list<int> starting_points;
            // Ezekeből a csúcsokból fogunk elkezdeni fogat növeszteni

            vector<int> status(n, 0); // 0: authorized, 1 odd, 2 forbidden, 3 extrem csúcs
            for (int i = 0; i < n; ++i)
            {
                if (H.is_odd[i] > 0)
                {
                    starting_points.push_back(i);
                    status[i] = 1;
                }else
                {
                    if (H.pos[i] == 3)
                    {
                        status[i] = 3;
                    }else
                    {
                        status[i] = 0;
                    }
                }
            }

            while ((num_of_found_teeth < target_num_of_teeth) && (starting_points.size() != 0))
                // Legyenek a fogak a kifele 1 él csúcsai
            {
                int t_1 = starting_points.back();
                starting_points.pop_back();
                int t_2 = H.pointer_to_outside[t_1];
                teeth[num_of_found_teeth].set(t_1);
                teeth[num_of_found_teeth].set(t_2);
                status[t_1] = status[t_2] = 2;
                ++num_of_found_teeth;
            }

            //Ha nincs elég ilyen csúcs vegyünk be maxback növesztett fogakat.
            int num_of_nodes_available = H.size - num_of_found_teeth; //Vigyázzunk ne fogyanak el a H-beli csúcsok

            while (num_of_found_teeth < target_num_of_teeth &&
                (target_num_of_teeth - num_of_found_teeth) < num_of_nodes_available)
                //Ha kell még fog és van rá csúcs
            {
                // Keressünk olyan nem tiltott csúcsot aminek minimális a maxback értéke H-ban
                int i = 0;
                for (; i < n; ++i)
                {
                    if (H.pos[i] == 0 && status[i] == 0)
                    {
                        break;
                    }
                }

                for (int j = i+1; j < n; ++j)
                {
                    if (H.maxbackval[i] > H.maxbackval[j] && H.pos[j] == 0 && status[j] == 0)
                    {
                        i = j;
                    }
                }


                // Kezdjünk el ebből max_back növelni
                vector<int> tooth_max_back(n, 0);
                double coboundary = 0;

                for (ListDigraph::OutArcIt a(G, V[i]); a != INVALID; ++a)
                {
                    double x_a = 0;
                    if (!IsCol[a] || (x_a = coreLP.primal(Col[a])) == 0)
                    {
                        continue; //Ha $a$ nem oszlop vagy x_a = 0 marad 0
                    }
                    const int t = Label[G.target(static_cast<ListDigraphBase::Arc>(a))];
                    tooth_max_back[t] += x_a;
                    coboundary += x_a;
                }

                for (ListDigraph::OutArcIt a(G, V[i]); a != INVALID; ++a)
                {
                    double x_a = 0;
                    if (!IsCol[a] || (x_a = coreLP.primal(Col[a])) == 0)
                    {
                        continue; //Ha $a$ nem oszlop vagy x_a = 0 marad 0
                    }
                    const int t = Label[G.target(static_cast<ListDigraphBase::Arc>(a))];
                    tooth_max_back[t] += x_a;
                    coboundary += x_a;
                }
            }
        }


        bool Comb_separation(const vector<int>& IN,const vector<int>& OUT) //Returns true if comb ineqs are found
        {
            bool found = false;
            int root = 0;
            vector<vector<int>> potential_handles{}; // 1esekből álló út egyik, másik vége párok vannak benne
            while (root < n)
            {
                if (IN[root] == -1 || OUT[root] == -1)
                {
                    ++root;
                    continue;
                }

                int node = OUT[root];

                while (OUT[node] != -1)
                {
                    node = OUT[node];
                };
                potential_handles.emplace_back(root, node);
                potential_handles.emplace_back(node, root);
            }

            while (!potential_handles.empty())
            {
                maxback_Handle_growing H_factory = potential_handles.back();
                potential_handles.pop_back();

                while (H_factory.produce())
                {
                    //lookforteeth(H_factory);
                }
            }

            return found;
        }
#pragma endregion

#pragma region Branching
        void Branch(size_t i, node_BnCnPtree &X)
        {
            ListDigraph::Node t = G.target(Columns[i].first);
            if constexpr (BRANCHING == 1)
            {
                logger.log("Branching on InArcs of node ", Label[t]);
                vector<size_t> Candidate_arcs;
                vector<size_t> Unused_tarcs;
                for (size_t j = 0; j < Columns.size(); j++)
                {
                    if (X.Edges_setto_0[j] || X.Edges_setto_1[j])
                    {
                        continue;
                    }
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
                if (Candidate_arcs.size() != 0)
                {
                    nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver); // Where all the arcs with high value are set to 0.
                }else
                {
                    Candidate_arcs = Unused_tarcs;
                }
                for (size_t k : Unused_tarcs)
                {
                    new0s.set(k);
                }

                for (size_t j : Candidate_arcs) // Where a Single high value arc is set to 1, and all others to 0.
                {
                    new1s.set(j);
                    new0s.reset(j);
                    nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver);
                    new1s.reset(j);
                    new0s.set(j);
                }
            }
            else if constexpr (BRANCHING == 2)
            {
                boost::dynamic_bitset<> new0s = X.Edges_setto_0;
                boost::dynamic_bitset<> new1s = X.Edges_setto_1;
                new1s.set(i);
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver);
                new1s.reset(i);
                new0s.set(i);
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver);
            }
            else if constexpr (BRANCHING == 3)
            {
                int j = i;
                for (; i < Columns.size(); i++)
                {
                    if (std::abs(coreLP.primal(Columns[i].second) - 0.5) < 1)
                    {
                        break;
                    }
                }

                if (i == Columns.size())
                {
                    i = j;
                }
                boost::dynamic_bitset<> new0s = X.Edges_setto_0;
                boost::dynamic_bitset<> new1s = X.Edges_setto_1;
                new1s.set(i);
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver);
                new1s.reset(i);
                new0s.set(i);
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver);
            }
            else if constexpr (BRANCHING == 4)
            {
                vector<size_t> S_1;
                vector<size_t> S_2;
                int S_1_used = 0;
                int S_1_unused = 0;
                int S_2_used = 0;
                int S_2_unused = 0;
                for (size_t j = 0; j < Columns.size(); j++)
                {
                    if (X.Edges_setto_0[j] || X.Edges_setto_1[j])
                    {
                        continue;
                    }
                    if (G.target(Columns[j].first) == t)
                    {
                        if (coreLP.primal(Columns[j].second) > 0.3)
                        {
                            if (S_1_used < S_2_used)
                            {
                                S_1.push_back(j);
                                S_1_used++;
                            }
                            else
                            {
                                S_2.push_back(j);
                                S_2_used++;
                            }
                        }
                        else
                        {
                            if (S_1_unused < S_2_unused)
                            {
                                S_1.push_back(j);
                                S_1_unused++;
                            }
                            else
                            {
                                S_2.push_back(j);
                                S_2_unused++;
                            }
                        }
                    }
                }
                boost::dynamic_bitset<> new0s = X.Edges_setto_0;
                boost::dynamic_bitset<> new1s = X.Edges_setto_1;
                new0s.resize(Columns.size());
                new1s.resize(Columns.size());

                for (size_t k : S_1)
                {
                    new0s.set(k);
                }
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver); // Where all the arcs in S_2 set are set to 0.
                for (size_t k : S_1)
                {
                    new0s.reset(k);
                }
                for (size_t k : S_2)
                {
                    new0s.set(k);
                }
                nodes.emplace_front(*this, X.ID, new0s, new1s, X.LB, ++idgiver); // Where all the arcs in S_2 set are set to 0.
            }
        }
#pragma endregion

#pragma region Processing_Node
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

            // Column generation step
            bool Not_Optimally_solved = true;
            while (Not_Optimally_solved){
                coreLP.solve();
                Not_Optimally_solved = Pricing();
            }

            // Now we can Compute a Lowerbound
            coreLP.solve();
            X.LB = coreLP.primal();
            if (!(X.LB <= Best_Val)) // Cutting step
            {
                return;
            }

            // sec_separation step
            if (sec_separation())
            {
                nodes.emplace_front(*this, X.ID, X.Edges_setto_0, X.Edges_setto_1, X.LB, ++idgiver);
                return;
            }

            //Making the sol intiger
            // If the solution is feasablie constraintwise we start looking for combs then fixing non intiger arcs as integers
            size_t i = 0;
            vector<int> IN(n, -1), OUT(n,-1);
            size_t j = Columns.size();
            while (i < Columns.size())
            {
                if (std::abs(coreLP.primal(Columns[i].second) - std::round(coreLP.primal(Columns[i].second))) < 1e-7)
                {
                    if (std::abs(coreLP.primal(Columns[i].second) - 1) < 1e-7) // Keressünk egyesekből álló utakat
                    {
                        int s = Label[G.source(Columns[i].first)];
                        int t = Label[G.target(Columns[i].first)];
                        IN[t] = s;
                        OUT[s] = t;
                    }
                }
                else
                {
                    j = i;
                }
                i++;
            }
            if (j < Columns.size())
            {
                Branch(j, X);
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
            const double timelimt = 10,
            const double _Upper_bound = std::numeric_limits<double>::max(),
            const int sep = 1000,
            const int col = 10,
            const int exploring_steps = 10,
            const size_t _initarcs = 15) : G(_G),
                                           Label(_Label),
                                           weight(_weight),
                                           IsCol(_G, false),
                                           Lookup(G),
                                           Upper_bound(_Upper_bound),
                                           Col(_G,INVALID),
                                           sepequality_add_number(sep),
                                           col_add_number(col),
                                           exploring_steps(exploring_steps),
                                           initarcs(_initarcs),
                                           timeout(timelimt)
        {
            n = countNodes(_G);
            if (initarcs > static_cast<size_t>(n-1)){initarcs = static_cast<size_t>(n-1);}
            V.resize(n);
            for (ListDigraph::NodeIt i(G); i != INVALID; ++i)
            {
                V[_Label[i]] = static_cast<ListDigraphBase::Node>(i);
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
            cout << " Arcs left out: " << Arcs_notincluded.size() << endl;
            solved = 1;
            return Best_Val;
        }

#pragma endregion

#pragma region Tour_Found
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

#pragma region Query_Functions
        bool get_OPTsolved() const { return OPTsolved; }

        double OPTval() const
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
                    O = static_cast<ListDigraphBase::Node>(v);
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