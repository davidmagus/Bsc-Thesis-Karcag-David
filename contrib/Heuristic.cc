#include "Heuristic.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/maps.h>
using namespace std;
using namespace lemon;

namespace Heuristic
{
#pragma region Logging Tools
    struct Logging
    {
        std::ofstream logFile;
        int counter;
        Logging() : logFile("Heu.log") { counter = 1500; }

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

#pragma region Nearest Neighbour
    Nearest_Neighbour::Nearest_Neighbour(
        const ListDigraph &_G,
        const ListDigraph::NodeMap<int> &_Label,
        const ListDigraph::ArcMap<double> &_weight,
        const int root) : G(_G), Label(_Label), weight(_weight)
    {
        size_t n = countNodes(_G);
        Length = 0;
        ListDigraph::NodeMap<int> seen(G, 0);
        ListDigraph::Node firstnode;
        for (ListDigraph::NodeIt a(G); a != INVALID; ++a)
        {
            if (Label[a] == root)
            {
                firstnode = a;
                break;
            }
        }

        ListDigraph::Node currentnode = firstnode;
        seen[currentnode] = 1;
        while (Route.size() < n - 1)
        {
            ListDigraph::OutArcIt a(G, currentnode);
            while (seen[(G).target(a)])
            {
                ++a;
            }
            ListDigraph::Arc Bestout = a;
            ++a;
            while (a != INVALID)
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

        for (ListDigraph::OutArcIt a(G, G.target(Route.back())); a != INVALID; ++a)
        {
            if (G.target(a) == firstnode)
            {
                lastarc = a;
                break;
            }
        }
        Route.push_back(lastarc);
        Length += weight[lastarc];
    };
#pragma endregion

#pragma region Repetitive Nearest Neighbour
    Repetitive_Nearest_Neighbour::Repetitive_Nearest_Neighbour(
        const ListDigraph &_G,
        const ListDigraph::NodeMap<int> &_Label,
        const ListDigraph::ArcMap<double> &_weight) : G(_G), Label(_Label), weight(_weight)
    {
        n = countNodes(G);
        Length = std::numeric_limits<double>::max();
        for (int i = 0; i < n - 1; i++)
        {
            Heuristic::Nearest_Neighbour NN{G, Label, weight, i};
            if (NN.Length < Length)
            {
                Length = NN.Length;
                Route = NN.Route;
            }
        }
    }
#pragma endregion

#pragma region Greedy algorithm
    Greedy::Greedy(
        const ListDigraph &_G,
        const ListDigraph::NodeMap<int> &_Label,
        const ListDigraph::ArcMap<double> &_weight) : G(_G), Label(_Label), weight(_weight)
    {
        Logging logger;
        Length = 0;
        for (ListDigraph::ArcIt a(G); a != INVALID; ++a)
        {
            Orderedarcs.push_back(a);
        }
        std::sort(Orderedarcs.begin(), Orderedarcs.end(),
                  [&_weight](const ListDigraph::Arc &a, const ListDigraph::Arc &b)
                  {
                      return _weight[a] < _weight[b];
                  });

        ListDigraph::NodeMap<bool> IN(G, false);
        ListDigraph::NodeMap<bool> OUT(G, false);

        ListDigraph::NodeMap<ListDigraph::Arc> Out(G);
        ListDigraph::NodeMap<ListDigraph::Arc> In(G);

        ListDigraph::NodeMap<ListDigraph::Node> end(G);

        for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
        {
            end[v] = v;
        }

        size_t i = 0;
        int cnt = 0;
        while (i < Orderedarcs.size())
        {
            ListDigraph::Arc a = Orderedarcs[i];
            ListDigraph::Node s = G.source(a);
            ListDigraph::Node t = G.target(a);

            if (!IN[t] && !OUT[s] && (end[s] != t || cnt == countNodes(G) - 1))
            {
                cnt++;
                end[t] = end[s];
                ListDigraph::Node next = t;
                while (OUT[next])
                {
                    end[G.target(Out[next])] = end[next];
                    next = G.target(Out[next]);
                    if (next == t)
                    {
                        break;
                    }
                }
                IN[t] = OUT[s] = true;
                In[t] = Out[s] = a;
            }
            ++i;
        }

        ListDigraph::NodeIt v(G);
        ListDigraph::Node r = v;
        ListDigraph::Arc Last = Out[r];
        Route.push_back(Last);
        Length += weight[Last];
        int counter = countNodes(G) * 10;
        do
        {
            Last = Out[G.target(Last)];
            Route.push_back(Last);
            Length += weight[Last];
            counter--;
        } while (G.target(Last) != r && counter);
        if (!counter)
        {
            cout << "\n no good \n";
        }
    };
#pragma endregion

#pragma region Insertion Methods

    Insertion_Methods::Insertion_Methods(
        const ListDigraph &_G,
        const ListDigraph::NodeMap<int> &_Label,
        const ListDigraph::ArcMap<double> &_weight)
        : n(countNodes(_G)), G(_G), Label(_Label), weight(_weight), Length(0)
    {
        for (ListDigraph::NodeIt a(G); a != INVALID; ++a)
        {
            Unseen.push_front(a);
        }
    }

    void Max_insert::Insert(lemon::ArcLookUp<ListDigraph> &lookUp)
    {
        if (Route.empty())
        {
            double Maxdist = -1; // minnél átálítani
            ListDigraph::Node s;
            ListDigraph::Node t;
            for (ListDigraph::Node u : Unseen)
            {
                for (ListDigraph::Node v : Unseen)
                {
                    if (v != u)
                    {
                        double uvinsertweight = weight[lookUp(u, v)] + weight[lookUp(v, u)];
                        if (Maxdist < uvinsertweight)
                        {
                            Maxdist = uvinsertweight;
                            s = u;
                            t = v;
                        }
                    }
                }
            }
            Length = Maxdist;
            Route.push_back(lookUp(s, t));
            Route.push_back(lookUp(t, s));
            Unseen.remove(s);
            Unseen.remove(t);
        }
        else
        {
            double Maxdist = -1;
            int maxplace = 0;
            ListDigraph::Node s;

            for (ListDigraph::Node u : Unseen)
            {
                int minplace = 0;
                int mininsertweight = Length - weight[Route[0]] + weight[lookUp(G.source(Route[0]), u)] + weight[lookUp(u, G.target(Route[0]))];
                for (size_t i = 1; i < Route.size(); i++)
                {
                    int insertweight = Length - weight[Route[i]] + weight[lookUp(G.source(Route[i]), u)] + weight[lookUp(u, G.target(Route[i]))];
                    if (insertweight < mininsertweight)
                    {
                        mininsertweight = insertweight;
                        minplace = i;
                    }
                }
                if (Maxdist < mininsertweight)
                {
                    s = u;
                    maxplace = minplace;
                    Maxdist = mininsertweight;
                }
            }

            Length = Maxdist;
            ListDigraph::Node t = G.target(Route[maxplace]);
            Route[maxplace] = lookUp(G.source(Route[maxplace]), s);
            auto it = Route.begin() + maxplace + 1;
            Route.insert(it, lookUp(s, t));
            Unseen.remove(s);
        }
    }

    double Max_insert::Run()
    {
        lemon::ArcLookUp<ListDigraph> lookUp(G);
        while (!Unseen.empty())
        {
            Insert(lookUp);
        }
        return Length;
    }

    void Min_insert::Insert(lemon::ArcLookUp<ListDigraph> &lookUp)
    {
        if (Route.empty())
        {
            double mindist = std::numeric_limits<double>::max();
            ListDigraph::Node s;
            ListDigraph::Node t;
            for (ListDigraph::Node u : Unseen)
            {
                for (ListDigraph::Node v : Unseen)
                {
                    if (v != u)
                    {
                        double uvinsertweight = weight[lookUp(u, v)] + weight[lookUp(v, u)];
                        if (mindist > uvinsertweight)
                        {
                            mindist = uvinsertweight;
                            s = u;
                            t = v;
                        }
                    }
                }
            }
            Length = mindist;
            Route.push_back(lookUp(s, t));
            Route.push_back(lookUp(t, s));
            Unseen.remove(s);
            Unseen.remove(t);
        }
        else
        {
            double mindist = 10000000;
            int maxplace = 0;
            ListDigraph::Node s;
            for (ListDigraph::Node u : Unseen)
            {
                int minplace = 0;
                int mininsertweight = Length - weight[Route[0]] + weight[lookUp(G.source(Route[0]), u)] + weight[lookUp(u, G.target(Route[0]))];
                for (size_t i = 1; i < Route.size(); i++)
                {
                    int insertweight = Length - weight[Route[i]] + weight[lookUp(G.source(Route[i]), u)] + weight[lookUp(u, G.target(Route[i]))];
                    if (insertweight < mininsertweight)
                    {
                        mininsertweight = insertweight;
                        minplace = i;
                    }
                }
                if (mindist > mininsertweight)
                {
                    s = u;
                    maxplace = minplace;
                    mindist = mininsertweight;
                }
            }
            Length = mindist;
            ListDigraph::Node t = G.target(Route[maxplace]);
            Route[maxplace] = lookUp(G.source(Route[maxplace]), s);
            auto it = Route.begin() + maxplace + 1;
            Route.insert(it, lookUp(s, t));
            Unseen.remove(s);
        }
    }

    double Min_insert::Run()
    {
        lemon::ArcLookUp<ListDigraph> lookUp(G);
        while (!Unseen.empty())
        {
            Insert(lookUp);
        }
        return Length;
    }

    void Rand_insert::Insert(lemon::ArcLookUp<ListDigraph> &lookUp)
    {
        if (Route.empty())
        {
            ListDigraph::Node u = Unseen.front();
            Unseen.pop_front();
            ListDigraph::Node v = Unseen.front();
            Unseen.pop_front();
            ListDigraph::Arc uv = lookUp(u, v);
            ListDigraph::Arc vu = lookUp(v, u);
            Route.push_back(uv);
            Route.push_back(vu);
            Length = weight[uv] + weight[vu];
        }
        else
        {
            ListDigraph::Node u = Unseen.front();
            Unseen.pop_front();

            int minplace = 0;
            int mininsertweight = Length - weight[Route[0]] + weight[lookUp(G.source(Route[0]), u)] + weight[lookUp(u, G.target(Route[0]))];
            for (size_t i = 1; i < Route.size(); i++)
            {
                int insertweight = Length - weight[Route[i]] + weight[lookUp(G.source(Route[i]), u)] + weight[lookUp(u, G.target(Route[i]))];
                if (insertweight < mininsertweight)
                {
                    mininsertweight = insertweight;
                    minplace = i;
                }
            }
            Length = mininsertweight;
            ListDigraph::Node t = G.target(Route[minplace]);
            Route[minplace] = lookUp(G.source(Route[minplace]), u);
            auto it = Route.begin() + minplace + 1;
            Route.insert(it, lookUp(u, t));
        }
    }

    double Rand_insert::Run()
    {
        lemon::ArcLookUp<ListDigraph> lookUp(G);
        while (!Unseen.empty())
        {
            Insert(lookUp);
        }
        return Length;
    }
#pragma endregion

}