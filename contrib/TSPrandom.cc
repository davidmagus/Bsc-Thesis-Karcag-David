#pragma region inculdok_namespaces
#include "TSP_generator.h"
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>
#include <random>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
using namespace std;
using namespace lemon;
#pragma endregion

namespace RandTSP
{
#pragma region Helper Tools
    static std::mt19937 global_gen;
    void init_random(int user_seed)
    {
        global_gen.seed(static_cast<unsigned int>(user_seed) * 1234567);
    }
    int get_next_random(int min, int max)
    {
        std::uniform_int_distribution<> dist(min, max);
        return dist(global_gen);
    }
    struct Point
    {
        int x, y;
    };
#pragma endregion

#pragma region Not Complete TSP
    void Make_non_completeTSP(const int n, const double mratio, const int seed)
    {
        int user_seed = n + seed;
        init_random(user_seed);
        int target_m = static_cast<int>(n * mratio);
        if (target_m < n)
        {
            cout << "Warning: For connectivity, m must be at least n. Setting m = n.\n";
            target_m = n;
        }

        // 1. Pontok generálása
        vector<Point> coords(n);
        for (int i = 0; i < n; ++i)
        {
            coords[i].x = get_next_random(0, 1000);
            coords[i].y = get_next_random(0, 1000);
        }

        ListDigraph G;
        ListDigraph::ArcMap<int> weight(G);
        vector<ListDigraph::Node> nodes;

        for (int i = 0; i < n; ++i)
        {
            nodes.push_back(G.addNode());
        }

        // Távolság számító lambda függvény a tisztább kódért
        auto get_dist = [&](int u, int v)
        {
            double dx = coords[u].x - coords[v].x;
            double dy = coords[u].y - coords[v].y;
            return static_cast<int>(ceil(sqrt(dx * dx + dy * dy)));
        };

        // 2. Összefüggőség biztosítása: Egy kör létrehozása (0->1->2->...->n-1->0)
        int current_m = 0;
        for (int i = 0; i < n; ++i)
        {
            int next = (i + 1) % n;
            ListDigraph::Arc a = G.addArc(nodes[i], nodes[next]);
            weight[a] = get_dist(i, next);
            current_m++;
        }

        // 3. További élek generálása a cél m eléréséig
        int max_possible_arcs = n * (n - 1);
        if (target_m > max_possible_arcs)
            target_m = max_possible_arcs;

        while (current_m < target_m)
        {
            int u_idx = get_next_random(0, n - 1);
            int v_idx = get_next_random(0, n - 1);

            // Csak akkor adjuk hozzá, ha nem hurokél és még nem létezik
            if (u_idx != v_idx && findArc(G, nodes[u_idx], nodes[v_idx]) == INVALID)
            {
                ListDigraph::Arc a = G.addArc(nodes[u_idx], nodes[v_idx]);
                weight[a] = get_dist(u_idx, v_idx);
                current_m++;
            }
        }

        // 4. Kiírás
        std::ofstream f("digraph_tsp.lgf");
        DigraphWriter<ListDigraph>(G, f)
            .arcMap("weight", weight)
            .run();

        cout << "Connected TSP graph generated successfully (" << current_m << " arcs).\n";
    }
#pragma endregion

#pragma region TSP
    void Make_completeTSP(const int n, const int seed)
    {
        int user_seed = n + seed;
        init_random(user_seed);
        int current_m = 0;
        vector<Point> coords(n);
        for (int i = 0; i < n; ++i)
        {
            coords[i].x = get_next_random(0, 100);
            coords[i].y = get_next_random(0, 100);
        }

        ListDigraph G;
        ListDigraph::ArcMap<int> weight(G);
        vector<ListDigraph::Node> nodes;

        for (int i = 0; i < n; ++i)
        {
            nodes.push_back(G.addNode());
        }

        // Távolság számító lambda függvény a tisztább kódért
        auto get_dist = [&](int u, int v)
        {
            double dx = coords[u].x - coords[v].x;
            double dy = coords[u].y - coords[v].y;
            return static_cast<int>(ceil(sqrt(dx * dx + dy * dy)));
        };

        for (size_t u_idx = 0; u_idx < nodes.size(); ++u_idx)
        {
            for (size_t v_idx = 0; v_idx < nodes.size(); ++v_idx)
            {
                if (u_idx != v_idx)
                {
                    ListDigraph::Arc a = G.addArc(nodes[u_idx], nodes[v_idx]);
                    weight[a] = get_dist(u_idx, v_idx);
                    current_m++;
                }
            }
        }

        // 4. Kiírás
        std::ofstream f("digraph_tsp.lgf");
        DigraphWriter<ListDigraph>(G, f)
            .arcMap("weight", weight)
            .run();

        cout << "Connected TSP graph generated successfully (" << current_m << " arcs).\n";
    }
#pragma endregion
}