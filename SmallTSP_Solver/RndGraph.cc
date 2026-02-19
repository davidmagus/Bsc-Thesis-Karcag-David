#pragma region inculdok_namespaces
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


#pragma region segdfuggvenyek
// min–max intervallumban véletlen int:
int random_int(int min, int max) {
    static std::random_device rd;   // seed
    static std::mt19937 gen(rd());  // mersenne twister
    std::uniform_int_distribution<> dist(min, max);
    return dist(gen);
}
struct Point {
    int x, y;
};
#pragma endregion


int main()
{
    // weight generation range: [a,b]
    const int a = -100;
    const int b = 100;



    std::string typeofgraph;
    std::cout << "Input the Graph type: DG (Digraph), TSP , G (Graph default), Ship_problem SP (int k: capacity, int n: number of cities, vector<vector<int>> C: c_ij == cost of a c_ij ticket, vector<vector<int>> D: d_ij == Max tickets that can be sold Btween ij): \n";
    std::cin >> typeofgraph;

#pragma region TSP
    if (typeofgraph == "TSP") {
        int n;
        double mratio;
        cout << "Number of nodes (n): ";
        cin >> n;
        cout << "Edge ratio (mratio): ";
        cin >> mratio;

        int target_m = static_cast<int>(n * mratio);
        if (target_m < n) {
            cout << "Warning: For connectivity, m must be at least n. Setting m = n.\n";
            target_m = n;
        }

        // 1. Pontok generálása
        vector<Point> coords(n);
        for (int i = 0; i < n; ++i) {
            coords[i].x = random_int(0, 1000);
            coords[i].y = random_int(0, 1000);
        }

        ListDigraph G;
        ListDigraph::ArcMap<int> weight(G);
        vector<ListDigraph::Node> nodes;

        for (int i = 0; i < n; ++i) {
            nodes.push_back(G.addNode());
        }

        // Távolság számító lambda függvény a tisztább kódért
        auto get_dist = [&](int u, int v) {
            double dx = coords[u].x - coords[v].x;
            double dy = coords[u].y - coords[v].y;
            return static_cast<int>(ceil(sqrt(dx*dx + dy*dy)));
        };

        // 2. Összefüggőség biztosítása: Egy kör létrehozása (0->1->2->...->n-1->0)
        int current_m = 0;
        for (int i = 0; i < n; ++i) {
            int next = (i + 1) % n;
            ListDigraph::Arc a = G.addArc(nodes[i], nodes[next]);
            weight[a] = get_dist(i, next);
            current_m++;
        }

        // 3. További élek generálása a cél m eléréséig
        int max_possible_arcs = n * (n - 1);
        if (target_m > max_possible_arcs) target_m = max_possible_arcs;

        while (current_m < target_m) {
            int u_idx = random_int(0, n - 1);
            int v_idx = random_int(0, n - 1);

            // Csak akkor adjuk hozzá, ha nem hurokél és még nem létezik
            if (u_idx != v_idx && findArc(G, nodes[u_idx], nodes[v_idx]) == INVALID) {
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
        #pragma endregion
        #pragma region egyeb

    } else if (typeofgraph == "DG") {
        cout << "Directed graph selected \n";
        int n = 5;
        int m = 10;
        cout << "n (5):";
        cin  >> n;
        cout << "m (10):";
        cin  >> m;

        ListDigraph G;
        vector<ListDigraph::Node> nodes;
        for(int j = 0; j < n;j++)
        {
            ListDigraph::Node v = G.addNode();
            nodes.push_back(v);
        }
        for(int j = 0; j < m;j++)
        {
            int u = random_int(0, n-1);
            int v = random_int(0, n-1);
            while (u == v)
            {
                v = random_int(0, n-1);
            }
            G.addArc(nodes[u],nodes[v]);
        }
        ListDigraph::ArcMap<int> weigth(G);
        for (ListDigraph::ArcIt e(G); e != INVALID ; ++e)
        {
            weigth[e] = random_int(a, b);
        }
        
        std::ofstream f("digraph.lgf");
        DigraphWriter<ListDigraph> A(G, f);
        A
            .arcMap("weight",weigth)
            .run();
    } else if (typeofgraph == "SP")
    {
    cout << "Shipping project\n";
    int k = 5;
    int n = 10;

    cout << "k: ";
    cin >> k;
    cout << "n: ";
    cin >> n;

    ofstream fout("ship.txt");   // létrehozza, ha nem létezik
    if (!fout) {
        cerr << "Error: ship.txt could not be created.\n";
        return 1;
    }

    // --- Random generátor ---
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(1, 100);   // pozitív egészek

    // --- Kiírás n és k ---
    fout << n << " " << k << "\n";

    // --- Mátrix C ---
    vector<vector<int>> C(n, vector<int>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = dist(gen);

    // C sorfolytonosan
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fout << C[i][j] << " ";
        fout << "\n";
    }

    // --- Mátrix D ---
    vector<vector<int>> D(n, vector<int>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            D[i][j] = dist(gen);
    int m = 6*n;
    int __k = 0;
    // D sorfolytonosan
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            if ( i  < j && __k < m)
            {
            fout << D[i][j] << " ";
            __k++;
            }else{
                fout << 0 << " ";
            }
        fout << "\n";
    }

    fout.close();
    cout << "ship.txt generated successfully.\n";
    }

    else {
        // Default: Graph
        int n = 5;
        int m = 10;
        cout << "n (5):";
        cin  >> n;
        cout << "m (10):";
        cin  >> m;

        ListGraph G;
        vector<ListGraph::Node> nodes;
        for(int j = 0; j < n;j++)
        {
            ListGraph::Node v = G.addNode();
            nodes.push_back(v);
        }
        for(int j = 1; j < m;j++)
        {
            int u = random_int(0, n-1);
            int v = random_int(0, n-1);
            while (u == v)
            {
                v = random_int(0, n-1);
            }
            G.addEdge(nodes[u],nodes[v]);
        }
        ListGraph::EdgeMap<int> weight(G);
        for(ListGraph::EdgeIt e(G); e != INVALID; ++e){
            weight[e] = random_int(a, b);
        }
        std::ofstream f("graph.lgf");
        GraphWriter<ListGraph> writer(G, f);
        writer
            .edgeMap("weight",weight)
            .run();

    }
    #pragma endregion
}