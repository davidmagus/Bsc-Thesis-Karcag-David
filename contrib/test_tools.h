#ifndef Test_tool
#define Test_tool
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include "Heuristic.h"
#include "Heldkarp.h"
#include "TSP_generator.h"
#include "BnC.h"
#include <string>
using namespace std;
using namespace lemon;

namespace test_tools
{

    struct entry
    {
        int n;

        // Nearest Neighbour details
        double NN_val;
        double NN_time;
        double NN_percentage;

        // Repetitive Nearest Neighbour details
        double RNN_val;
        double RNN_time;
        double RNN_percentage;

        // Greedy algorithm details
        double Grdy_val;
        double Grdy_time;
        double Grdy_percentage;

        // Mininsert algorithm details
        double Mininsert_val;
        double Mininsert_time;
        double Mininsert_percentage;

        // Maxinsert algorithm details
        double Maxinsert_val;
        double Maxinsert_time;
        double Maxinsert_percentage;

        // Randinsert algorithm details
        double Randinsert_val;
        double Randinsert_time;
        double Randinsert_percentage;

        // Branch and Count details
        bool BnC_solved;
        double BnC_val;
        double BnC_Time;

        // Held-Karp details
        bool HK_solved;
        double HK_val;
        double HK_Time;

        entry(
            const int n,

            // Closest Neighbour details
            const double NN_val = 0,
            double NN_time = 0,

            // Repetitive Nearest Neighbour details
            double RNN_val = 0,
            double RNN_time = 0,

            // Greedy algorithm details
            double Grdy_val = 0,
            double Grdy_time = 0,

            // Maxinsert algorithm details
            double Maxinsert_val = 0,
            double Maxinsert_time = 0,

            // Mininsert algorithm details
            double Mininsert_val = 0,
            double Mininsert_time = 0,

            // Randinsert algorithm details
            double Randinsert_val = 0,
            double Randinsert_time = 0,

            // Branch and Count details
            bool BnC_solved = false,
            double BnC_val = 0,
            double BnC_Time = 0,

            // Held-Karp details
            bool HK_solved = false,
            double HK_val = 0,
            double HK_Time = 0) : n(n), NN_val(NN_val), NN_time(NN_time), Grdy_val(Grdy_val), Grdy_time(Grdy_time), BnC_solved(BnC_solved), BnC_val(BnC_val), BnC_Time(BnC_Time), HK_solved(HK_solved), HK_val(HK_val), HK_Time(HK_Time)
        {
            finalize();
        }
        void finalize()
        {
            if (BnC_val && BnC_solved)
            {
                NN_percentage = NN_val / BnC_val;
                RNN_percentage = RNN_val / BnC_val;
                Grdy_percentage = Grdy_val / BnC_val;
                Maxinsert_percentage = Maxinsert_val / BnC_val;
                Mininsert_percentage = Mininsert_val / BnC_val;
                Randinsert_percentage = Randinsert_val / BnC_val;
            }
            else
            {
                NN_percentage = 0;
                Grdy_percentage = 0;
                RNN_percentage = 0;
                Maxinsert_percentage = 0;
                Mininsert_percentage = 0;
                Randinsert_percentage = 0;
            }
        };
    };

    std::vector<entry> Results;

    int one_run(int num, int seed = 0, vector<string> what_to_do = {"BnC", "Heu", "HK"})
    {
        Timer timer;
#pragma region Generating Examples
        RandTSP::Make_completeTSP(num, seed);
        // cout << "Graph Made in " << timer.realTime() << endl;
        timer.reset();

#pragma endregion

#pragma region Reading Input
        ListDigraph G;
        ListDigraph::NodeMap<int> Label(G);
        ListDigraph::ArcMap<double> weight(G);
        string filename = "digraph_tsp.lgf";

        // cout << "Reading Input..." << "\n";

        ifstream f(filename);
        if (!f)
        {
            cerr << "Hiba: A " << filename << " nem talalhato!" << "\n";
            return 1;
        }

        try
        {
            DigraphReader<ListDigraph> reader(G, f);
            reader
                .nodeMap("label", Label)
                .arcMap("weight", weight)
                .run();
        }
        catch (const Exception &e)
        { // LEMON saját kivétel osztálya
            cerr << "LEMON Hiba: " << e.what() << "\n";
            return 1;
        }
#pragma endregion
        Results.emplace_back(num);
        {

            int Heu = count(what_to_do.begin(), what_to_do.end(), "Heu");
            if (Heu)
            {
#pragma region Nearest Neighbour
                {
                    timer.restart();
                    cout << "\nNearest Neighbour...\n";
                    Heuristic::Nearest_Neighbour NN{G, Label, weight};
                    // cout << Label[(G.source(NN.Route[0]))];
                    Results.back().NN_time = timer.realTime();
                    Results.back().NN_val = NN.Length;
                    std::cout << "Tour value: " << Results.back().NN_val << endl;
                    // for (size_t i = 0; i < NN.Route.size(); i++)
                    // {
                    //     cout << Label[(G.target(NN.Route[i]))];
                    // }
                    // cout << " Tour Length: " << NN.Length << endl;
                    // cout << "Time: " << timer.realTime() << endl;
                }
#pragma endregion

#pragma region Repetitive Nearest Neighbour
                {
                    timer.restart();
                    cout << "\nRepetitive Nearest Neighbour...\n";
                    Heuristic::Repetitive_Nearest_Neighbour RNN{G, Label, weight};
                    // cout << Label[(G.source(NN.Route[0]))];
                    Results.back().RNN_time = timer.realTime();
                    Results.back().RNN_val = RNN.Length;
                    std::cout << "Tour value: " << Results.back().RNN_val << endl;
                    // for (size_t i = 0; i < NN.Route.size(); i++)
                    // {
                    //     cout << Label[(G.target(NN.Route[i]))];
                    // }
                    // cout << " Tour Length: " << NN.Length << endl;
                    // cout << "Time: " << timer.realTime() << endl;
                }
#pragma endregion

#pragma region Greedy
                {
                    timer.restart();
                    cout << "\nGreedy method...\n";
                    Heuristic::Greedy Grdy{G, Label, weight};
                    // cout << Label[(G.source(Grdy.Route[0]))];
                    Results.back().Grdy_time = timer.realTime();
                    Results.back().Grdy_val = Grdy.Length;
                    std::cout << "Tour value: " << Results.back().Grdy_val << endl;
                    // for (size_t i = 0; i < Grdy.Route.size(); i++)
                    // {
                    //     cout << Label[(G.target(Grdy.Route[i]))];
                    // }
                    // cout << " Tour Length: " << Grdy.Length << endl;
                    // cout << "Time: " << timer.realTime() << endl;
                }
#pragma endregion

#pragma region Inserting Methods
                {
                    timer.restart();
                    cout << "\nMaxinsert method...\n";
                    Heuristic::Max_insert Maxinsert{G, Label, weight};
                    Maxinsert.Run();
                    // cout << Label[(G.source(Grdy.Route[0]))];
                    Results.back().Maxinsert_time = timer.realTime();
                    Results.back().Maxinsert_val = Maxinsert.Length;
                    std::cout << "Tour value: " << Results.back().Maxinsert_val << endl;
                    // for (size_t i = 0; i < Grdy.Route.size(); i++)
                    // {
                    //     cout << Label[(G.target(Grdy.Route[i]))];
                    // }
                    // cout << " Tour Length: " << Grdy.Length << endl;
                    // cout << "Time: " << timer.realTime() << endl;
                }
                {
                    timer.restart();
                    cout << "\nMininsert method...\n";
                    Heuristic::Min_insert Minins{G, Label, weight};
                    Minins.Run();
                    // cout << Label[(G.source(Grdy.Route[0]))];
                    Results.back().Mininsert_time = timer.realTime();
                    Results.back().Mininsert_val = Minins.Length;
                    std::cout << "Tour value: " << Results.back().Mininsert_val << endl;
                    // for (size_t i = 0; i < Grdy.Route.size(); i++)
                    // {
                    //     cout << Label[(G.target(Grdy.Route[i]))];
                    // }
                    // cout << " Tour Length: " << Grdy.Length << endl;
                    // cout << "Time: " << timer.realTime() << endl;
                }
                {
                    timer.restart();
                    cout << "\nRandinsert method...\n";
                    Heuristic::Rand_insert Randinsert{G, Label, weight};
                    Randinsert.Run();
                    // cout << Label[(G.source(Grdy.Route[0]))];
                    Results.back().Randinsert_time = timer.realTime();
                    Results.back().Randinsert_val = Randinsert.Length;
                    std::cout << "Tour value: " << Results.back().Randinsert_val << endl;
                    // for (size_t i = 0; i < Grdy.Route.size(); i++)
                    // {
                    //     cout << Label[(G.target(Grdy.Route[i]))];
                    // }
                    // cout << " Tour Length: " << Grdy.Length << endl;
                    // cout << "Time: " << timer.realTime() << endl;
                }
#pragma endregion
            }
        }
#pragma region Branch and Cut
        int BnC = count(what_to_do.begin(), what_to_do.end(), "BnC");
        if (BnC)
        {
            {
                timer.restart();
                std::cout << "\nBranch and Cut and Price..." << "\n";
                BnCnP::Algorithm<BnCnP::Logging> ALG(G, Label, weight);
                Results.back().BnC_val = ALG.solve();
                Results.back().BnC_Time = timer.realTime();
                Results.back().BnC_solved = ALG.get_OPTsolved();
                std::cout << "Opt value: " << Results.back().BnC_val << endl;
                vector<ListDigraph::Arc> Tour = ALG.OPTroute();
                // for (size_t i = 0; i < Tour.size(); i++)
                // {
                //     cout << "->" << Label[(G.target(Tour[i]))];
                // }
                // cout << endl;
                // cout << "Time: " << timer.realTime() << endl;
            }
        }
#pragma endregion

#pragma region HeldKarp
        int HeldKarp = count(what_to_do.begin(), what_to_do.end(), "HK");
        if (HeldKarp)
        {
            {
                // ListDigraph::Node o;
                // for (ListDigraph::NodeIt v(G); v != INVALID; ++v)
                // {
                //     if (Label[v] == 0)
                //     {
                //         o = v;
                //         break;
                //     }
                // }

                // ListDigraph::Node osink = G.addNode();
                // Label[osink] = countNodes(G) - 1;
                // for (ListDigraph::InArcIt i(G, o); i != INVALID; ++i)
                // {
                //     ListDigraph::Node source = G.source(i);
                //     ListDigraph::Arc a = G.addArc(source, osink);
                //     weight[a] = weight[i];
                // }

                std::cout << "\nHeld-Karp SST..." << endl;
                // timer.restart();
                // Heldkarp::Heldkarp<Heldkarp::Silent, bound::SST> SST(G, Label, weight);
                // std::cout << "Opt value: " << SST.solve() << " Best route: ";
                // SST.printroute();
                // cout << "Time: " << timer.realTime() << endl;
            }
        }
#pragma endregion

        Results.back().finalize();
        // cout << "NN approximation percentage for: " << num << ": " << Results.back().NN_percentage << " Greedy percentage: " << Results.back().Grdy_percentage << " Greedy time: " << Results.back().Grdy_time << endl;
        return 0;
    }
}
#endif
