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
#include <fstream>  // Fájlba íráshoz (ofstream)
#include <iomanip>  // Formázáshoz (std::fixed, std::setprecision)
#include <ctime>    // Rendszeridő lekéréséhez (std::time, std::strftime)
using namespace std;
using namespace lemon;

namespace test_tools
{

    struct entry
    {
        int n;
        int seed;

        double NN_val, NN_time, NN_percentage;
        double RNN_val, RNN_time, RNN_percentage;
        double Grdy_val, Grdy_time, Grdy_percentage;
        double Maxinsert_val, Maxinsert_time, Maxinsert_percentage;
        double Mininsert_val, Mininsert_time, Mininsert_percentage;
        double Randinsert_val, Randinsert_time, Randinsert_percentage;

        bool BnC_solved;
        double BnC_val;
        double BnC_Time;

        bool HK_solved;
        double HK_val;
        double HK_Time;

        entry(
            const int n,
            const int seed = 0,

            const double NN_val = 0,
            double NN_time = 0,

            double RNN_val = 0,
            double RNN_time = 0,

            double Grdy_val = 0,
            double Grdy_time = 0,

            double Maxinsert_val = 0,
            double Maxinsert_time = 0,

            double Mininsert_val = 0,
            double Mininsert_time = 0,

            double Randinsert_val = 0,
            double Randinsert_time = 0,

            bool BnC_solved = false,
            double BnC_val = 0,
            double BnC_Time = 0,

            bool HK_solved = false,
            double HK_val = 0,
            double HK_Time = 0)
            : n(n), seed(seed),
              NN_val(NN_val), NN_time(NN_time),
              RNN_val(RNN_val), RNN_time(RNN_time),
              Grdy_val(Grdy_val), Grdy_time(Grdy_time),
              Maxinsert_val(Maxinsert_val), Maxinsert_time(Maxinsert_time),
              Mininsert_val(Mininsert_val), Mininsert_time(Mininsert_time),
              Randinsert_val(Randinsert_val), Randinsert_time(Randinsert_time),
              BnC_solved(BnC_solved), BnC_val(BnC_val), BnC_Time(BnC_Time),
              HK_solved(HK_solved), HK_val(HK_val), HK_Time(HK_Time)
        {
            finalize();
        }

        void finalize()
        {
            if (BnC_val > 1e-9 && BnC_solved) // Osztás 0-val elleni védelem
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
                NN_percentage = RNN_percentage = Grdy_percentage =
                    Maxinsert_percentage = Mininsert_percentage = Randinsert_percentage = 0;
            }
        }

        void save() const
        {
            std::string filename = "results.txt";

            std::ifstream testFile(filename);
            bool isEmpty = !testFile.is_open() || testFile.peek() == std::ifstream::traits_type::eof();
            testFile.close();

            std::ofstream file(filename, std::ios::app);
            if (!file.is_open())
            {
                std::cerr << "Hiba: A results.txt nem nyithato meg!" << std::endl;
                return;
            }

            if (isEmpty)
            {
                // A fejlécbe is bekerült a Seed
                file << "Date\tn\tSeed\tNN_val\tNN_time\tNN_%\tRNN_val\tRNN_time\tRNN_%\t"
                     << "Grdy_val\tGrdy_time\tGrdy_%\tMinI_val\tMinI_time\tMinI_%\t"
                     << "MaxI_val\tMaxI_time\tMaxI_%\tRandI_val\tRandI_time\tRandI_%\t"
                     << "BnC_S\tBnC_val\tBnC_Time\tHK_S\tHK_val\tHK_Time" << std::endl;
            }

            std::time_t now = std::time(nullptr);
            char date_buf[20];
            std::strftime(date_buf, sizeof(date_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

            file << std::fixed << std::setprecision(4);
            file << date_buf << "\t"
                 << n << "\t"
                 << seed << "\t" // Seed mentése
                 << NN_val << "\t" << NN_time << "\t" << NN_percentage << "\t"
                 << RNN_val << "\t" << RNN_time << "\t" << RNN_percentage << "\t"
                 << Grdy_val << "\t" << Grdy_time << "\t" << Grdy_percentage << "\t"
                 << Mininsert_val << "\t" << Mininsert_time << "\t" << Mininsert_percentage << "\t"
                 << Maxinsert_val << "\t" << Maxinsert_time << "\t" << Maxinsert_percentage << "\t"
                 << Randinsert_val << "\t" << Randinsert_time << "\t" << Randinsert_percentage << "\t"
                 << (BnC_solved ? "YES" : "NO") << "\t" << BnC_val << "\t" << BnC_Time << "\t"
                 << (HK_solved ? "YES" : "NO") << "\t" << HK_val << "\t" << HK_Time
                 << std::endl;

            file.close();
        }
    };

    inline std::vector<entry> Results;

    inline int one_run(int num, int seed = 0, vector<string> what_to_do = {"BnC"},int timelimit = 30, bool noisy = false)
    {
        Timer timer;
#pragma region Generating Examples
        RandTSP::Make_completeTSP(num, seed);
        // cout << "Graph Made in " << timer.realTime() << endl;
        timer.reset();

#pragma endregion
        std::cout << "Timelimit: " << timelimit << " seconds " << " Running: ";
        for (string x : what_to_do)
        {
            std::cout << x << ", ";
        }
        std::cout << std::endl;
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
        Results.emplace_back(num, seed);
        {
            if (count(what_to_do.begin(), what_to_do.end(), "Heu"))
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
        if (count(what_to_do.begin(), what_to_do.end(), "BnC"))
        {
            {
                timer.restart();
                std::cout << "\nBranch and Cut and Price 1..." << "\n";
                BnCnP::Algorithm<BnCnP::Silent, 1> ALG(G, Label, weight);
                Results.back().BnC_val = ALG.solve();
                Results.back().BnC_Time = timer.realTime();
                Results.back().BnC_solved = ALG.get_OPTsolved();
                std::cout << "Opt value: " << Results.back().BnC_val << " Opt: " << ALG.get_OPTsolved() << endl;
                if (noisy)
                {
                    ALG.printroute();
                }
                // cout << endl;
                // cout << "Time: " << timer.realTime() << endl;
            }
            {
                timer.restart();
                std::cout << "\nBranch and Cut and Price 2..." << "\n";
                BnCnP::Algorithm<BnCnP::Silent, 2> ALG(G, Label, weight, timelimit);
                Results.back().BnC_val = ALG.solve();
                Results.back().BnC_Time = timer.realTime();
                Results.back().BnC_solved = ALG.get_OPTsolved();
                std::cout << "Opt value: " << Results.back().BnC_val << " Opt: " << ALG.get_OPTsolved() << endl;
                if (noisy)
                {
                    ALG.printroute();
                }
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
                // // timer.restart();
                Heldkarp::Heldkarp<Heldkarp::Silent> HK(G, Label, weight);
                std::cout << "Opt value: " << HK.solve() << " Best route: ";
                if (noisy)
                {
                    HK.printroute();
                }
                cout << endl
                     << endl;
                // // HK.printroute();
                // cout << "Time: " << timer.realTime() << endl;
            }
        }
#pragma endregion

        Results.back().finalize();
        Results.back().save();
        // cout << "NN approximation percentage for: " << num << ": " << Results.back().NN_percentage << " Greedy percentage: " << Results.back().Grdy_percentage << " Greedy time: " << Results.back().Grdy_time << endl;
        return 0;
    }
}
#endif