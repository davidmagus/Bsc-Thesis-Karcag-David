#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include "Heuristic.h"
#include "Heldkarp.h"
#include "TSP_generator.h"
#include "test_tools.h"
#include "BnC.h"
#include <string>
using namespace std;
using namespace lemon;
using namespace test_tools;

template <typename T>
double get_average(const std::vector<entry> &Results, T entry::*field)
{
    double sum = 0;
    int cnt = 0;

    for (const auto &a : Results)
    {
        // Elérjük az adott mezőt az adattag-mutató segítségével
        if (a.*field)
        {
            sum += a.*field;
            cnt++;
        }
    }

    return (cnt > 0) ? (sum / cnt) : 0.0;
}

int main(int argc, char *argv[])
{
    vector<string> what_to_do = {"BnC", "Heu", "HK"};
    std::vector<int> runs;

    for (int i = 1; i < argc; i++)
    {
        // argv[i]-ből stringet csinálunk, majd int-té alakítjuk
        runs.push_back(std::stoi(argv[i]));
    }

    for (size_t i = 0; i < runs.size(); ++i)
    {
        cout << endl
             << "Run number: " << i << ", number of nodes: " << runs[i] << "\n";
        one_run(runs[i], i, what_to_do);
    }

    double avarage_BnCruntime = get_average(Results, &entry::BnC_Time);
    double avarage_NN_percentage = get_average(Results, &entry::NN_percentage);
    double avarage_RNN_percentage = get_average(Results, &entry::RNN_percentage);
    double avarage_Greedy_percentage = get_average(Results, &entry::Grdy_percentage);
    double avarage_Maxins_percentage = get_average(Results, &entry::Maxinsert_percentage);
    double avarage_Minins_percentage = get_average(Results, &entry::Mininsert_percentage);
    double avarage_Randins_percentage = get_average(Results, &entry::Randinsert_percentage);

    cout << "Avarages: " << endl
         << " Average Branch and Cut time: " << avarage_BnCruntime << endl
         << " Average NN approx: " << avarage_NN_percentage << endl
         << " Average RNN approx: " << avarage_RNN_percentage << endl
         << " Average Greedy approx: " << avarage_Greedy_percentage << endl
         << " Average Max Insert approx: " << avarage_Maxins_percentage << endl
         << " Average Min Insert approx: " << avarage_Minins_percentage << endl
         << " Average Random Insert approx: " << avarage_Randins_percentage << endl;
}