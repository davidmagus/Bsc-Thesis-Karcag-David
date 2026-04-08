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

int main(int argc, char *argv[])
{
    std::vector<string> modes;
    int num = stoi(argv[1]);
    int seed;
    if (argc < 3)
    {
        seed = 0;
    }
    else
    {
        seed = stoi(argv[2]);
    }
    if (argc < 4)
    {
        one_run(num, seed);
    }
    else
    {

        for (int i = 3; i < argc; i++)
        {
            modes.push_back(argv[i]);
        }

        one_run(num, seed, modes);
    }
}