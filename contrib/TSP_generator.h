#ifndef TSP_generator
#define TSP_generator

namespace RandTSP
{
    void Make_non_completeTSP(const int n, const double mratio, const int seed = 0);
    void Make_completeTSP(const int n, const int seed = 0);
}
#endif