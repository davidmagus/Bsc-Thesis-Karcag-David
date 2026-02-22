#include <lemon/list_graph.h>
using namespace lemon;

#ifndef subTask_bounds_H
#define subTask_bounds_H

typedef double (*BoundFunc)(const lemon::ListDigraph&, const ListDigraph::ArcMap<double>&);

double SST(const ListDigraph& G, const ListDigraph::ArcMap<double>& weight);

#endif