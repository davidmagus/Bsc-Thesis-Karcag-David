#ifndef Heuristic_H
#define Heuristic_H
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
using namespace std;
using namespace lemon;

namespace Heuristic
{
    class Nearest_Neighbour
    {
    public:
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        const ListDigraph::ArcMap<double> &weight;
        vector<ListDigraph::Arc> Route;
        double Length;
        Nearest_Neighbour(
            const ListDigraph &_G,
            const ListDigraph::NodeMap<int> &_Label,
            const ListDigraph::ArcMap<double> &_weight,
            const int root = 0);
    };

    class Repetitive_Nearest_Neighbour
    {
    public:
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        const ListDigraph::ArcMap<double> &weight;
        vector<ListDigraph::Arc> Route;
        double Length;
        Repetitive_Nearest_Neighbour(
            const ListDigraph &_G,
            const ListDigraph::NodeMap<int> &_Label,
            const ListDigraph::ArcMap<double> &_weight);
    };
    class Greedy
    {
    public:
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        const ListDigraph::ArcMap<double> &weight;
        vector<ListDigraph::Arc> Orderedarcs;
        vector<ListDigraph::Arc> Route;
        double Length;
        Greedy(
            const ListDigraph &_G,
            const ListDigraph::NodeMap<int> &_Label,
            const ListDigraph::ArcMap<double> &_weight);
        ~Greedy() = default;
    };

    class Insertion_Methods
    {
    public:
        int n;
        const ListDigraph &G;
        const ListDigraph::NodeMap<int> &Label;
        const ListDigraph::ArcMap<double> &weight;
        list<ListDigraph::Node> Unseen;
        vector<ListDigraph::Arc> Route;
        double Length;

        Insertion_Methods(
            const ListDigraph &_G,
            const ListDigraph::NodeMap<int> &_Label,
            const ListDigraph::ArcMap<double> &_weight);

        virtual ~Insertion_Methods() = default;

        
        virtual void Insert(lemon::ArcLookUp<ListDigraph>& lookUp) = 0;
        virtual double Run() = 0;
    };

    class Max_insert : public Insertion_Methods
    {
    public:
        using Insertion_Methods::Insertion_Methods;
        void Insert(lemon::ArcLookUp<ListDigraph>& lookUp) override;
        double Run() override;
    };

    class Min_insert : public Insertion_Methods
    {
    public:
        using Insertion_Methods::Insertion_Methods;
        void Insert(lemon::ArcLookUp<ListDigraph>& lookUp) override;
        double Run() override;
    };

    class Rand_insert : public Insertion_Methods
    {
    public:
        using Insertion_Methods::Insertion_Methods;
        void Insert(lemon::ArcLookUp<ListDigraph>& lookUp) override;
        double Run() override;
    };
}
#endif