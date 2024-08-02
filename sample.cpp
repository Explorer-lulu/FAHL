//#include "labeling.hpp"
#include <head.h>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <random>
#include <chrono>
#include <omp.h>


using namespace FAHLquery;

class Labeling {
    std::vector< std::vector< std::vector<Vertex> > > label_v;     
    std::vector< std::vector< std::vector<Distance> > > label_d;   
    Vertex n;

public:
    Labeling(size_t n = 0) :
        label_v(n, std::vector< std::vector<Vertex> >(2)),
        label_d(n, std::vector< std::vector<Distance> >(2)),
        n(n) {}

    Distance shortest_spatial_distance_query(Vertex u, Vertex v, bool f = true) {
        Distance r = infty;
        for (size_t i=0, j=0; i < label_v[u][f].size() && j < label_v[v][!f].size();) {
            if (label_v[u][f][i] == label_v[v][!f][j]) {
                assert(label_d[u][f][i] < infty - label_d[v][!f][j]);
                r = std::min(r, label_d[u][f][i++] + label_d[v][!f][j++]);
            } else if (label_v[u][f][i] < label_v[v][!f][j]) ++i;
            else ++j;
        }
        return r;
    }
};

void generateRandomQueries(std::vector<std::pair<Vertex, Vertex>>& queries, int numberOfRandomQueries, int numberOfVertices) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distr(0, numberOfVertices-1); // define the range
	queries.reserve(numberOfRandomQueries);

	for (int i(0); i < numberOfRandomQueries; ++i) queries.push_back(std::make_pair(Vertex(distr(gen)), Vertex(distr(gen))));
}

int main(int argc, char *argv[]) {
	char DesFile="./roadsample/";
    char orderfile=DesFile+"flow_order";
	char queryfile=DesFile+"query";
	char *graph_file = NULL;
	char *label_file = NULL;
	int numberOfRandomQueries = 0;
	int num_threads = omp_get_max_threads();
	int argi;
	for (argi = 1; argi < argc; ++argi) {
		if (argv[argi][0] == '-') {
			if (!strcmp("--", argv[argi])) { ++argi; break; }
			else if (!strcmp("-l", argv[argi])) { if (++argi >= argc) usage(argv); label_file = argv[argi]; }
			else if (!strcmp("-n", argv[argi])) { if (++argi >= argc) usage(argv); numberOfRandomQueries = std::stoi(argv[argi]); }
			else usage(argv);
		} else if (graph_file == NULL) graph_file = argv[argi];
		else break;
	}
	if (argi != argc || !graph_file || !label_file) usage(argv);
    
	graph_file = DesFile + "graph_file";

	FRNGraph g;
	if (!g.read(graph_file)) {
		std::cerr << "read graph from file false" << graph_file << std::endl;
		std::exit(1);
	}
	std::cout << "graph has " << g.get_n() << " vertices and " << g.get_m() << " arcs" << std::endl;

	Labeling labels(g.get_n());
	std::vector<Vertex> order;

	if (!labels.read(label_file, g.get_n())) {
		std::cerr << "read labels from file false" << label_file << std::endl;
		std::exit(1);
	}
    
    tb1=std::chrono::high_resolution_clock::now();
	g.FAHLconOrderMT(orderfile);
	tb2=std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(tb2-tb1);
	runT= time_span.count();
	cout<<"FAHL build at "<<runT<<endl;
	g.CorCheckFAHL();

	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	std::vector<std::pair<Vertex, Vertex>> queries;
	generateRandomQueries(queries, numberOfRandomQueries, g.get_n());
	long double totalTime(0);
	for (std::pair<Vertex, Vertex> vpair : queries) {
		auto tq1 = high_resolution_clock::now();
		g.QueryFAHL(vpair.first, vpair.second);
		auto tq2 = high_resolution_clock::now();
		duration<double, std::milli> ms_double = tq2 - tq1;
		totalTime += ms_double.count();
	}
	std::cout  << numberOfRandomQueries <<  totalTime <<  (double) (totalTime / numberOfRandomQueries);
	g.FAHLindexsize();
	g.EffiCheckFAHL(queryfile);

}