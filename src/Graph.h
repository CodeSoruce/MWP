#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Utility.h"
#include "LinearHeap.h"
using namespace std;

class Graph {
private:
    std::string dir; 
    ui q;
    ui n; 
    ept m; 
    ui K; 
    ui max_weight;
    ui max_color;
    ui mw;

    ept *pstart; 
    ui *edges; 
    ui *edgelist_pointer;

    ui *weight;
    ui *old_w;
    ui *color;
    ui *uw;

    std::vector<ui> kplex;

public:
    Graph(const char *_dir, const int _K, const int _Q); ;
    ~Graph() ;

    void read_graph() ;

    void exact(ui mode) ;
    ui count_w(std::vector<ui> set, ui num);

private:
    void reorganize_adjacency_lists(ui n, ui *peel_sequence, ui *rid, ui *pstart, ui *pend, ui *edges) ;

    void heuristic_kplex_max_degree(ui processed_threshold) ;
    void d_uw_prun(ui &n, ept &m, ui *peel_sequence, ui *out_mapping, ui *rid, ui *&rid_old, ept *&pstart,  ui *&edges, bool output) ;

    ui upperbound(ui u, std::vector<ui> ids, ept *pstart);

    void coloring();

    ui findMaxWeightSum(vector<ui> ids, ui* color, ui size);

    void output_one_kplex();

    void output_one_kplex_old(ui *rid_old);

    void verify_kplex();


    ui count_uw(ui u, ept *pstart);
    void extract_subgraph_with_prune(ui u, const ui *p_rid, ui *degree, std::vector<ui> &ids, ui *rid, std::vector<std::pair<ui, ui>> &vp, std::vector<ui> &
    vw, vector<ui> &vc, ui *exists, ept *pstart, ui *pend,ui *edges, ui& m_vw) ;
};
#endif
