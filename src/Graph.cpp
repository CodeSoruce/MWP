#include "Graph.h"

#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include "KPlex_BB_matrix.h"
#include <iomanip>

using namespace std;

Graph::Graph(const char* _dir, const int _K, const int _Q) {
    dir = string(_dir);
    K = _K;
    q = _Q;
    mw = 0;
    n = m = 0;
    max_weight = 0;
    max_color = 0;
    pstart = nullptr;
    edges = nullptr;
    edgelist_pointer = nullptr;
    weight = nullptr;
    old_w = nullptr;
    color = nullptr;
    uw = nullptr;

    kplex.clear();
}

Graph::~Graph() {
    if (pstart != nullptr) {
        delete[] pstart;
        pstart = nullptr;
    }
    if (edges != nullptr) {
        delete[] edges;
        edges = nullptr;
    }
    if (edgelist_pointer != nullptr) {
        delete[] edgelist_pointer;
        edgelist_pointer = nullptr;
    }
    if (weight != nullptr) {
        delete[] weight;
        weight = nullptr;
    }
    if (color != nullptr) {
        delete[] color;
        color = nullptr;
    }
    if (uw != nullptr) {
        delete[] uw;
        uw = nullptr;
    }
    if (old_w != nullptr) {
        delete[] old_w;
        old_w = nullptr;
    }
}


void Graph::read_graph() {
    printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
    FILE* f = Utility::open_file((string("../dataset/") + dir + string("/b_degree.bin")).c_str(), "rb");

    ui tt;
    fread(&tt, sizeof(int), 1, f);
    if (tt != sizeof(int)) {
        printf("sizeof int is different: edge.bin(%d), machine(%d)\n", tt, (int)sizeof(int));
        return;
    }
    fread(&n, sizeof(int), 1, f);
    fread(&m, sizeof(int), 1, f);

    printf("\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());

    ui* degree = new ui[n];
    fread(degree, sizeof(int), n, f);

    fclose(f);

    f = Utility::open_file((string("../dataset/") + dir + string("/b_adj.bin")).c_str(), "rb");

    if (pstart == nullptr) pstart = new ept[n + 1];
    if (edges == nullptr) edges = new ui[m];

    pstart[0] = 0;
    for (ui i = 0; i < n; i++) {
        if (degree[i] > 0) {
            fread(edges + pstart[i], sizeof(int), degree[i], f);

            ui* buff = edges + pstart[i];
            sort(buff, buff + degree[i]);
            ui idx = 0;

            std::unordered_set<ui> processed_edges;
            for (ui j = 0; j < degree[i]; j++) {
                if (buff[j] >= n) printf("vertex id %u wrong\n", buff[j]);

                if (buff[j] == i) continue;
                if (j > 0 && buff[j] == buff[j - 1]) {
                    printf("Duplicate edge detected between vertices %u and %u\n", i, buff[j]);
                    continue;
                }
                if (processed_edges.count(buff[j])) {
                    printf("Repeated edge between vertices %u and %u\n", i, buff[j]);
                    continue;
                }

                buff[idx++] = buff[j];
                processed_edges.insert(buff[j]);
            }
            degree[i] = idx;
        }
        pstart[i + 1] = pstart[i] + degree[i];
    }
  
    

    fclose(f);

    delete[] degree;

    f = Utility::open_file((string("../dataset/") + dir + string("/b_weight.bin")).c_str(), "rb");
    weight = new ui[n];
    old_w = new ui[n];
    fread(weight, sizeof(ui), n, f);
    for (ui i = 0; i < n; i++) {
        old_w[i] = weight[i];
        if (mw < weight[i]) mw = weight[i];
    }
    fclose(f);
}


void Graph::exact(ui mode) {
    auto start = std::chrono::high_resolution_clock::now();
    assert(K > 0);
    if (K <= 1) {
        printf("\tFor k <= 1, please invoke clique computation algorithms\n");
        return;
    }
    kplex.clear();
    coloring();
    heuristic_kplex_max_degree(50);
    
    auto end_h = std::chrono::high_resolution_clock::now();
    auto elapsed_h = std::chrono::duration_cast<std::chrono::microseconds>(end_h - start);
    printf("Time: %lld (microseconds)\n", elapsed_h.count());
    if(mode == 1 || mode == 2){
        ui new_q = (max_weight + mw - 1) / mw;
        if (new_q > q) q = new_q;
    }
    ui* rid = new ui[n];
    ui* out_mapping = new ui[n];
    ui* peel_sequence = new ui[n];
    for (ui i = 0; i < n; i++) peel_sequence[i] = i;

    uw = new ui[n];
    ui old = max_weight;
    for (ui i = 0; i < n; i++) uw[i] = count_uw(i, pstart);
    ui* rid_old = new ui[n];
    d_uw_prun(n, m, peel_sequence, out_mapping, rid, rid_old, pstart, edges, true);
    auto end_p = std::chrono::high_resolution_clock::now();
    auto elapsed_p = std::chrono::duration_cast<std::chrono::microseconds>(end_p - end_h);
    printf("Prun Time: %lld (microseconds)\n", elapsed_p.count());

    vector<ui> ids;
    vector<ui> vw;
    vector<ui> vc;
    vector<pair<ui, ui> > vp;
    ui* peel_sequence_rid = new ui[n];
    for (ui i = 0; i < n; i++) peel_sequence_rid[peel_sequence[i]] = i;
    KPLEX_BB_MATRIX* kplex_solver_m = new KPLEX_BB_MATRIX();
    kplex_solver_m->allocateMemory(n);
    ui* pend = new ept[n + 1];
    reorganize_adjacency_lists(n, peel_sequence, rid, pstart, pend, edges);

    ui* degree = new ui[n];
    ui* vis = new ui[n];
    memset(vis, 0, sizeof(ui) * n);

    for (ui i = 0; i < n; i++) {
        ui u = peel_sequence[i];
        fflush(stdout);
        ui m_vw = 0, q_s = q;
        extract_subgraph_with_prune(u, peel_sequence_rid, degree, ids, rid, vp, vw, vc, vis, pstart, pend, edges, m_vw);
        if(mode == 1 || mode == 2){
            if(m_vw != mw){
                ui new_q = (max_weight + m_vw - 1) / m_vw;
                if (new_q > q_s) q_s = new_q;
            }
            if (ids.size() < q_s) continue;
            ui UB = count_w(ids, ids.size());
            if (UB <= max_weight) continue;
            ui UB1 = upperbound(u, ids, pstart);
            if (UB1 <= max_weight) continue;
            ui UB2 = findMaxWeightSum(ids, color, ids.size());
            if (UB2 <= max_weight) continue;
        }
        ui t_old = max_weight;
        kplex_solver_m->load_graph(ids.size(), vp, vw, vc);
        kplex_solver_m->mwp(K, kplex, true, max_weight, kplex.size(), q_s, mode);

        if (max_weight > t_old) {
            for (ui j = 0; j < kplex.size(); j++) kplex[j] = ids[kplex[j]];
            if(mode == 1 | mode == 2){
                ui new_q = (max_weight + mw - 1) / mw;
                if (new_q > q) q = new_q;
            }
        }
    }
    delete kplex_solver_m;
    auto end_s = std::chrono::high_resolution_clock::now();
    auto elapsed_s = std::chrono::duration_cast<std::chrono::microseconds>(end_s - end_p);
    printf("*** Search time: %lld (microseconds), Total Time: %lld (microseconds)\n", elapsed_s, std::chrono::duration_cast<std::chrono::microseconds>(end_s - start));
    if (max_weight > old) {
        output_one_kplex();
    }
    else {
        output_one_kplex_old(rid_old);
    }
    verify_kplex();
    if (max_weight > old) {
        for (ui i = 0; i < kplex.size(); i++) {
            assert(kplex[i] < n);
            kplex[i] = out_mapping[kplex[i]];
        }
    }
    cout << "===== k-MWP Result =====" << endl;
    cout << left << setw(10) << "Vertex" << setw(10) << "Weight" << endl;
    cout << "------------------------" << endl;

    int total_weight = 0;
    for(ui i = 0; i < kplex.size(); i++) {
        int w = old_w[kplex[i]];
        total_weight += w;
        cout << left << setw(10) << kplex[i] 
            << setw(10) << w << endl;
    }

    cout << "------------------------" << endl;
    cout << left << setw(10) << "k-MWP Size: " << setw(10) << kplex.size() << endl;
    cout << left << setw(10) << "Total Weight: " << setw(10) << total_weight << endl;
    cout << "========================" << endl;
    delete[] out_mapping;
    delete[] rid;
    delete[] rid_old;
    delete[] pend;
    delete[] peel_sequence;
    delete[] peel_sequence_rid;
    delete[] vis;
    delete[] degree;

}

ui Graph::count_w(vector<ui> set, ui num) {
    if (num == 0) return 0;
    ui w = weight[set[0]];
    for (ui i = 1; i < num; i++) w += weight[set[i]];
    return w;
}

ui Graph::count_uw(ui u, ept* pstart) {
    ui ub = weight[u];

    for (ui i = pstart[u]; i < pstart[u + 1]; i++) {
        ub += weight[edges[i]];
    }
    ub += (K - 1) * mw;
    return ub;
}

void Graph::heuristic_kplex_max_degree(ui processed_threshold) {
    assert(kplex.empty());
    ui* head = new ui[n];
    ui* next = new ui[n];
    ui* degree = new ui[n];
    ui* p = new ui[n];
    ui* vis = new ui[n];
    memset(vis, 0, sizeof(ui) * n);
    memset(p, 0, sizeof(ui) * n);

    int max_degree = 0;
    for (ui i = 0; i < n; i++) head[i] = n;
    for (ui i = 0; i < n; i++) {
        degree[i] = pstart[i + 1] - pstart[i];
        if (degree[i] > max_degree) max_degree = degree[i];
        next[i] = head[degree[i]];
        head[degree[i]] = i;
    }
    ui now_size = 0;
    int index = 0;
    for (ui processed_vertices = 0; max_degree >= q - K && processed_vertices < processed_threshold; processed_vertices++) {
        ui now_weight = 0;
        ui u = n;
        while (max_degree >= q - K && u == n) {
            for (ui v = head[max_degree]; v != n;) {
                ui tmp = next[v];
                if (degree[v] == max_degree) {
                    u = v;
                    head[max_degree] = tmp;
                    break;
                }
                if (degree[v] >= q - K) {
                    next[v] = head[degree[v]];
                    head[degree[v]] = v;
                }
                v = tmp;
            }
            if (u == n) {
                head[max_degree] = n;
                --max_degree;
            }
        }
        if (u == n) break;
        vis[u] = 1;
        vector<ui> vs, vv;
        int vv_size = 0;
        int vv_weight = 0;
        int vs_size = 0;
        for (ui i = pstart[u]; i < pstart[u + 1]; i++) if (!vis[edges[i]]) {
                --degree[edges[i]];
                vs.pb(edges[i]);
                vs_size++;
            }

        for (ui i = 0; i < vs_size; i++) {
            ui neighbor = vs[i];
            if (!vis[neighbor]) {
                vis[neighbor] = 2;
                vv.pb(neighbor);
                vv_size++;
                vv_weight += weight[neighbor];
            }
            for (ui j = pstart[neighbor]; j < pstart[neighbor + 1]; j++) {
                if (!vis[edges[j]]) {
                    vv.pb(edges[j]);
                    vis[edges[j]] = 2;
                    vv_size++;
                    vv_weight += weight[edges[j]];
                }
            }
        }
        for (ui i = 0; i < vv.size(); i++) {
            if (vis[vv[i]] == 2) vis[vv[i]] = 0;
        }
        vector<ui> vs_deg(vs_size);
        for (ui i = 0; i < vs.size(); i++) vis[vs[i]] = 2;
        for (ui i = 0; i < vs.size(); i++) {
            ui v = vs[i], d = 0;
            for (ui j = pstart[v]; j < pstart[v + 1]; j++) {
                if (vis[edges[j]] == 2) ++d;
            }
            vs_deg[i] = d;
        }
        for (ui i = 0; i < vs.size(); i++) vis[vs[i]] = 0;
        vector<ui> res; res.pb(u); now_weight += weight[u]; p[u]++;
        ui res_size = 1;
        while (vs_size > 0 && vv_weight + now_weight > max_weight && vv_size + res_size >= q) {
            ui idx = 0;
            for (ui i = 1; i < vs_size; i++) {
                if (vs_deg[i] > vs_deg[idx]) idx = i;
                else if (vs_deg[i] == vs_deg[idx] && weight[vs[i]] > weight[vs[idx]]) idx = i;
                else if (weight[vs[i]] == weight[vs[idx]] && degree[vs[i]] > degree[vs[idx]]) idx = i;
                else if (degree[vs[i]] == degree[vs[idx]] && p[vs[i]] < p[vs[idx]]) idx = i;
            }
            u = vs[idx];
            ui new_size = 0;
            for (ui i = pstart[u]; i < pstart[u + 1]; i++) if (!vis[edges[i]]) vis[edges[i]] = 2;
            for (ui i = 0; i < vs_size; i++) if (vis[vs[i]] == 2) {
                    if (i != new_size) swap(vs[new_size], vs[i]);
                    vs_deg[new_size] = vs_deg[i];
                    ++new_size;
                }
            for (ui i = pstart[u]; i < pstart[u + 1]; i++) if (vis[edges[i]] == 2) vis[edges[i]] = 0;

            res.pb(u); res_size++; now_weight += weight[u]; p[u]++;
            for (auto it = vv.begin(); it != vv.end(); ) {
                if (*it == u) {
                    *it = vv.back();
                    vv.pop_back();
                    vv_size--;
                    vv_weight -= weight[*it];
                    break;
                }
                ++it;
            }
            for (ui k = 0; k < new_size; k++) vis[vs[k]] = k + 2;
            for (ui j = new_size; j < vs_size; j++) {
                ui v = vs[j];
                for (ui k = pstart[v]; k < pstart[v + 1]; k++) {
                    if (vis[edges[k]] >= 2) --vs_deg[vis[edges[k]] - 2];
                }
            }
            for (ui k = 0; k < new_size; k++) vis[vs[k]] = 0;

            vs_size = new_size;
        }
        vector<ui> res_deg;
        res_deg.clear();
        for (ui j = 0; j < res_size; j++) vis[res[j]] = 2;
        for (auto it = vv.begin(); it != vv.end();) {
            ui v = *it, d = 0;
            for (ui k = pstart[v]; k < pstart[v + 1]; k++) {
                if (vis[edges[k]] == 2) ++d;
            }
            if (d <= res_size - K) {
                vv_weight -= weight[*it];
                *it = vv.back();  
                vv.pop_back();   
                vv_size--;
            }
            else {
                res_deg.pb(d);
                ++it;
            }
        }
        for (ui j = 0; j < res_size; j++) vis[res[j]] = 0;
        vis[res[0]] = 1;
        ui cut = K - 1;
        while (cut > 0 && vv_size + res_size >= q && now_weight + vv_weight > max_weight) {
            int idx = -1;
            for (ui j = 0; j < vv_size; j++) {
                if (res_deg[j] > res_size - K) {
                    idx = j; break;
                }
            }
            if (idx == -1) break;
            for (ui i = idx + 1; i < vv_size; i++) {
                if (res_deg[i] > res_deg[idx]) idx = i;
                else if (res_deg[i] == res_deg[idx] && weight[vv[i]] > weight[vv[idx]]) idx = i;
                else if (weight[vv[i]] == weight[vv[idx]] && degree[vv[i]] > degree[vv[idx]]) idx = i;
                else if (degree[vv[i]] == degree[vv[idx]] && p[vv[i]] < p[vv[idx]]) idx = i;
            }
            u = vv[idx];
            if (res_deg[idx] < res_size) cut--;
            res_size++; res.pb(u);  now_weight += weight[u]; p[u]++;
            if (K - (res_size - res_deg[idx]) < cut) cut = K - (res_size - res_deg[idx]);
            if (cut == 0) break;
            vv[idx] = vv.back();
            res_deg[idx] = res_deg[vv_size - 1];
            vv.pop_back();
            vv_size--;

            ui* rid_v = new ui[n];
            for (ui j = 0; j < vv_size; j++) {
                vis[vv[j]] = 2;
                rid_v[vv[j]] = j;
            }
            for (ui j = pstart[u]; j < pstart[u + 1]; j++) {
                if (vis[edges[j]] == 2) res_deg[rid_v[edges[j]]]++;
            }
            for (ui j = 0; j < vv_size; j++) vis[vv[j]] = 0;
            delete[]rid_v;
        }
        for (ui j = 0; j < res.size(); j++) vis[res[j]] = 0;
        vis[res[0]] = 1;

        if (now_weight > max_weight && res_size >= q) {
            kplex = res;
            now_size = res_size;
            max_weight = now_weight;
            index = processed_vertices;
        }
    }

    delete[] vis;
    delete[] head;
    delete[] next;
    delete[] degree;
    delete[] p;

    printf("*** Heuristic kplex size: %lu weight: %u, processed: %u, ", now_size, max_weight, index);
}

void Graph::d_uw_prun(ui& n, ept& m, ui* peel_sequence, ui* out_mapping, ui* rid, ui*& rid_old, ept*& pstart, ui*& edges, bool output) {
ui* degree = new ui[n];
    ui* vis = new ui[n];
    queue<int> Q;
    for (ui i = 0; i < n; i++) vis[i] = 0;
    for (ui i = 0; i < n; i++) {
        degree[i] = pstart[i + 1] - pstart[i];
        if (degree[i] < q - K || uw[i] <= max_weight) {
            vis[i] = 1;
            Q.push(i);
        }
    }
    while (!Q.empty()) {
        int cur = Q.front();
        Q.pop();
        for (int i = pstart[cur]; i < pstart[cur + 1]; i++) {
            int neighbor = edges[i];
            if (!vis[neighbor]) {
                uw[neighbor] -= weight[cur];
                if (--degree[neighbor] < q - K || uw[neighbor] <= max_weight) {
                    Q.push(neighbor);
                    vis[neighbor] = 1;
                }
            }
        }
    }
    ui cnt = 0;

    for (ui i = 0; i < n; i++) if (!vis[i]) {
        rid[i] = cnt;
        rid_old[i] = rid[i];
        out_mapping[cnt] = i;
        ++cnt;
    }

    if (cnt != n) {
        cnt = 0;
        ui pos = 0;
        for (ui i = 0; i < n; i++)
            if (!vis[i]) {
                ept t_start = pstart[i];
                pstart[cnt] = pos;
                for (ept j = t_start; j < pstart[i + 1]; j++)
                    if (!vis[edges[j]]) {
                        edges[pos++] = rid[edges[j]];
                    }
                ++cnt;
            }
        pstart[cnt] = pos;

        ui new_mw = 0;
        assert(core[peel_sequence[n - cnt - 1]] == 0 || core[peel_sequence[n - cnt - 1]] + K <= kplex.size());
        assert(cnt == 0 || core[peel_sequence[n - cnt]] + K > kplex.size());
        if (cnt > 0) {
            ui* weight_new = new ui[cnt + 1];
            ui* color_new = new ui[cnt + 1];
            for (ui i = 0; i < cnt; i++) {
                weight_new[i] = weight[out_mapping[i]];
                color_new[i] = color[out_mapping[i]];
                if (weight_new[i] > new_mw) new_mw = weight_new[i];
            }
            if (new_mw < mw) mw = new_mw;
            ui index = 0;
            for (ui i = 0; i < n; i++) if (!vis[i]) {
                peel_sequence[index++] = rid[i];
            }


            std::sort(peel_sequence, peel_sequence + cnt, [&](ui a, ui b) {
                if (degree[out_mapping[a]] != degree[out_mapping[b]])
                    return degree[out_mapping[a]] < degree[out_mapping[b]];
                return uw[out_mapping[a]] < uw[out_mapping[b]];
                });


            delete[] weight;
            weight = weight_new;
            delete[] color;
            color = color_new;


        }
        if (pos > 0 && pos < m / 2) {
            ept* pstart_new = new ept[cnt + 1];
            ui* edges_new = new ui[pos];
            memcpy(pstart_new, pstart, sizeof(ept) * (cnt + 1));
            memcpy(edges_new, edges, sizeof(ui) * pos);
            delete[] pstart;
            pstart = pstart_new;
            delete[] edges;
            edges = edges_new;

        }
        if (cnt == 0) {
            for (ui i = 0; i < n; i++) {
                rid[i] = i;
                rid_old[i] = i;
                out_mapping[i] = i;
            }
        }
        n = cnt;
        m = pos;
    }


    if (output) printf("*** After degree shrink: n = %s, m = %s (undirected) ", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m / 2).c_str());
    delete[] vis;
    delete[] degree;
}

void Graph::reorganize_adjacency_lists(ui n, ui* peel_sequence, ui* rid, ui* pstart, ui* pend, ui* edges) {
    for (ui i = 0; i < n; i++) rid[peel_sequence[i]] = i;
    for (ui i = 0; i < n; i++) {
        ui& end = pend[i] = pstart[i];
        for (ui j = pstart[i]; j < pstart[i + 1]; j++) if (rid[edges[j]] > rid[i]) edges[end++] = edges[j];
    }
    for (ui i = n; i > 0; i--) {
        ui u = peel_sequence[i - 1];
        for (ui j = pstart[u]; j < pend[u] && rid[edges[j]] > rid[u]; j++) {
            ui v = edges[j];
            edges[pend[v]++] = u;
            assert(pend[v] <= pstart[v + 1]);
        }
    }
    for (ui i = 0; i < n; i++) {
        ui& end = pend[i] = pstart[i];
        while (end < pstart[i + 1] && rid[edges[end]] > rid[i]) ++end;
    }
}


void Graph::extract_subgraph_with_prune(ui u, const ui* p_rid, ui* degree, vector<ui>& ids, ui* rid, vector<pair<ui, ui> >& vp, vector<ui>& vw, vector<ui>& vc, ui* exists, ept* pstart, ui* pend, ui* edges, ui &m_vw) {
    ids.clear(); vp.clear(); vw.clear(); vc.clear();
    int degree_threshold = q - K, triangle_threshold = q - 2 * K, cn_threshold = q + 2 - 2 * K;

    ids.push_back(u); exists[u] = 1;
    for (ept i = pstart[u]; i < pend[u]; i++) {
        assert(p_rid[edges[i]] > p_rid[u]);
        ids.push_back(edges[i]); exists[edges[i]] = 2;
    }
    assert(pend[u] >= pstart[u + 1] || p_rid[edges[pend[u]]] < p_rid[u]);

    ui* Q = rid;
    ui Q_n = 0;
    for (ui i = 1; i < ids.size(); i++) {
        ui v = ids[i];
        degree[v] = 0;
        for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++) {
            if (exists[edges[j]]) ++degree[v];
        }
        if (degree[v] - triangle_threshold < 0) Q[Q_n++] = v;
    }
    for (ui i = 0; i < Q_n; i++) {
        ui v = Q[i];
        exists[v] = 3;
        for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++) if (exists[edges[j]] == 2) {
                if (degree[edges[j]] - triangle_threshold == 0) Q[Q_n++] = edges[j];
                --degree[edges[j]];
            }
    }
    assert(Q_n < ids.size());
    if (ids.size() - Q_n - 1 - degree_threshold < 0) {
        for (ui i = 0; i < ids.size(); i++) exists[ids[i]] = 0;
        ids.clear();
        return;
    }
    ui old_size = ids.size();
    for (ui i = 1; i < old_size; i++) if (exists[ids[i]] == 2) {
            ui v = ids[i];
            for (ept j = pstart[v]; j < pstart[v + 1] && p_rid[edges[j]] > p_rid[u]; j++) {
                if (!exists[edges[j]]) {
                    ids.push_back(edges[j]);
                    exists[edges[j]] = 1;
                    degree[edges[j]] = 1;
                }
                else ++degree[edges[j]];
            }
        }

    ui new_size = 1;
    for (ui i = 1; i < old_size; i++) {
        if (exists[ids[i]] == 3) exists[ids[i]] = 0;
        else ids[new_size++] = ids[i];
    }
    assert(new_size + Q_n == old_size);
    for (ui i = old_size; i < ids.size(); i++) {
        if (degree[ids[i]] - cn_threshold < 0) exists[ids[i]] = 0;
        else ids[new_size++] = ids[i];
    }
    ids.resize(new_size);
    for (ui i = 0; i < ids.size(); i++) rid[ids[i]] = i;
    for (ui i = 0; i < ids.size(); i++) {
        ui v = ids[i];
        vw.push_back(weight[v]);
        vc.push_back(color[v]);
        if(weight[v] > m_vw) m_vw = weight[v];
        for (ept j = pstart[v]; j < pstart[v + 1]; j++) if (p_rid[edges[j]] > p_rid[v] && exists[edges[j]]) {
                vp.push_back(make_pair(rid[v], rid[edges[j]]));
            }
    }
    for (ui i = 0; i < ids.size(); i++) exists[ids[i]] = 0;
}

ui Graph::upperbound(ui u, vector<ui> ids, ept* pstart) {
    ui* vis_i = new ui[n];
    ui ub = weight[u];
    std::priority_queue<int, std::vector<int>, std::greater<int>> minHeap;
    ui size = ids.size();
    memset(vis_i, 0, sizeof(ui) * n);
    for (ui i = pstart[u]; i < pstart[u + 1]; i++) vis_i[edges[i]] = 1;
    for (ui i = 1; i < size; i++) {
        ui v = ids[i];
        if (vis_i[v]) {
            ub += weight[v];
            if (ub > max_weight) {
                delete[] vis_i;
                while (!minHeap.empty()) {
                    minHeap.pop();
                }
                return ub;
            }
        }
        else {
            if (minHeap.size() < K - 1) {
                minHeap.push(weight[v]);
            }
            else if (weight[v] > minHeap.top()) {
                if (weight[v] * (K - 1) + ub > max_weight) {
                    while (!minHeap.empty()) {
                        minHeap.pop();
                    }
                    delete[]vis_i;
                    ub += weight[v] * (K - 1);
                    return ub;
                }
                minHeap.pop();
                minHeap.push(weight[v]);
            }

        }
    }

    while (!minHeap.empty()) {
        ub += minHeap.top();
        minHeap.pop();
    }

    delete[]vis_i;

    return ub;
}

void Graph::coloring() {
    int* cvis = new int[n];
    color = new ui[n];
    int* head = new int[n];
    int* nxt = new int[n];
    int max_degree = 0;
    int* degree = new int[n];


    for (int i = 0; i < n; i++) {
        degree[i] = pstart[i + 1] - pstart[i];
        head[i] = n;
    }

    for (int i = 0; i < n; i++) {
        nxt[i] = head[degree[i]];
        head[degree[i]] = i;
        if (degree[i] > max_degree) max_degree = degree[i];
    }
    delete[] degree;

    for (int i = 0; i < n; i++) cvis[i] = 0;
    for (int i = 0; i < n; i++) color[i] = n;
    max_color = 0;
    for (int ii = max_degree; ii >= 1; ii--) {
        for (int jj = head[ii]; jj != n; jj = nxt[jj]) {
            int u = jj;
            for (int j = pstart[u]; j < pstart[u + 1]; j++) {
                int c = color[edges[j]];
                if (c != n) {
                    cvis[c] = 1;
                }
            }

            for (int j = 0;; j++) {
                if (!cvis[j]) {
                    color[u] = j;
                    if (j > max_color) max_color = j;
                    break;
                }
            }
            for (int j = pstart[u]; j < pstart[u + 1]; j++) {
                int c = color[edges[j]];
                if (c != n) cvis[c] = 0;
            }
        }
    }
    max_color++;
    delete[] cvis;
    delete[] head;
    delete[] nxt;
}

ui Graph::findMaxWeightSum(vector<ui> ids, ui* color, ui size) {
    unordered_map<ui, vector<ui>> colorTopWeights;

    for (int i = 0; i < size; ++i) {
        ui u = ids[i];
        ui colorId = color[u];
        ui nodeWeight = weight[u];

        auto& topWeights = colorTopWeights[colorId];
        if (topWeights.size() < K) {
            topWeights.push_back(nodeWeight);
            std::push_heap(topWeights.begin(), topWeights.end(), std::greater<ui>());
        }
        else if (nodeWeight > topWeights.front()) {
            std::pop_heap(topWeights.begin(), topWeights.end(), std::greater<ui>());
            topWeights.back() = nodeWeight;
            std::push_heap(topWeights.begin(), topWeights.end(), std::greater<ui>());
        }
    }

    ui totalWeightSum = 0;
    for (const auto& pair : colorTopWeights) {
        const auto& topWeights = pair.second;
        for (ui weight : topWeights) {
            totalWeightSum += weight;
            if (totalWeightSum > max_weight) {
                colorTopWeights.clear();
                return totalWeightSum;
            }
        }
    }
    return totalWeightSum;
}


void Graph::output_one_kplex() {
    FILE* fout = Utility::open_file("kplexes.txt", "w");
    fprintf(fout, "%lu\n", kplex.size());
    sort(kplex.begin(), kplex.end());
    for (ui i = 0; i < kplex.size(); i++) fprintf(fout, " %u", kplex[i]);
    fprintf(fout, "\n");
    fclose(fout);
}

void Graph::output_one_kplex_old(ui* rid_old) {
    FILE* fout = Utility::open_file("kplexes.txt", "w");
    fprintf(fout, "%lu\n", kplex.size());
    sort(kplex.begin(), kplex.end());
    for (ui i = 0; i < kplex.size(); i++) fprintf(fout, " %u", rid_old[kplex[i]]);
    fprintf(fout, "\n");
    fclose(fout);
}

void Graph::verify_kplex() {
    char* vis = new char[n];
    memset(vis, 0, sizeof(char) * n);

    FILE* fin = Utility::open_file("kplexes.txt", "r");

    ui kplex_size = n, kplex_n, idx = 0;
    char ok = 1;
    while (fscanf(fin, "%u", &kplex_n) == 1) {
        ++idx;
        if (kplex_size == n) {
            kplex_size = kplex_n;
        }
        if (kplex_n != kplex_size) printf("!!! WA k-plex size: %u!\n", kplex_n);
        vector<ui> kplex;
        for (ui i = 0; i < kplex_n; i++) {
            ui tmp;
            fscanf(fin, "%u", &tmp);
            kplex.pb(tmp);
        }

        for (ui i = 0; i < kplex.size(); i++) {
            if (vis[kplex[i]]) {
                //printf("WA k-plex! Duplicate vertex: %u\n", idx);
                ok = 0;
                break;
            }
            vis[kplex[i]] = 1;
        }
        for (ui i = 0; i < kplex.size(); i++) {
            ui d = 0;
            for (ui j = pstart[kplex[i]]; j < pstart[kplex[i] + 1]; j++) if (vis[edges[j]]) ++d;
            //for(ui i = pstart[1]; i < pstart[2]; i++) {cout<<edges[i]<<" ";}
            if (d + K < kplex.size()) {
                ok = 0;
                //printf("WA k-plex! Not enough neighbors!\n");
            }
        }
        for (ui i = 0; i < kplex.size(); i++) vis[kplex[i]] = 0;
    }
    if (ok) printf("Correct k-MWP!\n");
    if (!ok) printf("k-MWP not correctly read!\n");
    fclose(fin);

    delete[] vis;
}
