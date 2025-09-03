#ifndef _KPLEX_BB_MATRIX_
#define _KPLEX_BB_MATRIX_

#include "Utility.h"
#include "LinearHeap.h"

using namespace std;

class KPLEX_BB_MATRIX {
private:
    ui n;

    char* matrix;
    long long matrix_size;

    ui* degree;
    ui* degree_in_S;
    ui* weight_now;
    ui* color_now;
    ui* color_cnt;
    ui mw;
    ui mc;

    ui K;
    ui q;
    ui* best_solution;
    ui best_solution_size;
    ui best_solution_weight;

    ui* neighbors;
    ui* nonneighbors;

    ui* SR; 
    ui* SR_rid; 
    std::queue<ui> Qv;
    ui* level_id;

    std::vector<std::pair<ui, ui> > vp;
    std::vector<ui> non_adj;

public:
    KPLEX_BB_MATRIX() {
        n = 0;
        matrix = nullptr;
        matrix_size = 0;

        degree = degree_in_S = nullptr;

        best_solution = nullptr;
        q = K = best_solution_weight = mc = mw = 0;

        neighbors = nonneighbors = nullptr;
        weight_now = nullptr;
        color_now = nullptr;
        color_cnt = nullptr;

        SR = SR_rid = nullptr;
        level_id = nullptr;
    }

    ~KPLEX_BB_MATRIX() {
        if (matrix != NULL) {
            delete[] matrix;
            matrix = NULL;
        }
        if (degree != NULL) {
            delete[] degree;
            degree = NULL;
        }
        if (degree_in_S != NULL) {
            delete[] degree_in_S;
            degree_in_S = NULL;
        }
        if (weight_now != NULL) {
            delete[] weight_now;
            weight_now = NULL;
        }
        if (color_now != NULL) {
            delete[] color_now;
            color_now = NULL;
        }
        if (color_cnt != NULL) {
            delete[] color_cnt;
            color_cnt = NULL;
        }
        if (best_solution != NULL) {
            delete[] best_solution;
            best_solution = NULL;
        }
        if (SR != NULL) {
            delete[] SR;
            SR = NULL;
        }
        if (SR_rid != NULL) {
            delete[] SR_rid;
            SR_rid = NULL;
        }
        if (neighbors != NULL) {
            delete[] neighbors;
            neighbors = NULL;
        }
        if (nonneighbors != NULL) {
            delete[] nonneighbors;
            nonneighbors = NULL;
        }
        if (level_id != NULL) {
            delete[] level_id;
            level_id = NULL;
        }
    }

    void allocateMemory(ui n) {
        if (n <= 0) return;

        matrix_size = 1;
        mc = mw = 0;
        matrix = new char[matrix_size];

        degree = new ui[n];
        degree_in_S = new ui[n];
        best_solution = new ui[n];
        SR = new ui[n];
        SR_rid = new ui[n];
        neighbors = new ui[n];
        nonneighbors = new ui[n];
        level_id = new ui[n];
        weight_now = new ui[n];
        color_now = new ui[n];
        color_cnt = new ui[n];
    }


    void load_graph(ui _n, const std::vector<std::pair<ui,ui> > &vp, std::vector<ui> &vw, std::vector<ui> &vc) {
        n = _n;
        if(((long long)n)*n > matrix_size) {
            do {
                matrix_size *= 2;
            } while(((long long)n)*n > matrix_size);
            delete[] matrix; matrix = new char[matrix_size];
        }

        memset(matrix, 0, sizeof(char)*((long long)n)*n);
        for(ui i = 0; i < n; i++) degree[i] = 0;
        for(ui i = 0;i < vp.size();i ++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
            ui a = vp[i].first, b = vp[i].second;
            degree[a] ++;
            degree[b] ++;
            long long ab = (long long)a * n + b;
            long long ba = (long long)b * n + a;
            if(matrix[ab]) cout<<a<<" "<<b<<endl;
            if(matrix[ab]) printf("Duplicate edge in KPLEX_BB_matrix.load_graph()\n");
            matrix[ab] = matrix[ba] = 1;
        }

        for(ui i = 0;i < n;i ++){
            weight_now[i] = vw[i];
            color_now[i] = vc[i];
            if(vw[i] > mw) mw = vw[i];
            if(vc[i] > mc) mc = vc[i];
        }
        mc++;
        color_cnt = new int[n];
        for (ui i = 0; i < n; i++) {
            color_cnt[i] = 0;
        }
        ui * flag_c = new int[mc];
        for(ui i = 0;i < n;i ++) {
            for (ui j = 0;j < mc;j ++) {flag_c[j] = 0;}
            long long un = (long long) i*n;
            char *t_matrix = matrix + un;
            for(ui j = 0;j < n;j ++) {
                if(t_matrix[j] ) {
                    if (flag_c[color_now[j]]==0) {
                        flag_c[color_now[j]] = weight_now[j];
                        color_cnt[i]+=weight_now[j];
                    } else {
                        if (weight_now[j] > flag_c[color_now[j]]) {
                            color_cnt[i] -= flag_c[color_now[j]];
                            flag_c[color_now[j]] = weight_now[j];
                            color_cnt[i]+=weight_now[j];
                        }
                    }
                }
            }
            for (ui j = 0;j < mc;j ++) {flag_c[j] = 0;}
        }
        delete[] flag_c;
    }

    void mwp(ui K_, std::vector<ui> &now_best, bool must_include_0, ui &old, ui old_size, ui size_limit, ui mode) {
        K = K_;
        q = size_limit;
        ui R_end = 0, R_w = 0, S_w = 0, S_end = 0;
        best_solution_size = old_size;
        best_solution_weight = old;
        if(mode == 0 || mode == 1) initialization_base(R_end, R_w, must_include_0);
        else initialization(R_end, R_w, must_include_0);
        if(R_end){
            if(mode == 0 || mode == 1) BB_search_base(S_end, R_end, S_w, R_w, 1, must_include_0, mode);
            else BB_search(S_end, R_end, S_w, R_w, 1, must_include_0);
        }
        if(best_solution_weight > old) {
            now_best.clear();now_best.shrink_to_fit();
            for(int i = 0;i < best_solution_size;i ++) now_best.push_back(best_solution[i]);
            old = best_solution_weight;
        }
    }

    ui count_ui(ui* set, ui begin, ui end) {
        if(end - begin == 0) return 0;
        ui w = 0;
        for(ui i = begin; i < end; i ++) w += weight_now[set[i]];
        return w;
    }

private:

    void initialization_base(ui &R_end, ui &R_w, bool must_include_0) {
        R_end = 0;
        for(ui i = 0;i < n;i ++) SR_rid[i] = n;
        for(ui i = 0;i < n;i ++) {
            SR[R_end] = i; SR_rid[i] = R_end;
            ++ R_end;
            R_w +=weight_now[i];
        }

        if(must_include_0&&SR_rid[0] == n) {
            R_end = 0;
            R_w = 0;
            return ;
        }

        for(ui i = 0;i < R_end;i ++) {
            ui u = SR[i];
            degree[u] = degree_in_S[u] = 0;
            long long un = (long long) u*n;
            char *t_matrix = matrix + un;
            for(ui j = 0;j < R_end;j ++) if(t_matrix[SR[j]]) ++ degree[u];
        }

        memset(level_id, 0, sizeof(ui)*n);
        for(ui i = 0;i < R_end;i ++) level_id[SR[i]] = n;
        if(!Qv.empty()) printf("!!! Something wrong. Qv must be empty in initialization\n");
        if(!remove_vertices_and_edges_with_prune(0, R_end, R_w, 0)) {
            R_end = 0;
            R_w = 0;
        }

    }

    void initialization(ui &R_end, ui &R_w, bool must_include_0) {
        ui *peel_sequence = neighbors;
        ui *core = nonneighbors;
        ui *vis = SR;
        memset(vis, 0, sizeof(ui)*n);
        ui max_core = 0, idx = n;
        for(ui i = 0;i < n;i ++) {
            ui u, min_degree = n, min_weight = mw + 1, min_cw = mc * mw + 1;
            for (ui j = 0;j < n;j ++) if(!vis[j]&&degree[j] < min_degree) min_degree = degree[j];
            for (ui j = 0;j < n;j ++) if(!vis[j]&&degree[j] == min_degree && weight_now[j] < min_weight) min_weight = weight_now[j];
            for (ui j = 0; j < n; j++) if (!vis[j] && degree[j] == min_degree && weight_now[j] == min_weight && color_cnt[j] < min_cw) {
                    min_cw = color_cnt[j];
                    u = j;
                }
            if(min_degree > max_core) max_core = min_degree;
            core[u] = max_core;
            peel_sequence[i] = u;
            vis[u] = 1;

            if(idx == n&&min_degree + K >= n - i) idx = i;
            long long un = (long long) u*n;
            char *t_matrix = matrix + un;
            for(ui j = 0;j < n;j ++) if(!vis[j]&&t_matrix[j]) -- degree[j];
        }
        ui now = count_ui(peel_sequence, idx, n);
        if(now > best_solution_weight && n - idx >= q) {
            delete []best_solution;
            best_solution = new ui[n-idx];
            best_solution_weight = now;
            best_solution_size = n - idx;
            for(ui i = idx;i < n;i ++) {
                best_solution[i-idx] = peel_sequence[i];
            }
            printf("*** CMS-Degeneracy: %u with size %u\n", best_solution_weight, best_solution_size);
            unsigned int new_q = (best_solution_weight + mw - 1) / mw;
            if (new_q > q) q = new_q;
        }

        R_end = 0;
        for(ui i = 0;i < n;i ++) SR_rid[i] = n;
        for(ui i = 0;i < n;i ++) if(core[i] >= q - K) {
                SR[R_end] = i; SR_rid[i] = R_end;
                ++ R_end;
                R_w +=weight_now[i];
            }

        if(must_include_0&&SR_rid[0] == n) {
            R_end = 0;
            R_w = 0;
            return ;
        }

        for(ui i = 0;i < R_end;i ++) {
            ui u = SR[i];
            degree[u] = degree_in_S[u] = 0;
            long long un = (long long) u*n;
            char *t_matrix = matrix + un;
            for(ui j = 0;j < R_end;j ++) if(t_matrix[SR[j]]) ++ degree[u];
        }

        memset(level_id, 0, sizeof(ui)*n);
        for(ui i = 0;i < R_end;i ++) level_id[SR[i]] = n;
        if(!Qv.empty()) printf("!!! Something wrong. Qv must be empty in initialization\n");
        if(!remove_vertices_and_edges_with_prune(0, R_end, R_w, 0)) {
            R_end = 0;
            R_w = 0;
        }

    }

    void store_solution_base(ui size, ui weight, ui mode) {
        delete []best_solution;
        best_solution = new ui[size];
        best_solution_size = size;
        best_solution_weight = weight;
        for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
        if(mode == 1){
            unsigned int new_q = (best_solution_weight + mw - 1) / mw;
            if (new_q > q) q = new_q;
        }
    }

    void store_solution(ui size, ui weight) {
        delete []best_solution;
        best_solution = new ui[size];
        best_solution_size = size;
        best_solution_weight = weight;
        for(ui i = 0;i < best_solution_size;i ++) best_solution[i] = SR[i];
        unsigned int new_q = (best_solution_weight + mw - 1) / mw;
        if (new_q > q) q = new_q;

    }

    void BB_search_base(ui S_end, ui R_end, ui S_w, ui R_w, ui level, ui choose_zero, ui mode) {
        assert(Qv.empty());
        ui old_S_end = S_end, old_R_end = R_end;

        if(choose_zero&&SR_rid[0] < R_end&&!move_u_to_S_with_prune_base(0, S_end, R_end, S_w, R_w, level)) {
            restore_SR_and_edges(S_end, R_end, R_w, S_w, old_S_end, old_R_end, level);
            return ;
        }
        if(R_end < q) {return;}
        if(mode == 1){
            ui m_vw = 0;
            for(ui i = S_end;i < R_end;i ++) if(weight_now[SR[i]] > m_vw) m_vw = weight_now[SR[i]];
            if(m_vw * (R_end - S_end) + S_w <= best_solution_weight) {return;}
        }
        else if(R_end < q) return;
        if(S_w > best_solution_weight && S_end >= 2 * K - 1) store_solution_base(S_end, S_w, mode);
        if(R_w > best_solution_weight&&is_kplex(R_end)) store_solution_base(R_end, R_w, mode);
        if(R_w <= best_solution_weight) return ;

        ui u = SR[S_end];
        ui  t_old_R_end = R_end;

        assert(Qv.empty());
        bool f = move_u_to_S_with_prune_base(u, S_end, R_end, S_w, R_w, level);
        if(f){
            BB_search_base(S_end, R_end, S_w, R_w, level+1, false, mode);
        }
        restore_SR_and_edges(S_end, R_end, R_w, S_w, S_end, t_old_R_end, level);
        bool succeed = remove_u_from_S_with_prune(S_end, R_end, R_w, S_w, level);
        if(succeed) succeed = remove_vertices_and_edges_with_prune(S_end, R_end, R_w, level);
        if(succeed) {
            BB_search_base(S_end, R_end, S_w, R_w, level+1, false, mode);
        }
        restore_SR_and_edges(S_end, R_end, R_w, S_w, old_S_end, old_R_end, level);

    }

    bool move_u_to_S_with_prune_base(ui u, ui &S_end, ui &R_end, ui &S_w, ui &R_w, ui level) {
        assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
        assert(degree_in_S[u] + K > S_end);
        if(SR_rid[u] != S_end) swap_pos(S_end, SR_rid[u]);
        ++ S_end;
        S_w +=weight_now[u];
        ui neighbors_n = 0, nonneighbors_n = 0;
        get_neighbors_and_nonneighbors(u, R_end, neighbors_n, nonneighbors_n);
        assert(neighbors_n + nonneighbors_n == R_end-1);
        for(ui i = 0;i < neighbors_n;i ++) ++ degree_in_S[neighbors[i]];
        while(!Qv.empty()) Qv.pop();
        for(ui i = 0;i < nonneighbors_n;i ++) {
            ui v = nonneighbors[i];
            if(SR_rid[v] >= S_end) {
                if(level_id[v] == level) continue;
                if(S_end - degree_in_S[v] >= K||S_end - degree_in_S[u] == K) {
                    level_id[v] = level;
                    Qv.push(v);
                }
            }
            else if(S_end - degree_in_S[v] == K) {
                long long nv = (long long) v*n;
                char *tt_matrix = matrix + nv;
                for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!tt_matrix[SR[j]]) {
                        level_id[SR[j]] = level;
                        Qv.push(SR[j]);
                    }
            }
        }
        return remove_vertices_and_edges_with_prune(S_end, R_end, R_w, level);
    }

    bool is_kplex(ui R_end) {
        for(ui i = 0;i < R_end;i ++) if(degree[SR[i]] + K < R_end) return false;
        return true;
    }

    void BB_search(ui S_end, ui R_end, ui S_w, ui R_w, ui level, ui choose_zero) {
        ui old_S_end = S_end, old_R_end = R_end;

        if(choose_zero&&SR_rid[0] < R_end&&!move_u_to_S_with_prune(0, S_end, R_end, S_w, R_w, level)) {
            restore_SR_and_edges(S_end, R_end, R_w, S_w, old_S_end, old_R_end, level);
            return ;
        }

        if(R_end < q) {return;}
        ui m_vw = 0;
        for(ui i = S_end;i < R_end;i ++) if(weight_now[SR[i]] > m_vw) m_vw = weight_now[SR[i]];
        if(m_vw * (R_end - S_end) + S_w <= best_solution_weight) {return;}

        if(!upper_1(S_end, R_end, S_w)) {return;}
        if(!upper_3(S_end, R_end, S_w)) {return;}
        if(S_w > best_solution_weight && S_end >= 2 * K - 1 ) store_solution(S_end, S_w);
        if(R_w > best_solution_weight&&is_kplex(R_end)) store_solution(R_end, R_w);
        if(R_w <= best_solution_weight) return ;

        ui u = choose_branch_vertex(S_end, R_end, R_w);
        ui  t_old_R_end = R_end;


        assert(Qv.empty());
        bool f1 = upper_2(u, S_end, R_end, S_w);
        bool f = move_u_to_S_with_prune(u, S_end, R_end, S_w, R_w, level);
        if(f && f1){
            BB_search(S_end, R_end, S_w, R_w, level+1, false);
        }
        restore_SR_and_edges(S_end, R_end, R_w, S_w, S_end, t_old_R_end, level);
        bool succeed = remove_u_from_S_with_prune(S_end, R_end, R_w, S_w, level);
        if(succeed) succeed = remove_vertices_and_edges_with_prune(S_end, R_end, R_w, level);
        if(succeed) {
            BB_search(S_end, R_end, S_w, R_w, level+1, false);
        }
        restore_SR_and_edges(S_end, R_end, R_w, S_w, old_S_end, old_R_end, level);
    }

    bool upper_1(ui S_end, ui R_end, ui S_w) {
        std::priority_queue<int, std::vector<int>, std::greater<int>> minHeap;
        ui * l_c = new ui[mc];
        for(ui i = 0;i < S_end;i ++){
            ui u = SR[i];
            long long un = (long long) u*n;
            char *t_matrix = matrix + un;
            ui non_neighbor = S_end - degree_in_S[u];
            ui ub = S_w;
            for(ui j= 0;j < mc;j ++) l_c[j] = K;
            for (ui j = 0; j < S_end; j++) {
                l_c[color_now[SR[j]]]--;
            }
            for(ui j = S_end;j < R_end;j ++)if(l_c[color_now[SR[j]]]){
                    if(t_matrix[SR[j]]) {
                        ub += weight_now[SR[j]];
                        if(ub > best_solution_weight){
                            delete []l_c;
                            return true;
                        }
                    }
                    else if(non_neighbor < K){
                        if (minHeap.size() < K - non_neighbor) {
                            minHeap.push(weight_now[SR[j]]);
                        } else if (weight_now[SR[j]] > minHeap.top()) {
                            minHeap.pop();
                            minHeap.push(weight_now[SR[j]]);
                        }
                    }
                }
            while (!minHeap.empty()) {
                ub += minHeap.top();
                minHeap.pop();
                if(ub > best_solution_weight){
                    while(!minHeap.empty()) minHeap.pop();
                    delete []l_c;
                    return true;
                }
            }
            if(ub <= best_solution_weight){
                delete []l_c;
                return false;
            }
        }
        delete[]l_c;
        return true;
    }

    bool upper_2(ui u, ui S_end, ui R_end, ui S_w) {
        ui * l_c = new ui[mc];
        for(ui i = 0;i < mc;i ++) l_c[i] = K;
        for(ui i = 0;i < S_end;i ++) --l_c[color_now[SR[i]]];
        l_c[color_now[u]]--;
        ui ub = S_w + weight_now[u];
        long long un = (long long) u*n;
        char *t_matrix = matrix + un;
        ui non_n = S_end - degree_in_S[u] + 1;
        priority_queue<int, std::vector<int>, std::greater<int>> minHeap;
        priority_queue<int>minHeap_n;
        unordered_map<int, std::vector<ui>> groupedData;
        for (ui i = S_end;i < R_end;i ++) if(l_c[color_now[SR[i]]]){
                if(t_matrix[SR[i]]){
                    groupedData[S_end - degree_in_S[SR[i]]].push_back(SR[i]);
                    minHeap_n.push(weight_now[SR[i]]);
                }
                else if (non_n < K) {
                    if (minHeap.size() < K - non_n) {
                        minHeap.push(weight_now[SR[i]]);
                    }
                    else if (weight_now[SR[i]] > minHeap.top()) {
                        minHeap.pop(); 
                        minHeap.push(weight_now[SR[i]]); 
                    }
                }
            }
        while (!minHeap.empty()) {
            ub += minHeap.top();
            if(ub > best_solution_weight) {
                while (!minHeap.empty()) minHeap.pop();
                while(!minHeap_n.empty()) minHeap_n.pop();
                delete []l_c;
                return true;
            }
            minHeap.pop();
        }
        ui *s = new ui[n];
        ui sup = 0;
        for(ui i = 0; i < n;i ++) s[i] = 0;

        for(ui i = 0;i < S_end;i ++) {
            s[SR[i]] = K - (S_end - degree_in_S[SR[i]]);
            if (!t_matrix[SR[i]]) s[SR[i]]--;
            sup +=s[SR[i]];
        }
        for(ui i = 0;i <= K - 1;i ++){
            if(groupedData[i].empty()) continue;
            for(ui j = 0;j < groupedData[i].size();j ++){
                if(sup >= i){
                    ui v = groupedData[i][j];
                    if(!l_c[color_now[v]]) continue;
                    if(i == 0){
                        ub += minHeap_n.top(); sup = sup - i;
                        minHeap_n.pop();
                        if(ub > best_solution_weight){
                            delete []l_c;
                            delete []s;
                            while(!minHeap_n.empty()) minHeap_n.pop();
                            return true;
                        }
                    }
                    else{
                        long long nv = (long long) n*v;
                        char *v_matrix = matrix + nv;
                        int idx = -1, min = K;
                        for(ui z = 0;z < S_end;z ++) if(!v_matrix[SR[z]]) {
                                if(s[SR[z]] < min) {
                                    idx = SR[z];
                                    min = s[idx];
                                }
                            }
                        if(min > 0){
                            ub += minHeap_n.top(); sup = sup - i;
                            minHeap_n.pop(); s[idx]--;
                            if(ub > best_solution_weight){
                                delete []l_c;
                                delete []s;
                                while(!minHeap_n.empty()) minHeap_n.pop();
                                return true;
                            }
                        }
                    }
                }
            }
        }
        delete []l_c;
        delete []s;
        while(!minHeap_n.empty()) minHeap_n.pop();
        if(ub > best_solution_weight) return true;
        return false;
    }

    bool upper_3(ui S_end, ui R_end, ui S_w) {
        unordered_map<ui, vector<ui>> colorTopWeights;
        ui * l_c = new ui[mc];
        for(ui i = 0;i < mc;i ++) l_c[i] = K;
        for (ui i = 0; i < S_end; i++) {
            l_c[color_now[SR[i]]]--;
        }
        for (int i = S_end; i < R_end; i++) {
            ui u = SR[i];
            ui colorId = color_now[u];
            ui nodeWeight = weight_now[u];
            auto& topWeights  = colorTopWeights[colorId];
            if (topWeights.size() < l_c[colorId]) {
                topWeights.push_back(nodeWeight);
                std::push_heap(topWeights.begin(), topWeights.end(), std::greater<ui>());
            } else if (nodeWeight > topWeights.front()) {
                std::pop_heap(topWeights.begin(), topWeights.end(), std::greater<ui>());
                topWeights.back() = nodeWeight;
                std::push_heap(topWeights.begin(), topWeights.end(), std::greater<ui>());
            }
        }
        ui totalWeightSum = S_w;
        for (const auto& pair : colorTopWeights) {
            const auto& topWeights = pair.second;
            for (ui weight : topWeights) {
                totalWeightSum += weight;
            }
        }
        delete[] l_c;
        colorTopWeights.clear();
        if(totalWeightSum <= best_solution_weight) return false;
        return true;

    }

    void get_neighbors_and_nonneighbors(ui u, ui R_end, ui &neighbors_n, ui &nonneighbors_n) {
        neighbors_n = 0; nonneighbors_n = 0;
        long long un = (long long) u*n;
        char *t_matrix = matrix + un;
        for(ui i = 0;i < R_end;i ++) if(SR[i] != u) {
                if(t_matrix[SR[i]] ) neighbors[neighbors_n++] = SR[i];
                else nonneighbors[nonneighbors_n++] = SR[i];
            }
    }

    bool move_u_to_S_with_prune(ui u, ui &S_end, ui &R_end, ui &S_w, ui &R_w, ui level) {
        assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end&&SR[SR_rid[u]] == u);
        assert(degree_in_S[u] + K > S_end);
        if(SR_rid[u] != S_end) swap_pos(S_end, SR_rid[u]);
        ++ S_end;
        S_w +=weight_now[u];
        int * l_c = new int[mc];
        for(ui i = 0;i < mc;i ++) l_c[i] = K;
        for(ui i = 0;i < S_end;i ++) --l_c[color_now[SR[i]]];

        ui neighbors_n = 0, nonneighbors_n = 0;
        get_neighbors_and_nonneighbors(u, R_end, neighbors_n, nonneighbors_n);
        assert(neighbors_n + nonneighbors_n == R_end-1);
        for(ui i = 0;i < neighbors_n;i ++) ++ degree_in_S[neighbors[i]];
        while(!Qv.empty()) Qv.pop();
        for(ui i = S_end;i < R_end;i ++) {
            if(l_c[color_now[SR[i]]] <= 0){
                level_id[SR[i]] = level;
                Qv.push(SR[i]);
            }
        }

        for(ui i = 0;i < nonneighbors_n;i ++) {
            ui v = nonneighbors[i];
            if(SR_rid[v] >= S_end) {
                if(level_id[v] == level) continue;
                if(S_end - degree_in_S[v] >= K||S_end - degree_in_S[u] == K) {
                    level_id[v] = level;
                    Qv.push(v);
                }
            }
            else if(S_end - degree_in_S[v] == K) {
                long long nv = (long long) v*n;
                char *tt_matrix = matrix + nv;
                for(ui j = S_end;j < R_end;j ++) if(level_id[SR[j]] > level&&!tt_matrix[SR[j]]) {
                        level_id[SR[j]] = level;
                        Qv.push(SR[j]);
                    }
            }
        }
        delete []l_c;
        return remove_vertices_and_edges_with_prune(S_end, R_end, R_w, level);
    }

    bool remove_vertices_and_edges_with_prune(ui S_end, ui &R_end, ui &R_w, ui level) {
        while(!Qv.empty()) {
            while(!Qv.empty()) {
                ui u = Qv.front(); Qv.pop(); 
                assert(SR[SR_rid[u]] == u);
                assert(SR_rid[u] >= S_end&&SR_rid[u] < R_end);
                -- R_end;
                R_w -=weight_now[u];
                swap_pos(SR_rid[u], R_end);
                bool terminate = false;
                ui neighbors_n = 0;
                long long un = (long long) u*n;
                char *t_matrix = matrix + un;
                for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
                        ui w = SR[i];
                        neighbors[neighbors_n++] = w;
                        -- degree[w];
                        if(degree[w] < q - K) {
                            if(i < S_end) terminate = true; 
                            else if(level_id[w] > level) { 
                                level_id[w] = level;
                                Qv.push(w);
                            }
                        }
                    }

                if(terminate) {
                    for(ui i = 0;i < neighbors_n;i ++) ++ degree[neighbors[i]];
                    level_id[u] = n;
                    ++ R_end;
                    R_w +=weight_now[u];
                    return false;
                }
            }
        }
        return true;
    }

    void restore_SR_and_edges(ui &S_end, ui &R_end, ui &R_w, ui &S_w, ui old_S_end, ui old_R_end, ui level) {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(level_id[u] == level);
            assert(SR_rid[u] < R_end);
            level_id[u] = n;
        }

        for(;R_end < old_R_end;R_end ++) { 
            ui u = SR[R_end];
            R_w +=weight_now[u];
            assert(level_id[u] == level&&SR_rid[u] == R_end);
            level_id[u] = n;

            ui neighbors_n = 0;
            long long un = (long long) u*n;
            char *t_matrix = matrix + un;
            degree[u] = degree_in_S[u] = 0;
            for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
                    ui w = SR[i];
                    neighbors[neighbors_n ++] = w;
                    ++ degree[w];
                    ++ degree[u];
                    if(i < S_end) ++ degree_in_S[u];
                }
        }
        for(;S_end > old_S_end;S_end --) { 
            ui u = SR[S_end-1];
            S_w -= weight_now[u];
            assert(SR_rid[u] == S_end-1);

            ui neighbors_n = 0;
            long long un = (long long) u*n;
            char *t_matrix = matrix + un;
            for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) {
                    ui w = SR[i];
                    neighbors[neighbors_n ++] = w;
                    -- degree_in_S[w];
                }
        }
    }

    bool remove_u_from_S_with_prune(ui &S_end, ui &R_end, ui &R_w, ui &S_w, ui level) {
        assert(S_end);
        ui u = SR[S_end-1];
        -- S_end; -- R_end;
        S_w -=weight_now[u];
        R_w -=weight_now[u];
        swap_pos(S_end, R_end);
        level_id[u] = level;

        bool terminate = false;
        ui neighbors_n = 0;
        long long un = (long long) u*n;
        char *t_matrix = matrix + un;
        for(ui i = 0;i < R_end;i ++) if(t_matrix[SR[i]]) neighbors[neighbors_n ++] = SR[i];
        for(ui i = 0;i < neighbors_n;i ++) {
            ui v = neighbors[i];
            -- degree_in_S[v];
            -- degree[v];
            if(degree[v] < q - K) {
                if(SR_rid[v] < S_end) terminate = true;
                else {
                    assert(level_id[v] > level);
                    level_id[v] = level;
                    Qv.push(v);
                }
            }
        }
        if(terminate) return false;
        return true;
    }

    void swap_pos(ui i, ui j) {
        std::swap(SR[i], SR[j]);
        SR_rid[SR[i]] = i;
        SR_rid[SR[j]] = j;
    }

    ui choose_branch_vertex(ui S_end, ui R_end, ui R_w) {
        assert(S_end < SR.size());
        assert(R_end <= SR.size());

        std::vector<ui> D;
        D.reserve(R_end - S_end);  

        for (ui i = S_end; i < R_end; i++) {
            if (R_end - degree[SR[i]] > K) {
                D.push_back(SR[i]);
            }
        }

        if (D.empty()) {
            return SR[S_end];
        }

        ui u = D[0];
        ui min_degree_in_S = degree_in_S[u];
        ui min_color = color_cnt[u];
        ui min_weight = weight_now[u];
        ui min_degree = degree[u];

        for (size_t i = 1; i < D.size(); ++i) {
            ui v = D[i];
            if (degree_in_S[v] < min_degree_in_S ||
                (degree_in_S[v] == min_degree_in_S && color_cnt[v] < min_color) ||
                (degree_in_S[v] == min_degree_in_S && color_cnt[v] == min_color && weight_now[v] < min_weight) ||
                (degree_in_S[v] == min_degree_in_S && color_cnt[v] == min_color && weight_now[v] == min_weight && degree[v] < min_degree)) {
                u = v;
                min_degree_in_S = degree_in_S[v];
                min_color = color_cnt[v];
                min_weight = weight_now[v];
                min_degree = degree[v];
            }
        }
        return u;
    }

};
#endif









