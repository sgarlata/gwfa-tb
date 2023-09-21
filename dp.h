#ifndef DP_H
#define DP_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
using namespace std;

//// DP MATRIX CELL TYPE
typedef struct dp_diag_t
{
    int32_t s;          //// edit distance
    vector<char> op;    //// edits array
    vector<int32_t> bl; //// edits number
    int32_t off;        //// diagonal offset
} dp_diag_t;

//// DIAGONAL-OFFSET -> ROW-COLUMN
inline int32_t
get_row(vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v, int32_t d)
{
    return diag_row_map[v][d];
}

//// DP matrix extension (within same vertex)
void dp_extend(int32_t s, vector<vector<dp_diag_t>> &wf, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v_dp, int32_t v, int32_t d, int32_t prev_k, int32_t k, FILE *out_debug)
{
    int32_t r, r_dp, c_dp;

    for (c_dp = prev_k; c_dp <= k; ++c_dp)
    {
        r_dp = d + c_dp;            //// row in the traditional dpd matrix
        if (r_dp >= 0 && c_dp >= 0) //// within bounds
        {
            r = get_row(diag_row_map, v, d);

            if (v == 0 && d == 0 && wf[v][r].s == INT32_MAX) //// STARTING MATCH
            {
                wf[v][r].s = 0;
                wf[v][r].op.push_back('=');
                wf[v][r].bl.push_back(1);
#ifdef DP_DEBUG
                fprintf(stdout, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, wf[v][r].s);
#endif
            }
            else //// central cells
            {
                //// if match and (first visit or update)
                if ((prev_k < c_dp) && (wf[v][r].off <= c_dp || wf[v][r].s > s))
                {

                    if (wf[v][r].s > s) //// update
                    {
#ifdef UPDATE_DEBUG
                        fprintf(stdout, "UPDATE: M\n"); //// this is likely the only update actually happening
#endif
                    }
#ifdef DP_DEBUG
                    fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_dp, r_dp - 1, c_dp - 1, wf[v][r].s, v_dp, r_dp, c_dp, s);
#endif

                    wf[v][r].s = s;

                    if (wf[v][r].op.back() == '=') //// if already coming from a match
                    {
                        wf[v][r].bl.back()++; //// just increase counter
                    }
                    else //// else, add new match
                    {
                        wf[v][r].op.push_back('=');
                        wf[v][r].bl.push_back(1);
                    }

                    wf[v][r].off++; //// increase the offset
                }
            }
        }
    }
}

//// DP matrix expansion (within same vertex)
void dp_expand(int32_t s, vector<vector<dp_diag_t>> &wf, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v_dp, int32_t v, int32_t v_len, int32_t d, int32_t k, char ed, FILE *out_debug)
{
    int32_t r, r_dp, c_dp, d_to;

    if (ed == 'D') //// deletion
    {
        d_to = d - 1;
        r_dp = k + d;
        c_dp = k + 1;
        if (r_dp < 0 || c_dp <= 0)
            return;
    }
    else if (ed == 'X') //// mismatch
    {
        d_to = d;
        r_dp = k + d + 1;
        c_dp = k + 1;
        if (r_dp < 0 || c_dp < 0)
            return;
    }
    else if (ed == 'I') //// insertion
    {
        d_to = d + 1;
        r_dp = k + d + 1;
        c_dp = k;
        if (r_dp <= 0 || c_dp < 0)
            return;
    }

    //// mismatch at very first base comparison
    if (ed == 'X' && v == 0 && d == 0 && k == -1)
    {
        r = 0;
        wf[v][r].s = 1; //// Right now, s == 0
        wf[v][r].op.push_back('X');
        wf[v][r].bl.push_back(1);
        wf[v][r].off++;

#ifdef DP_DEBUG
        fprintf(stdout, "[DEBUG] Starting mismatch (X): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, wf[v][r].s);
#endif
        return;
    }

    //// diagonal (row) setup
    if (ed != 'X' && diag_row_map[v].count(d_to) == 0) //// new diagonal
    {
        wf[v].push_back({INT32_MAX, {}, {}, c_dp - 1});              //// add the diagonal to wf[v]
        diag_row_map[v].insert({d_to, ((int32_t)wf[v].size() - 1)}); //// add mapping to diag_row_map[v]
    }

    //// coordinates
    r = get_row(diag_row_map, v, d_to);

    //// not visited yet
    if (wf[v][r].off < c_dp)
    {
        wf[v][r].s = s + 1;

        if (wf[v][r].op.empty() || wf[v][r].op.back() != ed) //// if empty or coming from different edit
        {
            wf[v][r].op.push_back(ed);
            wf[v][r].bl.push_back(1);
        }
        else if (wf[v][r].op.back() == ed) //// if same edit
        {
            wf[v][r].bl.back()++; //// increase length
        }

        wf[v][r].off++;

#ifdef DP_DEBUG
        fprintf(stdout, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", ed, v_dp, d + k, k, s, v_dp, r_dp, c_dp, s + 1);
#endif
    }
}

//// DP matrix data structures setup when a new vertex is visited
void dp_new_vd(unordered_map<int32_t, int32_t> &v_map, vector<vector<dp_diag_t>> &wf, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t w, int32_t w_len, int32_t d, int32_t ol, int32_t r_dp, int32_t &r)
{
    if (v_map.count(w) == 0) //// visiting the vertex for the first time
    {
        v_map[w] = v_map.size();
        wf.push_back(vector<dp_diag_t>(1, {.s = INT32_MAX, .op = {}, .bl = {}, .off = ol}));
        diag_row_map.push_back({{d, 0}}); //// new vertex, with diagonal $d to row 0
        r = 0;
    }
    else if (diag_row_map[v_map[w]].count(d) == 0) //// vertex already visited, but new diagonal
    {
        wf[v_map[w]].push_back({INT32_MAX, {}, {}, ol}); //// add row to wf[v]
        r = ((int32_t)wf[v_map[w]].size() - 1);
        diag_row_map[v_map[w]].insert({d, r}); //// add mapping to diag_row_map[v]
    }
    else //// vertex and diagonal already there
    {
        r = get_row(diag_row_map, v_map[w], d);
    }
}

//// print cigar string
void gwf_cigar(dp_diag_t cig)
{
    fprintf(stdout, "CIGAR:\t");
    for (int32_t i = 0; i < (int32_t)cig.op.size(); ++i)
    {
        fprintf(stdout, "%d%c", cig.bl[i], cig.op[i]);
    }
    fprintf(stdout, "\n");
}

#endif