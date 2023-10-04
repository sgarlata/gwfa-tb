#ifndef TB_H
#define TB_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
using namespace std;

//// DIAGONAL TYPE FOR TRACEBACK
typedef struct tb_diag_t
{
    int32_t s;          //// edit distance
    vector<char> op;    //// edits array
    vector<int32_t> bl; //// edits number
    int32_t off;        //// diagonal offset
} tb_diag_t;

//// DIAGONAL-OFFSET -> ROW-COLUMN
inline int32_t
get_row(vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v, int32_t d)
{
    return diag_row_map[v][d];
}

//// Discard last edit
void discard_ed(tb_diag_t &wfd)
{
    if (wfd.bl.back() > 0)
    {
        --wfd.bl.back();
    }
    else
    {
        wfd.op.pop_back();
        wfd.bl.pop_back();
    }
}

//// For DEBUG purposes: print current tracebacks
void print_tb(vector<vector<tb_diag_t>> wf)
{
    for (int32_t i = 0; i < (int32_t)wf.size(); i++)
    {
        fprintf(stdout, "\tNode %d\n", i);
        for (int32_t j = 0; j < (int32_t)wf[i].size(); j++)
        {
            fprintf(stdout, "\t\tDiagonal %d: ", j);
            for (int32_t k = 0; k < (int32_t)wf[i][j].op.size(); ++k)
            {
                fprintf(stdout, "%d%c", wf[i][j].bl[k], wf[i][j].op[k]);
            }
            fprintf(stdout, "\n");
        }
    }
}

//// Extension for traceback (within same vertex)
void tb_extend(int32_t s, vector<vector<tb_diag_t>> &wf, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v_dp, int32_t v, int32_t d, int32_t prev_k, int32_t k, FILE *out_debug)
{
    int32_t r, r_dp, c_dp;

    for (c_dp = prev_k; c_dp <= k; ++c_dp)
    {
        r_dp = d + c_dp;            //// row in the traditional dpd matrix
        if (r_dp >= 0 && c_dp >= 0) //// within bounds
        {
            r = get_row(diag_row_map, v, d);

            if (v == 0 && d == 0 && wf[v][r].s == INT32_MAX) //// MATCH OF FIRST ELEMENT
            {
                wf[v][r].s = 0;
                wf[v][r].op.push_back('=');
                wf[v][r].bl.push_back(1);
                //// offset already set to 0
#ifdef TB_DEBUG
                fprintf(stdout, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, wf[v][r].s);
#endif

#ifdef TB_PRINT
                print_tb(wf);
#endif
            }
            //// MATCH OF OTHER ELEMENTS
            else if ((prev_k < c_dp) && (wf[v][r].off <= c_dp || wf[v][r].s > s))
            {
                /* if (wf[v][r].s > s) //// UPDATE -> discard previous traceback and impose incoming one
                {
                    wf[v][r].s = s;
                    //// there's no way to recover the previous traceback, as it would have been already overwritten
                    //// so, an idea could be to remove backwards one by one the new edits
                    for (int32_t i = 0; i < wf[v][r].off - c_dp; ++i)
                    {
                        if (wf[v][r].bl.back() > 1)
                        {
                            --wf[v][r].bl.back();
                        }
                        else
                        {
                            wf[v][r].op.pop_back();
                            wf[v][r].bl.pop_back();
                        }
                    }
                } */

                if (!wf[v][r].op.empty() && wf[v][r].op.back() == '=') //// if already coming from a match
                {
                    wf[v][r].bl.back()++; //// just increase counter
                }
                else //// else, add new match
                {
                    wf[v][r].op.push_back('=');
                    wf[v][r].bl.push_back(1);
                }

                wf[v][r].off = c_dp;

#ifdef TB_DEBUG
                fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_dp, r_dp - 1, c_dp - 1, wf[v][r].s, v_dp, r_dp, c_dp, wf[v][r].s);
#endif

#ifdef TB_PRINT
                print_tb(wf);
#endif
            }
        }
    }
}

//// Expansion for traceback (within same vertex)
void tb_expand(int32_t s, vector<vector<tb_diag_t>> &wf, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v_dp, int32_t v, int32_t v_len, int32_t d, int32_t k, char ed, FILE *out_debug)
{
    int32_t r, r_from, r_dp, c_dp, d_to;

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
        //// offset already set to 0

#ifdef TB_DEBUG
        fprintf(stdout, "[DEBUG] Starting mismatch (X): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, wf[v][r].s);
#endif

#ifdef TB_PRINT
        print_tb(wf);
#endif
        return;
    }

    //// MISMATCH
    if (ed == 'X')
    {
        r = get_row(diag_row_map, v, d_to);
        //// not visited yet
        if (wf[v][r].off <= c_dp)
        {
            wf[v][r].s = s + 1;

            if (!wf[v][r].op.empty() && wf[v][r].op.back() == ed) //// if already coming from same edit
            {
                wf[v][r].bl.back()++; //// just increase counter
            }
            else //// else, add new edit
            {
                wf[v][r].op.push_back(ed);
                wf[v][r].bl.push_back(1);
            }

            wf[v][r].off = c_dp;

#ifdef TB_DEBUG
            fprintf(stdout, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", ed, v_dp, d + k, k, s, v_dp, r_dp, c_dp, s + 1);
#endif

#ifdef TB_PRINT
            print_tb(wf);
#endif

            return;
        }
    }

    //// NEW DIAGONAL
    if (ed != 'X' && diag_row_map[v].count(d_to) == 0)
    {
        r_from = get_row(diag_row_map, v, d);
        tb_diag_t tmp_diag = wf[v][r_from];

        wf[v].push_back(tmp_diag);                                   //// add the diagonal to wf[v]
        diag_row_map[v].insert({d_to, ((int32_t)wf[v].size() - 1)}); //// add mapping to diag_row_map[v]

        r = get_row(diag_row_map, v, d_to);

        //// not visited yet
        if (wf[v][r].off <= c_dp)
        {
            wf[v][r].s = s + 1;

            if (!wf[v][r].op.empty() && wf[v][r].op.back() == ed) //// if already coming from same edit
            {
                wf[v][r].bl.back()++; //// just increase counter
            }
            else //// else, add new edit
            {
                wf[v][r].op.push_back(ed);
                wf[v][r].bl.push_back(1);
            }

            wf[v][r].off = c_dp;

#ifdef TB_DEBUG
            fprintf(stdout, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", ed, v_dp, d + k, k, s, v_dp, r_dp, c_dp, s + 1);
#endif

#ifdef TB_PRINT
            print_tb(wf);
#endif
        }
    }
}

//// DP matrix data structures setup when a new vertex is visited
void tb_new_vd(unordered_map<int32_t, int32_t> &v_map, vector<vector<tb_diag_t>> &wf, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t w, int32_t w_len, int32_t d, int32_t ol, int32_t r_dp, int32_t &r, tb_diag_t wf_from)
{
    if (v_map.count(w) == 0) //// visiting the vertex for the first time
    {
        v_map[w] = v_map.size();
        wf.push_back(vector<tb_diag_t>(1, {.s = wf_from.s, .op = wf_from.op, .bl = wf_from.bl, .off = ol}));
        diag_row_map.push_back({{d, 0}}); //// new vertex, with diagonal $d to row 0
        r = 0;
    }
    else if (diag_row_map[v_map[w]].count(d) == 0) //// vertex already visited, but new diagonal
    {
        wf[v_map[w]].push_back({wf_from.s, wf_from.op, wf_from.bl, ol}); //// add row to wf[v]
        r = ((int32_t)wf[v_map[w]].size() - 1);
        diag_row_map[v_map[w]].insert({d, r}); //// add mapping to diag_row_map[v]
    }
    else //// vertex and diagonal already there
    {
        r = get_row(diag_row_map, v_map[w], d);
    }
}

//// print cigar string
void gwf_cigar(tb_diag_t cig)
{
    fprintf(stdout, "CIGAR:\t");
    for (int32_t i = 0; i < (int32_t)cig.op.size(); ++i)
    {
        fprintf(stdout, "%d%c", cig.bl[i], cig.op[i]);
    }
    fprintf(stdout, "\n");
}

#endif