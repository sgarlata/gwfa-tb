#ifndef DP_H
#define DP_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
using namespace std;

//// DP MATRIX CELL TYPE
typedef struct gwf_cigar_t
{
    int32_t s;   //// edit distance
    char *op;    //// edits array (TODO: understand whether it would be better to convert this to 'string')
    int32_t *bl; //// edits number (TODO: understand whether it would be better to convert this convert to 'vector')
    int32_t l;   ////  array length
} gwf_cigar_t;

//// DIAGONAL-OFFSET -> ROW-COLUMN
inline int32_t get_row(vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v, int32_t d)
{
    return diag_row_map[v][d];
}

inline int32_t get_col(vector<vector<int32_t>> &diag_off, int32_t v, int32_t r_dp, int32_t c_dp, int32_t r)
{
    return min(r_dp, c_dp) - diag_off[v][r];
}

//// Check if DP matrix cell is empty or actually stored in the data structure
int32_t
dp_check_cell(unordered_map<int32_t, int32_t> &v_map, vector<vector<vector<gwf_cigar_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, vector<vector<int32_t>> &diag_off, int32_t vtx, int32_t d, int32_t k)
{
    int32_t v, r, c;
    if (v_map.count(vtx)) //// vertex: ok
    {
        v = v_map[vtx];
        if (diag_row_map[v].count(d)) //// diagonal ($dpd's row): ok
        {
            r = get_row(diag_row_map, v, d);
            c = get_col(diag_off, v, d + k, k, r);

            if (0 <= c && c < (int32_t)dpd[v][r].size()) //// column: ok
            {
                return 1;
            }
        }
    }
    return 0;
}

//// DP matrix extension (within same vertex)
void dp_extend(vector<vector<vector<gwf_cigar_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, vector<vector<int32_t>> &diag_off, int32_t v_dp, int32_t v, int32_t d, int32_t prev_k, int32_t k, FILE *out_debug)
{
    int32_t r, c, r_dp, c_dp, r_prev, c_prev;

    for (c_dp = prev_k; c_dp <= k; ++c_dp) //// along the previous cells of the diagonal to extend //// TODO: check that starting from prev_k is right (before it started all the way from 0)
    {
        r_dp = d + c_dp;            //// row in the traditional dpd matrix
        if (r_dp >= 0 && c_dp >= 0) //// within bounds
        {
            r = get_row(diag_row_map, v, d);
            c = get_col(diag_off, v, r_dp, c_dp, r);

            if (v == 0 && d == 0 && c == 0 && dpd[v][r][c].s == INT32_MAX) //// dpd[0][0][0]: first match
            {
                dpd[v][r][c].s = 0;
#ifdef DP_DEBUG
                fprintf(out_debug, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif

                dpd[v][r][c].op = (char *)malloc(sizeof(char));
                dpd[v][r][c].bl = (int32_t *)malloc(sizeof(int32_t));
                if (dpd[v][r][c].op == NULL || dpd[v][r][c].bl == NULL)
                {
                    fprintf(stderr, "Allocation error at line number %d in file %s\n", __LINE__, __FILE__);
                    abort();
                }
                dpd[v][r][c].op[0] = '=';
                dpd[v][r][c].bl[0] = 1;
                dpd[v][r][c].l++; //// l is first used as index, while from now on as length (+1)
            }
            else if (c > 0) //// central cells
            {               //// Previous condition: "r_dp > 0 && c_dp > diag_off[v][r]"
                /* if ((int32_t)dpd[v][r].size() <= c) //// check if the column has still to be allocated
                {
                    dpd[v][r].push_back({.s = -1, .op = NULL, .bl = NULL, .l = 0});
                } */

                //// TODO: understand if here the for is needed
                int32_t c_start = (int32_t)dpd[v][r].size() - 1;
                for (int32_t i = 0; i < (c - c_start); ++i) //// add as many columns up to $c, if empty they will be later filled during $d_new's extension
                    dpd[v][r].push_back({.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0});

                if ((prev_k < c_dp) && (dpd[v][r][c].s == INT32_MAX || dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s))
                { //// match

                    if (dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s) //// update
                    {
                        free(dpd[v][r][c].op);
                        free(dpd[v][r][c].bl);
                    }

                    dpd[v][r][c].s = dpd[v][r][c - 1].s;

                    if (dpd[v][r][c - 1].op[dpd[v][r][c - 1].l - 1] == '=') //// if already coming from a match
                    {
                        dpd[v][r][c].l = dpd[v][r][c - 1].l - 1; //// keep same length, but decrese it for the moment to use it as index
                    }
                    else
                    {
                        dpd[v][r][c].l = dpd[v][r][c - 1].l; //// previous length + 1 to store the additional operation
                    }
                    dpd[v][r][c].op = (char *)malloc((dpd[v][r][c].l + 1) * sizeof(char)); //// + 1 as l still stands for the index at the moment
                    dpd[v][r][c].bl = (int32_t *)malloc((dpd[v][r][c].l + 1) * sizeof(int32_t));
                    if (dpd[v][r][c].op == NULL || dpd[v][r][c].bl == NULL)
                    {
                        fprintf(stderr, "Allocation error at line number %d in file %s\n", __LINE__, __FILE__);
                        abort();
                    }

                    for (int32_t i = 0; i < dpd[v][r][c - 1].l; ++i) //// copy
                    {
                        dpd[v][r][c].op[i] = dpd[v][r][c - 1].op[i];
                        dpd[v][r][c].bl[i] = dpd[v][r][c - 1].bl[i];
                    }

                    if (dpd[v][r][c].l == dpd[v][r][c - 1].l - 1) //// kept same length
                    {
                        dpd[v][r][c].bl[dpd[v][r][c].l]++; //// just increment
                    }
                    else
                    {
                        dpd[v][r][c].op[dpd[v][r][c].l] = '=';
                        dpd[v][r][c].bl[dpd[v][r][c].l] = 1;
                    }

                    dpd[v][r][c].l++; //// now l stands for the length

#ifdef DP_DEBUG
                    fprintf(out_debug, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_dp, r_dp - 1, c_dp - 1, dpd[v][r][c - 1].s, v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
                }
                else if ((prev_k == c_dp) && (dpd[v][r][c].s == INT32_MAX || dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s + 1)) //// TODO: make it >= if you want to give precedence to mismatches
                {                                                                                                                                        //// mismatch
                    if (dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s + 1)                                                       //// update
                    {
                        free(dpd[v][r][c].op);
                        free(dpd[v][r][c].bl);
                    }
                    dpd[v][r][c].s = dpd[v][r][c - 1].s + 1;

                    if (dpd[v][r][c - 1].op[dpd[v][r][c - 1].l - 1] == 'X') //// if already coming from a mismatch
                    {
                        dpd[v][r][c].l = dpd[v][r][c - 1].l - 1; //// keep same length, but decrese it for the moment to use it as index
                    }
                    else
                    {
                        dpd[v][r][c].l = dpd[v][r][c - 1].l; //// previous length + 1 to store the additional operation
                    }
                    dpd[v][r][c].op = (char *)malloc((dpd[v][r][c].l + 1) * sizeof(char)); //// + 1 as l still stands for the index at the moment
                    dpd[v][r][c].bl = (int32_t *)malloc((dpd[v][r][c].l + 1) * sizeof(int32_t));
                    if (dpd[v][r][c].op == NULL || dpd[v][r][c].bl == NULL)
                    {
                        fprintf(stderr, "Allocation error at line number %d in file %s\n", __LINE__, __FILE__);
                        abort();
                    }

                    for (int32_t i = 0; i < dpd[v][r][c - 1].l; ++i) //// copy
                    {
                        dpd[v][r][c].op[i] = dpd[v][r][c - 1].op[i];
                        dpd[v][r][c].bl[i] = dpd[v][r][c - 1].bl[i];
                    }

                    if (dpd[v][r][c].l == dpd[v][r][c - 1].l - 1) //// kept same length
                    {
                        dpd[v][r][c].bl[dpd[v][r][c].l]++; //// just increment
                    }
                    else
                    {
                        dpd[v][r][c].op[dpd[v][r][c].l] = 'X';
                        dpd[v][r][c].bl[dpd[v][r][c].l] = 1;
                    }

                    dpd[v][r][c].l++; //// now l stands for the length
#ifdef DP_DEBUG
                    fprintf(out_debug, "[DEBUG] Extension (X): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_dp, r_dp - 1, c_dp - 1, dpd[v][r][c - 1].s, v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
                }
            }
        }
    }
}

//// DP matrix expansion (within same vertex)
void dp_expand(vector<vector<vector<gwf_cigar_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, vector<vector<int32_t>> &diag_off, int32_t v_dp, int32_t v, int32_t d, int32_t k, char ed, FILE *out_debug)
{
    int32_t r, r_old, r_dp, c, c_old, c_dp, d_new, off = 0;

    if (ed == 'D') //// deletion
    {
        d_new = d - 1;
        r_dp = k + d;
        c_dp = k + 1;
        if (r_dp < 0 || c_dp <= 0)
            return;
        off = min(r_dp, c_dp);
    }
    else if (ed == 'X') //// mismatch
    {
        d_new = d;
        r_dp = k + d + 1;
        c_dp = k + 1;
        if (r_dp < 0 || c_dp < 0)
            return;
    }
    else if (ed == 'I') //// insertion
    {
        d_new = d + 1;
        r_dp = k + d + 1;
        c_dp = k;
        if (r_dp <= 0 || c_dp < 0)
            return;
        off = min(r_dp, c_dp);
    }

    if (ed != 'X' && diag_row_map[v].count(d_new) == 0) //// new diagonal
    {
        dpd[v].push_back(vector<gwf_cigar_t>(1, {.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0})); //// add row to dpd[v]
        diag_row_map[v].insert({d_new, ((int32_t)dpd[v].size() - 1)});                              //// add mapping to diag_row_map[v]
        diag_off[v].push_back(off);                                                                 //// add row's offset to diag_off[v]
    }
    r_old = get_row(diag_row_map, v, d);           //// row of the diagonal we are coming from
    c_old = get_col(diag_off, v, d + k, k, r_old); //// column of the diagonal we are coming from (possibly decreased by its row's offset)
    r = get_row(diag_row_map, v, d_new);
    c = get_col(diag_off, v, r_dp, c_dp, r);

    int32_t c_start = (int32_t)dpd[v][r].size() - 1;
    for (int32_t i = 0; i < (c - c_start); ++i) //// add as many columns up to $c, if empty they will be later filled during $d_new's extension
        dpd[v][r].push_back({.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0});

    //// mismatch only: mismatch at very first base comparison
    if (ed == 'X' && v == 0 && d == 0 && k == -1)
    {
        dpd[v][r][c].s = 1; //// Right now, s == 0

        dpd[v][r][c].op = (char *)malloc(sizeof(char));
        dpd[v][r][c].bl = (int32_t *)malloc(sizeof(int32_t));
        if (dpd[v][r][c].op == NULL || dpd[v][r][c].bl == NULL)
        {
            fprintf(stderr, "Allocation error at line number %d in file %s\n", __LINE__, __FILE__);
            abort();
        }
        dpd[v][r][c].op[dpd[v][r][c].l] = 'X'; //// Here, dpd[v][r][c].l = 0
        dpd[v][r][c].bl[dpd[v][r][c].l] = 1;
        dpd[v][r][c].l++; //// l is first used as index, while from now on as length (+1)

#ifdef DP_DEBUG
        fprintf(out_debug, "[DEBUG] Starting mismatch (X): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
    }
    else if (dpd[v][r][c].s == INT32_MAX || dpd[v][r][c].s > dpd[v][r_old][c_old].s + 1)
    {
        if (dpd[v][r][c].s > dpd[v][r_old][c_old].s + 1) //// update
        {
            free(dpd[v][r][c].op);
            free(dpd[v][r][c].bl);
        }

        dpd[v][r][c].s = dpd[v][r_old][c_old].s + 1;

        if (dpd[v][r_old][c_old].op[dpd[v][r_old][c_old].l - 1] == ed) //// if same edit
        {
            dpd[v][r][c].l = dpd[v][r_old][c_old].l - 1; //// keep same length, but decrease it for the moment to use it as index
        }
        else
        {
            dpd[v][r][c].l = dpd[v][r_old][c_old].l; //// previous length + 1 to store the additional operation
        }
        dpd[v][r][c].op = (char *)malloc((dpd[v][r][c].l + 1) * sizeof(char)); //// + 1 as l still stands for the index at the moment
        dpd[v][r][c].bl = (int32_t *)malloc((dpd[v][r][c].l + 1) * sizeof(int32_t));
        if (dpd[v][r][c].op == NULL || dpd[v][r][c].bl == NULL)
        {
            fprintf(stderr, "Allocation error at line number %d in file %s\n", __LINE__, __FILE__);
            abort();
        }

        for (int32_t i = 0; i < dpd[v][r_old][c_old].l; ++i) //// copy
        {
            dpd[v][r][c].op[i] = dpd[v][r_old][c_old].op[i];
            dpd[v][r][c].bl[i] = dpd[v][r_old][c_old].bl[i];
        }

        if (dpd[v][r][c].l == dpd[v][r_old][c_old].l - 1) //// kept same length
        {
            dpd[v][r][c].bl[dpd[v][r][c].l]++; //// just increment
        }
        else
        {
            dpd[v][r][c].op[dpd[v][r][c].l] = ed;
            dpd[v][r][c].bl[dpd[v][r][c].l] = 1;
        }

        dpd[v][r][c].l++; //// now l stands for the length

#ifdef DP_DEBUG
        fprintf(out_debug, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", ed, v_dp, d + k, k, dpd[v][r_old][c_old].s, v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
    }
}

//// DP matrix data structures setup when a new vertex is visited
void dp_new_vd(unordered_map<int32_t, int32_t> &v_map, vector<vector<vector<gwf_cigar_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, vector<vector<int32_t>> &diag_off, int32_t w, int32_t d, int32_t r_dp, int32_t c_dp, int32_t &r, int32_t &c)
{
    if (v_map.count(w) == 0) //// check whether the DP already contains the vertex
    {
        v_map[w] = v_map.size();
        dpd.push_back(vector<vector<gwf_cigar_t>>(1, vector<gwf_cigar_t>(1, {.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0})));
        diag_row_map.push_back({{d, 0}});          //// new vertex, with diagonal $d to row 0
        diag_off.push_back(vector<int32_t>(1, 0)); //// new vertex, with offset 0 (as we are directly extending to the vertex label) for its first row (0)
    }
    else if (diag_row_map[v_map[w]].count(d) == 0) //// new diagonal
    {
        dpd[v_map[w]].push_back(vector<gwf_cigar_t>(1, {.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0})); //// add row to dpd[v]
        diag_row_map[v_map[w]].insert({d, ((int32_t)dpd[v_map[w]].size() - 1)});                           //// add mapping to diag_row_map[v]
        diag_off[v_map[w]].push_back(0);
    }
    r = get_row(diag_row_map, v_map[w], d);
    c = get_col(diag_off, v_map[w], r_dp, c_dp, r); //// TODO: problem here with negative $c, maybe offset needs be reconsidered, also possible to define a diagonal DP type

    /* if ((int32_t)dpd[v_map[w]][r].size() <= c) //// check whether the column is already there or not
    {
#ifdef DP_DEBUG_COL
        if ((int32_t)dpd[v_map[w]][r].size() < c)
            fprintf(stdout, "[DEBUG] Column (dpd[%d][%d].size() = %d, c = %d\n", v_map[w], r, (int32_t)dpd[v_map[w]][r].size(), c);
#endif
        dpd[v_map[w]][r].push_back({.s = -1, .op = NULL, .bl = NULL, .l = 0});
    } */

    for (int32_t i = (int32_t)dpd[v_map[w]][r].size() - 1; i <= c; ++i)
        dpd[v_map[w]][r].push_back({.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0});
}

//// print cigar string
void gwf_cigar(gwf_cigar_t cig)
{
    fprintf(stdout, "CIGAR:\t");
    for (int32_t i = 0; i < cig.l; ++i)
    {
        fprintf(stdout, "%d%c", cig.bl[i], cig.op[i]);
    }
    fprintf(stdout, "\n");
}

#endif