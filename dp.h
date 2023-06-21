#ifndef DP_H
#define DP_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
using namespace std;

//// DP MATRIX CELL TYPE
typedef struct dp_cell_t
{
    int32_t s;   //// edit distance
    char *op;    //// edits array (TODO: understand whether it would be better to convert this to 'string')
    int32_t *bl; //// edits number (TODO: understand whether it would be better to convert this convert to 'vector')
    int32_t l;   ////  array length
} dp_cell_t;

////  DP DIAGONAL TYPE (TODO: refactor implementation to use this -> less DSs around)
typedef struct dp_diag_t
{
    int32_t v;
    int32_t d;
    vector<dp_cell_t> cell;
    int32_t row;
    int32_t off;
} dp_diag_t;

//// DIAGONAL-OFFSET -> ROW-COLUMN
inline int32_t
get_row(vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v, int32_t d)
{
    return diag_row_map[v][d];
}

inline int32_t get_col(int32_t r_dp, int32_t c_dp)
{
    /* int32_t col = min(r_dp, c_dp);
    int32_t off = diag_off[v][r];
    int32_t column = col - off;
    if (column < 0)
        fprintf(stderr, "[Negative column (%d)]: Vertex = %d, r_dp = %d, c_dp = %d, off = %d\n", column, v, r_dp, c_dp, off); */
    return min(r_dp, c_dp);
}

//// Check if DP matrix cell is empty or actually stored in the data structure
int32_t
dp_check_cell(unordered_map<int32_t, int32_t> &v_map, vector<vector<vector<dp_cell_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t vtx, int32_t d, int32_t k)
{
    int32_t v, r, c;
    if (v_map.count(vtx)) //// vertex: ok
    {
        v = v_map[vtx];
        if (diag_row_map[v].count(d)) //// diagonal ($dpd's row): ok
        {
            r = get_row(diag_row_map, v, d);
            c = get_col(d + k, k);

            if (0 <= c && c < (int32_t)dpd[v][r].size()) //// column: ok
            {
                return 1;
            }
        }
    }
    return 0;
}

//// DP matrix extension (within same vertex)
void dp_extend(vector<vector<vector<dp_cell_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v_dp, int32_t v, int32_t d, int32_t prev_k, int32_t k, FILE *out_debug)
{
    int32_t r, c, r_dp, c_dp;

    for (c_dp = prev_k; c_dp <= k; ++c_dp)
    {
        r_dp = d + c_dp;            //// row in the traditional dpd matrix
        if (r_dp >= 0 && c_dp >= 0) //// within bounds
        {
            r = get_row(diag_row_map, v, d);
            c = get_col(r_dp, c_dp);

            if (v == 0 && d == 0 && c == 0 && dpd[v][r][c].s == INT32_MAX) //// dpd[0][0][0]: first match
            {
                dpd[v][r][c].s = 0;
#ifdef DP_DEBUG
                fprintf(stdout, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, dpd[v][r][c].s);
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
                /* int32_t c_start = (int32_t)dpd[v][r].size() - 1;
                for (int32_t i = 0; i < (c - c_start); ++i) //// add as many columns up to $c, if empty they will be later filled during $d_to's extension
                    dpd[v][r].push_back({.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0}); */

                if ((prev_k < c_dp) && (dpd[v][r][c].s == INT32_MAX || (dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s)))
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
                    fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_dp, r_dp - 1, c_dp - 1, dpd[v][r][c - 1].s, v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
                }
                else if ((prev_k == c_dp) && (dpd[v][r][c].s == INT32_MAX || (dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s + 1))) //// TODO: make it >= if you want to give precedence to mismatches
                {                                                                                                                                          //// mismatch
                    if (dpd[v][r][c - 1].s < INT32_MAX && dpd[v][r][c].s > dpd[v][r][c - 1].s + 1)                                                         //// update
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
                    fprintf(stdout, "[DEBUG] Extension (X): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_dp, r_dp - 1, c_dp - 1, dpd[v][r][c - 1].s, v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
                }
            }
        }
    }
}

//// DP matrix expansion (within same vertex)
void dp_expand(vector<vector<vector<dp_cell_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t v_dp, int32_t v, int32_t v_len, int32_t d, int32_t k, char ed, FILE *out_debug)
{
    int32_t r, r_from, r_dp, c, c_from, c_dp, d_to;

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
        c = 0;
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
        fprintf(stdout, "[DEBUG] Starting mismatch (X): [%d][%d][%d] = %d\n", v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
        return;
    }

    //// diagonal (row) setup
    if (ed != 'X' && diag_row_map[v].count(d_to) == 0) //// new diagonal
    {
        dpd[v].push_back(vector<dp_cell_t>(v_len, {.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0})); //// add row to dpd[v]
        diag_row_map[v].insert({d_to, ((int32_t)dpd[v].size() - 1)});                                 //// add mapping to diag_row_map[v]
    }

    //// coordinates
    r_from = get_row(diag_row_map, v, d); //// row of the diagonal we are coming from
    c_from = get_col(d + k, k);           //// column of the diagonal we are coming from (possibly decreased by its row's offset)
    r = get_row(diag_row_map, v, d_to);
    c = get_col(r_dp, c_dp);

    //// DP and CIGAR update
    if (dpd[v][r][c].s == INT32_MAX || dpd[v][r][c].s > dpd[v][r_from][c_from].s + 1)
    {
        if (dpd[v][r][c].s > dpd[v][r_from][c_from].s + 1) //// update
        {
            free(dpd[v][r][c].op);
            free(dpd[v][r][c].bl);
        }

        dpd[v][r][c].s = dpd[v][r_from][c_from].s + 1;

        if (dpd[v][r_from][c_from].op[dpd[v][r_from][c_from].l - 1] == ed) //// if same edit
        {
            dpd[v][r][c].l = dpd[v][r_from][c_from].l - 1; //// keep same length, but decrease it for the moment to use it as index
        }
        else
        {
            dpd[v][r][c].l = dpd[v][r_from][c_from].l; //// previous length + 1 to store the additional operation
        }
        dpd[v][r][c].op = (char *)malloc((dpd[v][r][c].l + 1) * sizeof(char)); //// + 1 as l still stands for the index at the moment
        dpd[v][r][c].bl = (int32_t *)malloc((dpd[v][r][c].l + 1) * sizeof(int32_t));
        if (dpd[v][r][c].op == NULL || dpd[v][r][c].bl == NULL)
        {
            fprintf(stderr, "Allocation error at line number %d in file %s\n", __LINE__, __FILE__);
            abort();
        }

        for (int32_t i = 0; i < dpd[v][r_from][c_from].l; ++i) //// copy
        {
            dpd[v][r][c].op[i] = dpd[v][r_from][c_from].op[i];
            dpd[v][r][c].bl[i] = dpd[v][r_from][c_from].bl[i];
        }

        if (dpd[v][r][c].l == dpd[v][r_from][c_from].l - 1) //// kept same length
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
        fprintf(stdout, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", ed, v_dp, d + k, k, dpd[v][r_from][c_from].s, v_dp, r_dp, c_dp, dpd[v][r][c].s);
#endif
    }
}

//// DP matrix data structures setup when a new vertex is visited
void dp_new_vd(unordered_map<int32_t, int32_t> &v_map, vector<vector<vector<dp_cell_t>>> &dpd, vector<unordered_map<int32_t, int32_t>> &diag_row_map, int32_t w, int32_t w_len, int32_t d, int32_t ol, int32_t r_dp, int32_t c_dp, int32_t &r, int32_t &c)
{
    if (v_map.count(w) == 0) //// visiting the vertex for the first time
    {
        v_map[w] = v_map.size();
        dpd.push_back(vector<vector<dp_cell_t>>(1, vector<dp_cell_t>(w_len, {.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0})));
        diag_row_map.push_back({{d, 0}}); //// new vertex, with diagonal $d to row 0
        r = 0;
        c = ol;
    }
    else if (diag_row_map[v_map[w]].count(d) == 0) //// vertex already visited, but new diagonal
    {
        dpd[v_map[w]].push_back(vector<dp_cell_t>(w_len, {.s = INT32_MAX, .op = NULL, .bl = NULL, .l = 0})); //// add row to dpd[v]
        r = ((int32_t)dpd[v_map[w]].size() - 1);
        diag_row_map[v_map[w]].insert({d, r}); //// add mapping to diag_row_map[v]
        c = ol;
    }
    else //// vertex and diagonal already there
    {
        r = get_row(diag_row_map, v_map[w], d);
        c = get_col(r_dp, c_dp);
    }
}

//// print cigar string
void gwf_cigar(dp_cell_t cig)
{
    fprintf(stdout, "CIGAR:\t");
    for (int32_t i = 0; i < cig.l; ++i)
    {
        fprintf(stdout, "%d%c", cig.bl[i], cig.op[i]);
    }
    fprintf(stdout, "\n");
}

#endif