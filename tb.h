#ifndef TB_H
#define TB_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
using namespace std;

// Number of bits for each field
const int OP_BITS = 2;   // 2 bits for operations (M, X, D, I)
const int LEN_BITS = 30; // 30 bits for the length

// Maximum length for a CIGAR operation (2^30 - 1)
const uint32_t MAX_OP_LEN = (1 << LEN_BITS) - 1;

// Operations
enum class CigarOperation : uint32_t
{
    MATCH = 0,    // =
    MISMATCH = 1, // X
    DELETION = 2, // D
    INSERTION = 3 // I
};

// Function to pack a CIGAR operation and length into a single 32-bit integer
uint32_t packCigarOperation(CigarOperation operation, uint32_t length)
{
    if (length > MAX_OP_LEN)
    {
        fprintf(stderr, "CIGAR error: the length of the operation exceeds the maximum allowed.\n");
        return 1;
    }

    return static_cast<uint32_t>(operation) << LEN_BITS | (length & MAX_OP_LEN);
}

// Function to unpack a packed CIGAR operation into operation and length
void unpackCigarOperation(uint32_t packedCigarOperation, CigarOperation &operation, uint32_t &length)
{
    operation = static_cast<CigarOperation>(packedCigarOperation >> LEN_BITS);
    length = packedCigarOperation & MAX_OP_LEN;
}

//// DIAGONAL TYPE FOR TRACEBACK
typedef struct TB_DIAG
{
    int32_t s; //// edit distance
    vector<uint32_t> packedCigar;
} tb_diag_t;

/* //// For DEBUG purposes: print current tracebacks
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
} */

//// Extension for traceback (within same vertex)
void tb_extend(int32_t s, tb_diag_t &diag, int32_t v, int32_t d, int32_t k_old, int32_t k)
{
    int32_t r, c;
    CigarOperation op;
    uint32_t len;

    for (c = k_old; c <= k; ++c)
    {
        if (c > k_old || diag.s > s)
        {
            if (!diag.packedCigar.empty()) //// if already coming from a match
            {
                unpackCigarOperation(diag.packedCigar.back(), op, len);
                if (op == CigarOperation::MATCH)
                    len++;
                else //// else, add new match
                    diag.packedCigar.push_back(packCigarOperation(CigarOperation::MATCH, 1));
            }
            else //// else, add new match
                diag.packedCigar.push_back(packCigarOperation(CigarOperation::MATCH, 1));
#ifdef TB_DEBUG
            int32_t r = d + c; //// row in the DP matrix
            if (r == 0 && c == 0)
                fprintf(stdout, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v, r, c, diag.s);
            else
                fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v, r - 1, c - 1, diag.s, v, r, c, diag.s);
#endif

#ifdef TB_PRINT
            print_tb(wf);
#endif
        }
    }
}

//// Expansion for traceback (within same vertex)
void tb_expand(int32_t s, tb_diag_t &diag, int32_t v, int32_t d, int32_t k, CigarOperation op)
{
    int32_t r, r_from, r_dp, c_dp, d_to;

    switch (op)
    {
    case CigarOperation::DELETION:
        d_to = d - 1;
        r_dp = k + d;
        c_dp = k + 1;
        if (r_dp < 0 || c_dp <= 0)
            return;
        break;
    case CigarOperation::MISMATCH:
        d_to = d;
        r_dp = k + d + 1;
        c_dp = k + 1;
        if (r_dp < 0 || c_dp < 0)
            return;
        break;
    case CigarOperation::INSERTION:
        d_to = d + 1;
        r_dp = k + d + 1;
        c_dp = k;
        if (r_dp <= 0 || c_dp < 0)
            return;
        break;
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

//// print cigar string
void tb_cigar(vector<uint32_t> packedCigar)
{
    CigarOperation op;
    uint32_t len;
    char opChar;

    fprintf(stdout, "CIGAR:\t");
    for (const uint32_t packedCigarOperation : packedCigar)
    {
        unpackCigarOperation(packedCigarOperation, op, len);
        switch (op)
        {
        case CigarOperation::MATCH:
            opChar = '=';
            break;
        case CigarOperation::MISMATCH:
            opChar = 'X';
            break;
        case CigarOperation::DELETION:
            opChar = 'D';
            break;
        case CigarOperation::INSERTION:
            opChar = 'I';
            break;
        }
        fprintf(stdout, "%d%c", len, opChar);
    }
    fprintf(stdout, "\n");
}

#endif