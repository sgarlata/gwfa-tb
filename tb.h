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
    int32_t k;
    vector<uint32_t> packedCigar;
} tb_diag_t;

//// Extension for traceback (within same vertex)
void tb_extend(int32_t s, tb_diag_t &diag, int32_t v, int32_t d, int32_t k_old, int32_t k)
{
    int32_t c;
    CigarOperation op_old;
    uint32_t len;

    diag.k = k;
    for (c = k_old + 1; c <= k; c++)
    {
        if (!diag.packedCigar.empty()) //// if already coming from a match
        {
            unpackCigarOperation(diag.packedCigar.back(), op_old, len);
            if (op_old == CigarOperation::MATCH)
            {
                diag.packedCigar.pop_back();
                diag.packedCigar.push_back(packCigarOperation(CigarOperation::MATCH, ++len));
            }
            else //// else, add new match
                diag.packedCigar.push_back(packCigarOperation(CigarOperation::MATCH, 1));
        }
        else //// else, add new match
            diag.packedCigar.push_back(packCigarOperation(CigarOperation::MATCH, 1));

#ifdef TB_DEBUG
        int32_t r = d + c; //// row in the DP matrix
        if (r == 0 && c == 0)
            fprintf(stdout, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v, r, c, s);
        else
            fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v, r - 1, c - 1, s, v, r, c, s);
#endif
    }
}

//// Expansion for traceback (within same vertex)
void tb_expand(int32_t s, tb_diag_t &diag, CigarOperation op_new, int32_t v_old, int32_t v_new, int32_t r_new, int32_t c_old, int32_t c_new)
{
    CigarOperation op_old;
    uint32_t len;
    char opChar;
    int32_t r_old;

    diag.k = c_new;
    if (!diag.packedCigar.empty())
    {
        unpackCigarOperation(diag.packedCigar.back(), op_old, len);
        if (op_old == op_new)
        {
            diag.packedCigar.pop_back();
            diag.packedCigar.push_back(packCigarOperation(op_new, ++len));
        }
        else
            diag.packedCigar.push_back(packCigarOperation(op_new, 1));
    }
    else
        diag.packedCigar.push_back(packCigarOperation(op_new, 1));

    switch (op_new)
    {
    case CigarOperation::MATCH:
        opChar = 'M';
        r_old = r_new - 1;
        break;
    case CigarOperation::MISMATCH:
        opChar = 'X';
        r_old = r_new - 1;
        break;
    case CigarOperation::DELETION:
        opChar = 'D';
        r_old = r_new;
        break;
    case CigarOperation::INSERTION:
        opChar = 'I';
        r_old = r_new - 1;
        break;
    }

#ifdef TB_DEBUG
    if (op_new == CigarOperation::MATCH)
        fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_old, r_old, c_old, s, v_new, r_new, c_new, s);
    else if (op_new == CigarOperation::MISMATCH && r_old == -1 && c_old == -1)
        fprintf(stdout, "[DEBUG] Starting mismatch (X): [%d][%d][%d] = %d\n", v_old, r_new, c_new, s + 1);
    else
        fprintf(stdout, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", opChar, v_old, r_old, c_old, s, v_new, r_new, c_new, s + 1);
#endif
}

//// print cigar string
void tb_cigar(vector<uint32_t> packedCigar)
{
    CigarOperation op;
    uint32_t len;
    char opChar;
    FILE *fCigOut = fopen("cigar/cigar.txt", "w");
    if (fCigOut == NULL)
    {
        fprintf(stderr, "Error opening output CIGAR file.\n");
        abort();
    }

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
        fprintf(fCigOut, "%d%c", len, opChar);
    }
    fprintf(fCigOut, "\n");

    fclose(fCigOut);
}

#endif