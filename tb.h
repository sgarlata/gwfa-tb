#ifndef TB_H
#define TB_H

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
using namespace std;

// Number of bits for each field
const int OP_BITS = 2;   // 2 bits for operations (M, X, D, I)
const int LEN_BITS = 14; // 14 bits for the length

// Maximum length for a CIGAR operation (2^14 - 1)
const uint16_t MAX_OP_LEN = (1 << LEN_BITS) - 1;

// Operations
enum class CigarOperation : uint16_t
{
    MATCH = 0,    // =
    MISMATCH = 1, // X
    DELETION = 2, // D
    INSERTION = 3 // I
};

// Function to pack a CIGAR operation and length into a single 32-bit integer
uint16_t packCigarOperation(CigarOperation operation, uint16_t length)
{
    return static_cast<uint16_t>(operation) << LEN_BITS | (length & MAX_OP_LEN);
}

// Function to unpack a packed CIGAR operation into operation and length
void unpackCigarOperation(uint16_t packedCigarOperation, CigarOperation &operation, uint16_t &length)
{
    operation = static_cast<CigarOperation>(packedCigarOperation >> LEN_BITS);
    length = packedCigarOperation & MAX_OP_LEN;
}

//// Extension for traceback (within same vertex)
void tb_extend(int32_t s, vector<uint16_t> &diag, int32_t v, int32_t d, int32_t k_old, int32_t k)
{
    int32_t c;
    CigarOperation op_old;
    uint16_t len;

    for (c = k_old + 1; c <= k; c++)
    {
        if (!diag.empty()) //// if already coming from a match
        {
            unpackCigarOperation(diag.back(), op_old, len);
            if (op_old == CigarOperation::MATCH && len < MAX_OP_LEN)
            {
                diag.pop_back();
                diag.push_back(packCigarOperation(CigarOperation::MATCH, ++len));
            }
            else //// else, add new match
                diag.push_back(packCigarOperation(CigarOperation::MATCH, 1));
        }
        else //// else, add new match
            diag.push_back(packCigarOperation(CigarOperation::MATCH, 1));

#ifdef TB_DEBUG
        int32_t r = d + c; //// row in the DP matrix
        char op_char;

        if (r == 0 && c == 0)
            fprintf(stdout, "[DEBUG] Starting match (=): [%d][%d][%d] = %d\n", v, r, c, s);
        else
            fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v, r - 1, c - 1, s, v, r, c, s);
#endif
    }
}

//// Expansion for traceback (within same vertex)
void tb_expand(int32_t s, vector<uint16_t> &diag, CigarOperation op_new, int32_t v_old, int32_t v_new, int32_t r_new, int32_t c_old, int32_t c_new)
{
    CigarOperation op_old;
    uint16_t len;

    if (!diag.empty())
    {
        unpackCigarOperation(diag.back(), op_old, len);
        if (op_old == op_new && len < MAX_OP_LEN)
        {
            diag.pop_back();
            diag.push_back(packCigarOperation(op_new, ++len));
        }
        else
            diag.push_back(packCigarOperation(op_new, 1));
    }
    else
        diag.push_back(packCigarOperation(op_new, 1));

#ifdef TB_DEBUG
    char op_char;
    int32_t r_old;

    switch (op_new)
    {
    case CigarOperation::MATCH:
        op_char = 'M';
        r_old = r_new - 1;
        break;
    case CigarOperation::MISMATCH:
        op_char = 'X';
        r_old = r_new - 1;
        break;
    case CigarOperation::DELETION:
        op_char = 'D';
        r_old = r_new;
        break;
    case CigarOperation::INSERTION:
        op_char = 'I';
        r_old = r_new - 1;
        break;
    }

    if (op_new == CigarOperation::MATCH)
        fprintf(stdout, "[DEBUG] Extension (=): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", v_old, r_old, c_old, s, v_new, r_new, c_new, s);
    else if (op_new == CigarOperation::MISMATCH && r_old == -1 && c_old == -1)
        fprintf(stdout, "[DEBUG] Starting mismatch (X): [%d][%d][%d] = %d\n", v_old, r_new, c_new, s + 1);
    else
        fprintf(stdout, "[DEBUG] Expansion (%c): [%d][%d][%d] = %d -> [%d][%d][%d] = %d\n", op_char, v_old, r_old, c_old, s, v_new, r_new, c_new, s + 1);
#endif
}

//// print cigar string
void tb_cigar(vector<uint16_t> packedCigar)
{
    CigarOperation op;
    CigarOperation op_next;
    uint16_t len;
    uint16_t len_next;
    uint32_t len_tot = 0;
    char op_char;
    FILE *fCigOut = fopen("cigar/cigar.txt", "w");
    if (fCigOut == NULL)
    {
        fprintf(stderr, "Error opening output CIGAR file.\n");
        abort();
    }

    for (vector<uint16_t>::size_type i = 0; i < packedCigar.size(); i++)
    {
        unpackCigarOperation(packedCigar[i], op, len);
        len_tot += len;

        switch (op)
        {
        case CigarOperation::MATCH:
            op_char = '=';
            break;
        case CigarOperation::MISMATCH:
            op_char = 'X';
            break;
        case CigarOperation::DELETION:
            op_char = 'D';
            break;
        case CigarOperation::INSERTION:
            op_char = 'I';
            break;
        }

        if (i + 1 < packedCigar.size())
        {
            unpackCigarOperation(packedCigar[i + 1], op_next, len_next);
            if (op != op_next)
            {
                fprintf(fCigOut, "%d%c", len_tot, op_char);
                len_tot = 0;
            }
        }
        else
            fprintf(fCigOut, "%d%c", len_tot, op_char);
    }
    fprintf(fCigOut, "\n");

    fclose(fCigOut);
}

#endif