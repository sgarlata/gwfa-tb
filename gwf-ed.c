#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "gwfa.h"
#include "kalloc.h"
#include "ksort.h"

/**********************
 * Indexing the graph *
 **********************/

#define arc_key(x) ((x).a)
KRADIX_SORT_INIT(gwf_arc, gwf_arc_t, arc_key, 8)

// index the graph such that we can quickly access the neighbors of a vertex
void gwf_ed_index_arc_core(uint64_t *idx, uint32_t n_vtx, uint32_t n_arc, gwf_arc_t *arc)
{
	uint32_t i, st;
	radix_sort_gwf_arc(arc, arc + n_arc);
	for (st = 0, i = 1; i <= n_arc; ++i)
	{
		if (i == n_arc || arc[i].a >> 32 != arc[st].a >> 32)
		{
			uint32_t v = arc[st].a >> 32;
			assert(v < n_vtx);
			idx[v] = (uint64_t)st << 32 | (i - st);
			st = i;
		}
	}
}

void gwf_ed_index(void *km, gwf_graph_t *g)
{
	KMALLOC(km, g->aux, g->n_vtx);
	gwf_ed_index_arc_core(g->aux, g->n_vtx, g->n_arc, g->arc);
}

// free the index
void gwf_cleanup(void *km, gwf_graph_t *g)
{
	kfree(km, g->aux);
	g->aux = 0;
}

/**************************************
 * Graph WaveFront with edit distance *
 **************************************/

#include "khashl.h" // make it compatible with kalloc
#include "kdq.h"
#include "kvec.h"

#define GWF_DIAG_SHIFT 0x40000000 //// Hex2Bin: 1000000000000000000000000000000 (31 digits)

static inline uint64_t gwf_gen_vd(uint32_t v, int32_t d) //// combine v and d
{
	return (uint64_t)v << 32 | (GWF_DIAG_SHIFT + d);
}

/*
 * Diagonal interval
 */
typedef struct
{
	uint64_t vd0, vd1;
} gwf_intv_t;

typedef kvec_t(gwf_intv_t) gwf_intv_v;

#define intvd_key(x) ((x).vd0)
KRADIX_SORT_INIT(gwf_intv, gwf_intv_t, intvd_key, 8)

static int gwf_intv_is_sorted(int32_t n_a, const gwf_intv_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i - 1].vd0 > a[i].vd0)
			break;
	return (i == n_a);
}

void gwf_ed_print_intv(size_t n, gwf_intv_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i)
		printf("Z\t%d\t%d\t%d\n", (int32_t)(a[i].vd0 >> 32), (int32_t)a[i].vd0 - GWF_DIAG_SHIFT, (int32_t)a[i].vd1 - GWF_DIAG_SHIFT);
}

// merge overlapping intervals; input must be sorted
static size_t gwf_intv_merge_adj(size_t n, gwf_intv_t *a)
{
	size_t i, k;
	uint64_t st, en;
	if (n == 0)
		return 0;
	st = a[0].vd0, en = a[0].vd1;
	for (i = 1, k = 0; i < n; ++i)
	{
		if (a[i].vd0 > en)
		{
			a[k].vd0 = st, a[k++].vd1 = en;
			st = a[i].vd0, en = a[i].vd1;
		}
		else
			en = en > a[i].vd1 ? en : a[i].vd1;
	}
	a[k].vd0 = st, a[k++].vd1 = en;
	return k;
}

// merge two sorted interval lists
static size_t gwf_intv_merge2(gwf_intv_t *a, size_t n_b, const gwf_intv_t *b, size_t n_c, const gwf_intv_t *c)
{
	size_t i = 0, j = 0, k = 0;
	while (i < n_b && j < n_c)
	{
		if (b[i].vd0 <= c[j].vd0)
			a[k++] = b[i++];
		else
			a[k++] = c[j++];
	}
	while (i < n_b)
		a[k++] = b[i++];
	while (j < n_c)
		a[k++] = c[j++];
	return gwf_intv_merge_adj(k, a);
}

/*
 * Diagonal
 */
typedef struct
{				 // a diagonal
	uint64_t vd; // higher 32 bits: vertex ID; lower 32 bits: diagonal+0x4000000 //// Why not split them?
	int32_t k;	 //// wavefront position on the vertex
	uint32_t xo; // higher 31 bits: anti diagonal; lower 1 bit: out-of-order or not
	int32_t t;
} gwf_diag_t;

typedef kvec_t(gwf_diag_t) gwf_diag_v;

#define ed_key(x) ((x).vd)
KRADIX_SORT_INIT(gwf_ed, gwf_diag_t, ed_key, 8)

KDQ_INIT(gwf_diag_t)

void gwf_ed_print_diag(size_t n, gwf_diag_t *a) // for debugging only
{
	size_t i;
	for (i = 0; i < n; ++i)
	{
		int32_t d = (int32_t)a[i].vd - GWF_DIAG_SHIFT;
		printf("Z\t%d\t%d\t%d\t%d\t%d\n", (int32_t)(a[i].vd >> 32), d, a[i].k, d + a[i].k, a[i].xo >> 1);
	}
}

// push (v,d,k) to the end of the queue
static inline void gwf_diag_push(void *km, gwf_diag_v *a, uint32_t v, int32_t d, int32_t k, uint32_t x, uint32_t ooo, int32_t t)
{
	gwf_diag_t *p;
	kv_pushp(gwf_diag_t, km, *a, &p); //// push a pointer to the newly generated diagonal p
	p->vd = gwf_gen_vd(v, d), p->k = k, p->xo = x << 1 | ooo, p->t = t;
}

// determine the wavefront on diagonal (v,d)
static inline int32_t gwf_diag_update(gwf_diag_t *p, uint32_t v, int32_t d, int32_t k, uint32_t x, uint32_t ooo, int32_t t)
{
	uint64_t vd = gwf_gen_vd(v, d);
	if (p->vd == vd)
	{
		p->xo = p->k > k ? p->xo : x << 1 | ooo;
		p->t = p->k > k ? p->t : t;
		p->k = p->k > k ? p->k : k;
		return 0;
	}
	return 1;
}

static int gwf_diag_is_sorted(int32_t n_a, const gwf_diag_t *a)
{
	int32_t i;
	for (i = 1; i < n_a; ++i)
		if (a[i - 1].vd > a[i].vd)
			break;
	return (i == n_a);
}

// sort a[]. This uses the gwf_diag_t::ooo field to speed up sorting.
static void gwf_diag_sort(int32_t n_a, gwf_diag_t *a, void *km, gwf_diag_v *ooo)
{
	int32_t i, j, k, n_b, n_c;
	gwf_diag_t *b, *c;

	kv_resize(gwf_diag_t, km, *ooo, n_a);
	for (i = 0, n_c = 0; i < n_a; ++i)
		if (a[i].xo & 1)
			++n_c;
	n_b = n_a - n_c;
	b = ooo->a, c = b + n_b;
	for (i = j = k = 0; i < n_a; ++i)
	{
		if (a[i].xo & 1)
			c[k++] = a[i];
		else
			b[j++] = a[i];
	}
	radix_sort_gwf_ed(c, c + n_c);
	for (k = 0; k < n_c; ++k)
		c[k].xo &= 0xfffffffeU;

	i = j = k = 0;
	while (i < n_b && j < n_c)
	{
		if (b[i].vd <= c[j].vd)
			a[k++] = b[i++];
		else
			a[k++] = c[j++];
	}
	while (i < n_b)
		a[k++] = b[i++];
	while (j < n_c)
		a[k++] = c[j++];
}

// remove diagonals not on the wavefront
static int32_t gwf_diag_dedup(int32_t n_a, gwf_diag_t *a, void *km, gwf_diag_v *ooo)
{
	int32_t i, n, st;
	if (!gwf_diag_is_sorted(n_a, a))
		gwf_diag_sort(n_a, a, km, ooo);
	for (i = 1, st = 0, n = 0; i <= n_a; ++i)
	{
		if (i == n_a || a[i].vd != a[st].vd)
		{
			int32_t j, max_j = st;
			if (st + 1 < i)
				for (j = st + 1; j < i; ++j) // choose the far end (i.e. the wavefront)
					if (a[max_j].k < a[j].k)
						max_j = j;
			a[n++] = a[max_j];
			st = i;
		}
	}
	return n;
}

// use forbidden bands to remove diagonals not on the wavefront
static int32_t gwf_mixed_dedup(int32_t n_a, gwf_diag_t *a, int32_t n_b, gwf_intv_t *b)
{
	int32_t i = 0, j = 0, k = 0;
	while (i < n_a && j < n_b)
	{
		if (a[i].vd >= b[j].vd0 && a[i].vd < b[j].vd1)
			++i;
		else if (a[i].vd >= b[j].vd1)
			++j;
		else
			a[k++] = a[i++];
	}
	while (i < n_a)
		a[k++] = a[i++];
	return k;
}

/*
 * Traceback stack
 */
KHASHL_MAP_INIT(KH_LOCAL, gwf_map64_t, gwf_map64, uint64_t, int32_t, kh_hash_uint64, kh_eq_generic)

typedef struct
{
	int32_t v;
	int32_t pre;
} gwf_trace_t;

typedef kvec_t(gwf_trace_t) gwf_trace_v;

static int32_t gwf_trace_push(void *km, gwf_trace_v *a, int32_t v, int32_t pre, gwf_map64_t *h)
{
	uint64_t key = (uint64_t)v << 32 | (uint32_t)pre;
	khint_t k;
	int absent;
	k = gwf_map64_put(h, key, &absent);
	if (absent)
	{
		gwf_trace_t *p;
		kv_pushp(gwf_trace_t, km, *a, &p);
		p->v = v, p->pre = pre;
		kh_val(h, k) = a->n - 1;
		return a->n - 1;
	}
	return kh_val(h, k);
}

/*
 * Core GWFA routine
 */
KHASHL_INIT(KH_LOCAL, gwf_set64_t, gwf_set64, uint64_t, kh_hash_dummy, kh_eq_generic)

typedef struct
{
	void *km;			  //// chunk of memory, see "kalloc.c"
	gwf_set64_t *ha;	  // hash table for adjacency
	gwf_map64_t *ht;	  // hash table for traceback
	gwf_intv_v intv;	  //// dynamic array of diagonal intervals
	gwf_intv_v tmp, swap; //// same as above, used as aux support
	gwf_diag_v ooo;		  //// "out-of-order", dynamic array of diagonals, used to speed up sorting
	gwf_trace_v t;		  //// dynamic array of traces
} gwf_edbuf_t;

// remove diagonals not on the wavefront
static int32_t gwf_dedup(gwf_edbuf_t *buf, int32_t n_a, gwf_diag_t *a)
{
	if (buf->intv.n + buf->tmp.n > 0)
	{
		if (!gwf_intv_is_sorted(buf->tmp.n, buf->tmp.a))
			radix_sort_gwf_intv(buf->tmp.a, buf->tmp.a + buf->tmp.n);
		kv_copy(gwf_intv_t, buf->km, buf->swap, buf->intv);
		kv_resize(gwf_intv_t, buf->km, buf->intv, buf->intv.n + buf->tmp.n);
		buf->intv.n = gwf_intv_merge2(buf->intv.a, buf->swap.n, buf->swap.a, buf->tmp.n, buf->tmp.a);
	}
	n_a = gwf_diag_dedup(n_a, a, buf->km, &buf->ooo);
	if (buf->intv.n > 0)
		n_a = gwf_mixed_dedup(n_a, a, buf->intv.n, buf->intv.a);
	return n_a;
}

// remove diagonals that lag far behind the furthest wavefront
static int32_t gwf_prune(int32_t n_a, gwf_diag_t *a, uint32_t max_lag)
{
	int32_t i, j;
	uint32_t max_x = 0;
	for (i = 0; i < n_a; ++i)
		max_x = max_x > a[i].xo >> 1 ? max_x : a[i].xo >> 1;
	if (max_x <= max_lag)
		return n_a; // no filtering
	for (i = j = 0; i < n_a; ++i)
		if ((a[i].xo >> 1) + max_lag >= max_x)
			a[j++] = a[i];
	return j;
}

// reach the wavefront
static inline int32_t gwf_extend1(int32_t d, int32_t k, int32_t vl, const char *ts, int32_t ql, const char *qs)
{
	int32_t max_k = (ql - d < vl ? ql - d : vl) - 1; //// max wavefront position = min(query length - diagonal, label length) - 1
	const char *ts_ = ts + 1, *qs_ = qs + d + 1;
#if 1 //// unoptimized, but easier to understand
	// int32_t i = k + d; while (k + 1 < g->len[v] && i + 1 < ql && g->seq[v][k+1] == q[i+1]) ++k, ++i;
	while (k < max_k && *(ts_ + k) == *(qs_ + k)) //// LCP: extending along exact matches to find the furthest cell
		++k;
#else
	uint64_t cmp = 0;
	while (k + 7 < max_k) //// why 7? perhaps because of byte dimension?
	{
		uint64_t x = *(uint64_t *)(ts_ + k); // warning: unaligned memory access
		uint64_t y = *(uint64_t *)(qs_ + k);
		cmp = x ^ y;  //// bitwise-exclusive-OR (bit set to 0 if corresponding operands' bits are equal)
		if (cmp == 0) //// if x and y are bitwise equal
			k += 8;	  //// 1 byte
		else
			break;
	}
	if (cmp)
		k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr //// Bit Scan Reverse
	else if (k + 7 >= max_k)
		while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
			++k;
#endif
	return k;
}

// This is essentially Landau-Vishkin for linear sequences. The function speeds up alignment to long vertices. Not really necessary.
static void gwf_ed_extend_batch(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t n, gwf_diag_t *a, gwf_diag_v *B,
								kdq_t(gwf_diag_t) * A, gwf_intv_v *tmp_intv)
{
	int32_t j, m;
	int32_t v = a->vd >> 32;
	int32_t vl = g->len[v];
	const char *ts = g->seq[v];
	gwf_diag_t *b;

	// wfa_extend
	for (j = 0; j < n; ++j)
	{
		int32_t k;
		k = gwf_extend1((int32_t)a[j].vd - GWF_DIAG_SHIFT, a[j].k, vl, ts, ql, q);
		a[j].xo += (k - a[j].k) << 2;
		a[j].k = k;
	}

	// wfa_next
	kv_resize(gwf_diag_t, km, *B, B->n + n + 2);
	b = &B->a[B->n];
	b[0].vd = a[0].vd - 1;
	b[0].xo = a[0].xo + 2; // 2 == 1<<1
	b[0].k = a[0].k + 1;
	b[0].t = a[0].t;
	b[1].vd = a[0].vd;
	b[1].xo = n == 1 || a[0].k > a[1].k ? a[0].xo + 4 : a[1].xo + 2;
	b[1].t = n == 1 || a[0].k > a[1].k ? a[0].t : a[1].t;
	b[1].k = (n == 1 || a[0].k > a[1].k ? a[0].k : a[1].k) + 1;
	for (j = 1; j < n - 1; ++j)
	{
		uint32_t x = a[j - 1].xo + 2;
		int32_t k = a[j - 1].k, t = a[j - 1].t;
		x = k > a[j].k + 1 ? x : a[j].xo + 4;
		t = k > a[j].k + 1 ? t : a[j].t;
		k = k > a[j].k + 1 ? k : a[j].k + 1;
		x = k > a[j + 1].k + 1 ? x : a[j + 1].xo + 2;
		t = k > a[j + 1].k + 1 ? t : a[j + 1].t;
		k = k > a[j + 1].k + 1 ? k : a[j + 1].k + 1;
		b[j + 1].vd = a[j].vd, b[j + 1].k = k, b[j + 1].xo = x, b[j + 1].t = t;
	}
	if (n >= 2)
	{
		b[n].vd = a[n - 1].vd;
		b[n].xo = a[n - 2].k > a[n - 1].k + 1 ? a[n - 2].xo + 2 : a[n - 1].xo + 4;
		b[n].t = a[n - 2].k > a[n - 1].k + 1 ? a[n - 2].t : a[n - 1].t;
		b[n].k = a[n - 2].k > a[n - 1].k + 1 ? a[n - 2].k : a[n - 1].k + 1;
	}
	b[n + 1].vd = a[n - 1].vd + 1;
	b[n + 1].xo = a[n - 1].xo + 2;
	b[n + 1].t = a[n - 1].t;
	b[n + 1].k = a[n - 1].k;

	// drop out-of-bound cells
	for (j = 0; j < n; ++j)
	{
		gwf_diag_t *p = &a[j];
		if (p->k == vl - 1 || (int32_t)p->vd - GWF_DIAG_SHIFT + p->k == ql - 1)
			p->xo |= 1, *kdq_pushp(gwf_diag_t, A) = *p;
	}
	for (j = 0, m = 0; j < n + 2; ++j)
	{
		gwf_diag_t *p = &b[j];
		int32_t d = (int32_t)p->vd - GWF_DIAG_SHIFT;
		if (d + p->k < ql && p->k < vl)
		{
			b[m++] = *p;
		}
		else if (p->k == vl)
		{
			gwf_intv_t *q;
			kv_pushp(gwf_intv_t, km, *tmp_intv, &q);
			q->vd0 = gwf_gen_vd(v, d), q->vd1 = q->vd0 + 1;
		}
	}
	B->n += m;
}

// wfa_extend and wfa_next combined
static gwf_diag_t *gwf_ed_extend(gwf_edbuf_t *buf, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v1, uint32_t max_lag, int32_t traceback,
								 int32_t *end_v, int32_t *end_off, int32_t *end_tb, int32_t *n_a_, gwf_diag_t *a, int32_t s, int32_t ***dp)
{
	int32_t i, x, n = *n_a_, do_dedup = 1; //// do_dedup is a binary flag used to know when to remove diagonals not on the wavefront
	kdq_t(gwf_diag_t) * A;				   //// queue to keep track of the diagonals on which the wavefront can be further updated
	gwf_diag_v B = {0, 0, 0};			   //// dynamic array of diagonals: gwf_diag_v is a typedef of kvec_t(gwf_diag_t) (paper's Q)
	gwf_diag_t *b;						   //// array of diagonals to which B will be copied at the end

	*end_v = *end_off = *end_tb = -1;
	buf->tmp.n = 0;
	gwf_set64_clear(buf->ha);				 // hash table $h to avoid visiting a vertex twice
	for (i = 0, x = 1; i < 32; ++i, x <<= 1) //// x left-wise-shifted by 1 bit
		if (x >= n)							 //// $x is probably used later for the batch method used to speed up alignment -> not relevant for us
			break;
	if (i < 4)
		i = 4;
	A = kdq_init2(gwf_diag_t, buf->km, i);	  // $A is a queue
	kv_resize(gwf_diag_t, buf->km, B, n * 2); //// to properly resize the queue (if not large enough)
#if 1										  // unoptimized version without calling gwf_ed_extend_batch() at all. The final result will be the same.
	A->count = n;
	memcpy(A->a, a, n * sizeof(*a)); //// $a is copied to $A
#else								 // optimized for long vertices.
	for (x = 0, i = 1; i <= n; ++i)
	{
		if (i == n || a[i].vd != a[i - 1].vd + 1)
		{
			gwf_ed_extend_batch(buf->km, g, ql, q, i - x, &a[x], &B, A, &buf->tmp); //// To speed up alignment
			x = i;
		}
	}
	if (kdq_size(A) == 0) //// if the queue is empty (no more diagonals to further update the wavefront are available)
		do_dedup = 0;	  //// flag set to false as no more diagonals are there to be removed
#endif
	kfree(buf->km, a); // $a is not used as it has been copied to $A

	while (kdq_size(A)) //// while there are still diagonals to update the wavefront
	{
		gwf_diag_t t;				 //// single diagonal
		uint32_t x0;				 //// anti diagonal
		int32_t ooo, v, d, k, i, vl; //// $ooo: "out-of-order", $v: vertex ID, $d: diagonal (paper's k)

		t = *kdq_shift(gwf_diag_t, A);		//// store in $t the vertex+diagonal on the queue head
		ooo = t.xo & 1, v = t.vd >> 32;		// vertex //// bitwise AND with 1 to keep just the lower 1 bit (flag for out-of-order); right shift to keep just the higher 32 bits (vertex ID)
		d = (int32_t)t.vd - GWF_DIAG_SHIFT; // diagonal //// vd (uint64_t) structure: higher 32 bits: vertex ID; lower 32 bits: diagonal + GWF_DIAG_SHIFT
		k = t.k;							// wavefront position on the vertex //// for the given diagonal (paper's j)
		vl = g->len[v];						// $vl is the vertex length

		k = gwf_extend1(d, k, vl, g->seq[v], ql, q);
		i = k + d; // query position //// (paper's "i = Hvk + k")

		for (int l = 0; l <= k; ++l)
			if (d + l >= 0 && (dp[v][d + l][l] == -1 || dp[v][d + l][l] > s))
				dp[v][d + l][l] = s;

		x0 = (t.xo >> 1) + ((k - t.k) << 1); // current anti diagonal
		//// EXPANSION
		if (k + 1 < vl && i + 1 < ql)
		{ // the most common case: the wavefront is in the middle
			int32_t push1 = 1, push2 = 1;
			if (B.n >= 2)
				push1 = gwf_diag_update(&B.a[B.n - 2], v, d - 1, k + 1, x0 + 1, ooo, t.t);
			if (B.n >= 1)
				push2 = gwf_diag_update(&B.a[B.n - 1], v, d, k + 1, x0 + 2, ooo, t.t);
			if (push1)
			{
				gwf_diag_push(buf->km, &B, v, d - 1, k + 1, x0 + 1, 1, t.t); //// diagonal above
				if (i >= 0 && (dp[v][(k + 1) + (d - 1)][k + 1] == -1 || dp[v][(k + 1) + (d - 1)][k + 1] > s + 1))
					dp[v][(k + 1) + (d - 1)][k + 1] = s + 1;
			}
			if (push2 || push1)
			{
				gwf_diag_push(buf->km, &B, v, d, k + 1, x0 + 2, 1, t.t); //// current diagonal
				if (dp[v][(k + 1) + d][k + 1] == -1 || dp[v][(k + 1) + d][k + 1] > s + 1)
					dp[v][(k + 1) + d][k + 1] = s + 1;
			}
			gwf_diag_push(buf->km, &B, v, d + 1, k, x0 + 1, ooo, t.t); //// diagonal below
			if (i >= 0 && (dp[v][k + (d + 1)][k] == -1 || dp[v][k + (d + 1)][k] > s + 1))
				dp[v][k + (d + 1)][k] = s + 1;
		}
		else if (i + 1 < ql)
		{																				  // k + 1 == g->len[v]; reaching the end of the vertex but not the end of query
			int32_t ov = g->aux[v] >> 32, nv = (int32_t)g->aux[v], j, n_ext = 0, tw = -1; //// $nv: number of v's neighbors
			gwf_intv_t *p;
			kv_pushp(gwf_intv_t, buf->km, buf->tmp, &p);
			p->vd0 = gwf_gen_vd(v, d), p->vd1 = p->vd0 + 1;
			if (traceback)
				tw = gwf_trace_push(buf->km, &buf->t, v, t.t, buf->ht);
			for (j = 0; j < nv; ++j)
			{											 // traverse $v's neighbors
				uint32_t w = (uint32_t)g->arc[ov + j].a; // $w is next to $v
				int32_t ol = g->arc[ov + j].o;
				int absent;
				gwf_set64_put(buf->ha, (uint64_t)w << 32 | (i + 1), &absent); // test if ($w,$i) has been visited
				if (q[i + 1] == g->seq[w][ol])
				{			 // can be extended to the next vertex without a mismatch
					++n_ext; //// EXTENSION
					if (absent)
					{
						gwf_diag_t *p;
						p = kdq_pushp(gwf_diag_t, A);
						p->vd = gwf_gen_vd(w, i + 1 - ol), p->k = ol, p->xo = (x0 + 2) << 1 | 1, p->t = tw;
						if (dp[w][i + 1][ol] == -1 || dp[w][i + 1][ol] > s)
							dp[w][i + 1][ol] = s;
					}
				}
				else if (absent) //// EXPANSION
				{
					gwf_diag_push(buf->km, &B, w, i - ol, ol, x0 + 1, 1, tw); //// w's diagonal above wrt v's
					if (dp[w][i][ol] == -1 || dp[w][i][ol] > s + 1)
						dp[w][i][ol] = s + 1;

					gwf_diag_push(buf->km, &B, w, i + 1 - ol, ol, x0 + 2, 1, tw); //// w's current diagonal wrt v's
					if (dp[w][i + 1][ol] == -1 || dp[w][i + 1][ol] > s + 1)
						dp[w][i + 1][ol] = s + 1;
				}
			}
			if (nv == 0 || n_ext != nv) // add an insertion to the target; this *might* cause a duplicate in corner cases
				gwf_diag_push(buf->km, &B, v, d + 1, k, x0 + 1, 1, t.t);
			if (dp[v][d + 1 + k][k] == -1 || dp[v][d + 1 + k][k] > s + 1) //// probably redundant
				dp[v][d + 1 + k][k] = s + 1;
		}
		else if (v1 < 0 || (v == v1 && k + 1 == vl))
		{ // i + 1 == ql
			*end_v = v, *end_off = k, *end_tb = t.t, *n_a_ = 0;
			kdq_destroy(gwf_diag_t, A);
			kfree(buf->km, B.a);
			return 0;
		}
		else if (k + 1 < vl)
		{																   // i + 1 == ql; reaching the end of the query but not the end of the vertex
			gwf_diag_push(buf->km, &B, v, d - 1, k + 1, x0 + 1, ooo, t.t); // add an deletion; this *might* case a duplicate in corner cases
			if (dp[v][d - k][k + 1] == -1 || dp[v][d + 1 + k][k + 1] > s + 1)
				dp[v][d - k][k + 1] = s + 1;
		}
		else if (v != v1)
		{ // i + 1 == ql && k + 1 == g->len[v]; not reaching the last vertex $v1
			int32_t ov = g->aux[v] >> 32, nv = (int32_t)g->aux[v], j, tw = -1;
			if (traceback)
				tw = gwf_trace_push(buf->km, &buf->t, v, t.t, buf->ht);
			for (j = 0; j < nv; ++j)
			{
				uint32_t w = (uint32_t)g->arc[ov + j].a;
				int32_t ol = g->arc[ov + j].o;
				gwf_diag_push(buf->km, &B, w, i - ol, ol, x0 + 1, 1, tw); // deleting the first base on the next vertex
				if (dp[w][i][ol] == -1 || dp[w][i][ol] > s + 1)
					dp[w][i][ol] = s + 1;
			}
		}
		else
			assert(0); // should never come here
	}

	kdq_destroy(gwf_diag_t, A);
	*n_a_ = n = B.n, b = B.a; //// b <- B

	if (do_dedup)
		*n_a_ = n = gwf_dedup(buf, n, b);
	if (max_lag > 0)
		*n_a_ = n = gwf_prune(n, b, max_lag);
	return b;
}

static void gwf_traceback(gwf_edbuf_t *buf, int32_t end_v, int32_t end_tb, gwf_path_t *path)
{
	int32_t i = end_tb, n = 1;
	while (i >= 0 && buf->t.a[i].v >= 0)
		++n, i = buf->t.a[i].pre;
	KMALLOC(buf->km, path->v, n);
	i = end_tb, n = 0;
	path->v[n++] = end_v;
	while (i >= 0 && buf->t.a[i].v >= 0)
		path->v[n++] = buf->t.a[i].v, i = buf->t.a[i].pre;
	path->nv = n;
	for (i = 0; i < path->nv >> 1; ++i)
		n = path->v[i], path->v[i] = path->v[path->nv - 1 - i], path->v[path->nv - 1 - i] = n;
}
//// ALGORITHM CORE WRAPPER
int32_t gwf_ed(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v0, int32_t v1, uint32_t max_lag, int32_t traceback, gwf_path_t *path)
{
	int32_t s = 0, n_a = 1, end_tb; //// $s: edit distance, $n_a: number of diagonals on which WF can be updated, $end_tb: end traceback
	gwf_diag_t *a;					//// array of diagonals
	gwf_edbuf_t buf;				//// ??? perhaps a struct (buffer) to store temporary alignment information per single read
	int32_t ***dp;
	dp = (int32_t ***)malloc(g->n_vtx * sizeof(int32_t **)); //// dynamic programming matrix

	for (int v = 0; v < g->n_vtx; ++v)
	{
		dp[v] = (int **)malloc(ql * sizeof(int32_t *));

		for (int i = 0; i < ql; ++i)
		{
			dp[v][i] = (int *)malloc(g->len[v] * sizeof(int32_t));

			for (int j = 0; j < g->len[v]; ++j)
			{
				dp[v][i][j] = -1;
			}
		}
	}

	FILE *out = fopen("dp.csv", "w+");

	if (out == NULL)
	{
		fprintf(stderr, "Error opening DP file\n");
		return 1;
	}

	memset(&buf, 0, sizeof(buf)); //// buffer initialization
	buf.km = km;				  //// memory chunk, see "kalloc.c"
	buf.ha = gwf_set64_init2(km); //// initialization of hash table for adjacency
	buf.ht = gwf_map64_init2(km); //// initialization of hash table for traceback
	kv_resize(gwf_trace_t, km, buf.t, g->n_vtx + 16);
	KCALLOC(km, a, 1);
	a[0].vd = gwf_gen_vd(v0, 0), a[0].k = -1, a[0].xo = 0; // the initial state
	if (traceback)
		a[0].t = gwf_trace_push(km, &buf.t, -1, -1, buf.ht); //// traceback info for the initial state
	while (n_a > 0)
	{
		a = gwf_ed_extend(&buf, g, ql, q, v1, max_lag, traceback, &path->end_v, &path->end_off, &end_tb, &n_a, a, s, dp);
		if (path->end_off >= 0 || n_a == 0)
			break;
		++s; //// increase edit distance (alignment cost)
#ifdef GWF_DEBUG
		printf("[%s] dist=%d, n=%d, n_intv=%ld, n_tb=%ld\n", __func__, s, n_a, buf.intv.n, buf.t.n);
#endif
	}
	if (traceback)
		gwf_traceback(&buf, path->end_v, end_tb, path);
	gwf_set64_destroy(buf.ha);
	gwf_map64_destroy(buf.ht);
	kfree(km, buf.intv.a);
	kfree(km, buf.tmp.a);
	kfree(km, buf.swap.a);
	kfree(km, buf.t.a);

	fprintf(out, ",");

	for (int v = 0; v < g->n_vtx; ++v)
	{
		for (int l = 0; l < g->len[v]; ++l)
		{
			fprintf(out, "%c,", g->seq[v][l]);
		}
		fprintf(out, ",");
	}
	fprintf(out, "\n");

	for (int i = 0; i < ql; ++i)
	{
		fprintf(out, "%c,", q[i]);
		for (int v = 0; v < g->n_vtx; ++v)
		{
			for (int j = 0; j < g->len[v]; ++j)
			{
				if (0 <= dp[v][i][j] && dp[v][i][j] <= s)
					fprintf(out, "%d,", dp[v][i][j]);
				else
					fprintf(out, ",");
			}
			fprintf(out, ",");
		}
		if (i < ql - 1)
			fprintf(out, "\n");
	}

	//// FREE THE DP MATRIX
	for (int v = 0; v < g->n_vtx; ++v)
	{
		for (int i = 0; i < ql; ++i)
		{
			free(dp[v][i]);
		}
		free(dp[v]);
	}
	free(dp);
	fclose(out);

	path->s = path->end_v >= 0 ? s : -1;
	return path->s; // end_v < 0 could happen if v0 can't reach v1
}
