#ifndef GWFA_H
#define GWFA_H

#include <stdint.h>
#include <string>
using namespace std;

typedef struct
{
	uint64_t a; // first 32 bits for starting vertex and last 32 bits for ending vertex
	int32_t o;	// overlap between the two
} gwf_arc_t;

typedef struct
{
	uint32_t n_vtx, n_arc;
	uint32_t *len;
	uint32_t *src;
	char **seq;
	gwf_arc_t *arc;
	uint64_t *aux;
} gwf_graph_t;

typedef struct
{
	int32_t s;
	int32_t end_v, end_off;
	int32_t nv;
	int32_t *v;
} gwf_path_t;

#ifdef __cplusplus
extern "C"
{
#endif

	void gwf_ed_index(void *km, gwf_graph_t *g);
	void gwf_cleanup(void *km, gwf_graph_t *g);
	void gwf_ed(void *km, const gwf_graph_t *g, int32_t ql, const char *q, int32_t v0, int32_t v1, uint32_t max_lag, gwf_path_t *path, string &cig, int32_t &abl);

#ifdef __cplusplus
}
#endif

#endif
