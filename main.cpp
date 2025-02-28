#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include "gfa.h"
#include "gfa-priv.h"
#include "gwfa.h"
#include "ketopt.h"
#include "kalloc.h"
#include "kseq.h"
#include <vector>
#include <time.h>
#include <iostream>
KSEQ_INIT(gzFile, gzread)

void gwf_sub_free(gfa_sub_t *sub)
{
	kfree(sub->km, sub->v);
	kfree(sub->km, sub->a);
	// km_destroy(sub->km);
	kfree(0, sub);
}

gwf_graph_t *gwf_gfa2gwf(const gfa_t *gfa, uint32_t v0)
{
	int32_t i, k;
	gwf_graph_t *g;
	gfa_sub_t *sub;
	sub = gfa_sub_from(0, gfa, v0, 1 << 30);
	GFA_CALLOC(g, 1);
	g->n_vtx = sub->n_v;
	g->n_arc = sub->n_a;
	GFA_MALLOC(g->len, g->n_vtx);
	GFA_MALLOC(g->src, g->n_vtx);
	GFA_MALLOC(g->seq, g->n_vtx);
	GFA_MALLOC(g->arc, g->n_arc);
	for (i = k = 0; i < sub->n_v; ++i)
	{
		uint32_t v = sub->v[i].v, len = gfa->seg[v >> 1].len, j;
		const gfa_seg_t *s = &gfa->seg[v >> 1];
		g->len[i] = len;
		g->src[i] = v;
		GFA_MALLOC(g->seq[i], len + 1);
		if (v & 1)
		{
			for (j = 0; j < len; ++j)
				g->seq[i][j] = gfa_comp_table[(uint8_t)s->seq[len - j - 1]];
		}
		else
			memcpy(g->seq[i], s->seq, len);
		g->seq[i][len] = 0; // null terminated for convenience
		for (j = 0; j < (uint32_t)sub->v[i].n; ++j)
		{
			uint64_t a = sub->a[sub->v[i].off + j];
			g->arc[k].a = (uint64_t)i << 32 | a >> 32;
			g->arc[k].o = gfa->arc[(uint32_t)a].ow;
			++k;
		}
		assert(k <= (int32_t)g->n_arc);
	}

	// Here sub is no more useful, so it should be freed
	gwf_sub_free(sub);

	return g;
}

void gwf_free(gwf_graph_t *g)
{
	int32_t i;
	for (i = 0; i < (int32_t)g->n_vtx; ++i)
		free(g->seq[i]);
	free(g->len);
	free(g->seq);
	free(g->arc);
	free(g->src);
	free(g);
}

void gwf_graph_print(FILE *fp, const gwf_graph_t *g)
{
	uint32_t i;
	for (i = 0; i < g->n_vtx; ++i)
		fprintf(fp, "S\t%d\t%s\tLN:i:%d\n", i, g->seq[i], g->len[i]);
	for (i = 0; i < g->n_arc; ++i)
		fprintf(fp, "L\t%d\t+\t%d\t+\t%dM\n", (uint32_t)(g->arc[i].a >> 32), (uint32_t)g->arc[i].a, g->arc[i].o);
}

int main(int argc, char *argv[]) // kstring_t name, comment, seq, qual;
{
	gzFile fp;				  // file pointer
	kseq_t *ks;				  // fasta's reads, some of the fields of the kseq_t type are "kstring_t name, comment, seq, qual"
	ketopt_t o = KETOPT_INIT; // command-line argument parser
	gfa_t *gfa;				  // reference graph
	gwf_graph_t *g;			  // graph representation for the algorithm
	gwf_path_t path;		  // path in the graph
	int c, print_graph = 0;
	uint32_t v0 = 0 << 1 | 0; // first segment, forward strand
	uint32_t max_lag = 0;	  // max lag behind the furthest wavefront (related to pruning)
	void *km = 0;			  // chunk of memory managed with kalloc, see "kalloc.c"
	char *sname = 0;		  // starting segment name

	FILE *fGAF; // output GAF file
	char filename[32];
	struct tm *timenow;
	time_t now = time(NULL);
	timenow = gmtime(&now); // timestamp (UTC)
	strftime(filename, sizeof(filename), "gwfa-tb_%Y-%m-%d_%H:%M:%S.gaf", timenow);
	fGAF = fopen(filename, "w");
	if (fGAF == NULL)
	{
		fprintf(stderr, "Error opening output GAF file.\n");
		abort();
	}

	while ((c = ketopt(&o, argc, argv, 1, "ptl:s:", 0)) >= 0)
	{ // command-line arguments management
		if (c == 'p')
			print_graph = 1;
		else if (c == 'l')
			max_lag = atoi(o.arg);
		else if (c == 's')
			sname = o.arg;
	}
	if ((!print_graph && argc - o.ind < 2) || (print_graph && argc == o.ind)) // missing arguments
	{
		fprintf(stderr, "Usage: gwf-tb [options] <target.gfa|fa> <query.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -l INT    max lag behind the furthest wavefront; 0 to disable [0]\n");
		fprintf(stderr, "  -s STR    starting segment name [first]\n");
		fprintf(stderr, "  -p        output GFA in the forward strand\n");
		return 1;
	}

	km = km_init();

	gfa = gfa_read(argv[o.ind]); // input reference graph reading
	assert(gfa);
	if (sname)
	{								   // if a starting segment name has been provided
		int32_t sid;				   // segment ID
		sid = gfa_name2id(gfa, sname); // returns the ID of the segment starting from the name
		if (sid < 0)
			fprintf(stderr, "ERROR: failed to find segment '%s'\n", sname);
		else
			v0 = sid << 1 | 0; // TODO: also allow to change the orientation
	}
	g = gwf_gfa2gwf(gfa, v0); // convert the input gfa reference to the algorithm's internal representation
	if (print_graph)
	{
		gwf_graph_print(stdout, g); // graph printing to stdout
		return 0;					// free memory
	}
	gwf_ed_index(km, g); // indexing the graph to easily access vertices' neighbors, see gwf-ed.c

	fp = gzopen(argv[o.ind + 1], "r"); // open file for reading
	assert(fp);						   // check file has opened correctly
	ks = kseq_init(fp);				   // parse fasta's reads

	clock_t before, after;
	before = clock();
	while (kseq_read(ks) >= 0) // while the file still contains reads
	{
		int32_t abl; // alignment block length (total number of '='s, 'X's, 'D's, and 'I's in the alignment)
		string cig = "";
		gwf_ed(km, g, ks->seq.l, ks->seq.s, 0, -1, max_lag, &path, cig, abl);

		// GAF preparation
		int32_t i, last_len = -1, len = 0;
		fprintf(fGAF, "%s\t%ld\t0\t%ld\t+\t", ks->name.s, ks->seq.l, ks->seq.l);
		for (i = 0; i < path.nv; ++i)
		{
			uint32_t v = g->src[path.v[i]];
			fprintf(fGAF, "%c%s", "><"[v & 1], gfa -> seg[v >> 1].name);
			last_len = gfa->seg[v >> 1].len;
			len += last_len;
		}
		fprintf(fGAF, "\t%d\t0\t%d\t%d\t%d\t60\tcg:Z:%s\n", len, len - (last_len - path.end_off) + 1, path.s, abl, cig.c_str());
	}
	after = clock();
	cout << "Total alignment time: " << (double)(after - before) / CLOCKS_PER_SEC << " s" << endl;
	kseq_destroy(ks);
	gzclose(fp);
	gfa_destroy(gfa);

	gwf_cleanup(km, g);
	gwf_free(g);
	km_destroy(km);

	fclose(fGAF);
	return 0;
}
