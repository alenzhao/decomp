#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

/* #include "mt_rand.h" */
#include "parse_option.h"

/**************************************************************
  To accept a weighted (undirected) graph, and to detect possibly
  overlapping communities, by means of finding semi-cliques and patch
  up the rest.

***********************************************************************/
static void* Malloc_buf(size_t s, char* purpose) {
  void* r = NULL;
  r = malloc(s);
  if(r == NULL) {
    fprintf(stderr, "Error: Cannot allocate %d bytes for %s.\n", s, purpose);
    exit(1);
  }
  return r;
}
#define Malloc(n,t,p) ((t *)Malloc_buf(sizeof(t)*(n), p))

/****************************************************************/
/* only positive weights are included. edges with zero or negative
   weights are ignored.  */
typedef struct s_edge {
  int from, to; /* 0-based indices */
  int idx; /* id of component it belongs at the moment */
  double weight;
  double betweeness; /* basically the number of shortest paths between
			vertices that pass through this edge */
} edge;

static int n_edges=0, n_ignored_edges=0;
static int n_vertices=0; /* 1 + the maximum index seen in the non-ignored edges */
static int n_capacity_edges=0;
static edge* all_edges=NULL;

static void add_edge(int from, int to, double weight) {
  if(weight <= 0 || from==to) { /* ignore self loops */
    n_ignored_edges++;
    return;
  }

  if(n_edges >= n_capacity_edges) { /* expand buffer */
    n_capacity_edges = n_capacity_edges < 10 ? 10 : 2*n_capacity_edges;
    all_edges = (edge*)realloc(all_edges, n_capacity_edges*sizeof(edge));
    if(all_edges == NULL) {
      fprintf(stderr, "Not enough memory for %d edges in add_edge().\n", n_capacity_edges);
      exit(-1);
    }
  }

  all_edges[n_edges].from = from;
  all_edges[n_edges].to = to;
  all_edges[n_edges].weight = weight;
  n_edges++;

  if(from+1 > n_vertices) n_vertices = from+1;
  if(to+1 > n_vertices) n_vertices = to+1;
}

static void clear_all_edges() {
  if(all_edges) free(all_edges);
  all_edges = NULL;
  n_capacity_edges = 0;
  n_edges = 0;
  n_ignored_edges = 0;
}

static int read_word(FILE* f, char* buf, int L) {
  /* read next item and until next whitespace, read at most L-1
     characters, and put '\0' at the end.

     return EOF if end of file reached before reading any thing
     interesting.
     return 0 if end of line reached before reading other things.
     Otherwise returns the number of characters read.
  */
  int i=0, c=0;
  while((c=fgetc(f))==' ' || c=='\t' || c=='\r')
    ;
  if(c==EOF) return EOF;
  if(c=='\n') return 0;

  buf[i] = c;
  for(i=1; i<L-1; i++) {
    c = fgetc(f);
    if(c==' ' || c=='\t' || c=='\r' || c==EOF)
      break;
    if(c=='\n') {ungetc(c,f); break;}
    buf[i] = c;
  }
  buf[i] = '\0';

  /*
  printf("==== debug: read_word(): \"%s\"\n", buf);
  */

  return i;
}
static int read_one_edge(FILE* f, int* from, int* to, double* weight) {
  /* either
       from to weight
     or
       From: from To: to Coef: weight
       From: from To: to Effect: weight

       **** trying
       From: from To: to Score: weight

     If use the second form, the order of from, to and coef (effect)
     can be different, other things on the same line are ignored.

     return EOF if reached end of file before reading anything interesting.
     otherwise returns the number of items read, which on success should be 3.
 */
  int r=0, ifrom=0, ito=0, iweight=0;
#define MAX_BUF 64
  static char buf[MAX_BUF]; /* should be large enough for our use */
  r = read_word(f, buf, MAX_BUF);
  if(r == EOF) return EOF;
  if(r == 0) return 0; /* read nothing in this line */

  if('0'<=buf[0] && buf[0]<='9') { /* simple form */
    *from = atoi(buf);
    if((r = read_word(f, buf, MAX_BUF)) <= 0) return 1;
    *to = atoi(buf);
    if((r = read_word(f, buf, MAX_BUF)) <= 0) return 2;
    *weight = atof(buf);
    if((r = read_word(f, buf, MAX_BUF)) > 0) { /* eliminate the end of line */
      printf("Unexpected item after reading 'from', 'to' and 'weight': \"%s\"\n", buf);
    }
    return 3;
  }

  /* keyword form */
  while(1) {
    for(r=0; buf[r]!='\0'; r++) buf[r] = tolower(buf[r]);
    if(0==strcmp(buf, "from:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *from = atoi(buf);
      ifrom = 1;
    } else if(0==strcmp(buf, "to:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *to = atoi(buf);
      ito = 1;
    } else if(0==strcmp(buf, "coef:") || 0==strcmp(buf, "effect:") || 0==strcmp(buf, "score:")) {
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
      *weight = atof(buf);
      iweight = 1;
    } else { /* ignore other keywords and its associated value */
      if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
    }
    /* possibly next keyword */
    if((r = read_word(f, buf, MAX_BUF)) <= 0) break;
  }

  return ifrom+ito+iweight;
}
static void read_edges_from_file(char* name) {
  FILE* f=NULL;
  int from=0, to=0, n=0, i=0;
  double weight=0;

  f = fopen(name, "r");
  if(f == NULL) {
    fprintf(stderr, "Fail to read \"%s\" for edges.\n", name);
    exit(-1);
  }

  i=0;
  while(!feof(f)) {
    n = read_one_edge(f, &from, &to, &weight);
    if(n==EOF || (n == 0 && feof(f))) {
      break;
    } else if(n==0) {
      continue;
    } else if(n != 3) {
      fprintf(stderr, "Error in reading edge file, after reading %d edge(s).\nIn the edge file, each line represents an edge, with the 0-based 'from' index, followed by 0-based 'to' index, followed by the weight, all separated by space.\n", i);
      exit(-1);
    }
    add_edge(from, to, fabs(weight));
    i++;
  }

  fclose(f);
}

static void print_all_edges(FILE* f) { /* mainly for debug */
  int i=0;
  fprintf(f, "=========================================\n");
  fprintf(f, "%d (positive weight) edges, %d ignored:\n", n_edges, n_ignored_edges);
  for(i=0; i<n_edges; i++)
    fprintf(f, "%d:\tfrom: %d\tto: %d\tweight: %f\n", i, all_edges[i].from, all_edges[i].to, all_edges[i].weight);
  fprintf(f, "=========================================\n");
}

/***************************************************************/
/* to index the edges (0-based ids in all_edges) of each vertex */
typedef struct s_vertex_neis {
  int idx; /* for use during index construction. Later used to
	      indicate components in the graph when calculating
	      betweeness. */
  int n_parents;
  int attribute; /* use for temporary use */
  int* parents; /* not own this, but point into a long buffer.
		     vertice indices of parents sorted ascendingly */
  /* used in calculating the shortest path distance between pairs of vertices */
  int n_edges;
  int* edges; /* not own this, but point into a long buffer */
  int distance; /* -1 for unknown yet, >0 for the distance to a
                    particular vertex under consideration */
  int predecessor; /* for the common case of 1 predecessor. >=0
		      indicates the index of the predecessor edge in the
		      path. -1 indicates unknown yet.  <= -2 indicates
		      there are more than 1 predecessor, and the
		      absolute value is the number of predecessors,
		      need to look at the edges, those neighbors with
		      distance one less than the current distance is
		      one of the predecessors. */
  double bk; /* initially 1, to be propagated to count the paths to
		the given vertex, and the edges */
  double betweeness; /* the vertex betweeness, basically is the number
			of shortest paths of all pairs of vertices
			that pass through this vertex. */
} vertex_neis;

static int* vertex_edges_buf=NULL; /* length 2*n_edges, since each edge is incident to two vertices */
static int* vertex_neis_buf=NULL; /* length 2*n_edges, since each edge is incident to two vertices */
static vertex_neis* all_vertices=NULL;

static void clear_vertices() {
  if(all_vertices) {free(all_vertices); all_vertices = NULL;}
  if(vertex_edges_buf) {free(vertex_edges_buf); vertex_edges_buf = NULL;}
  if(vertex_neis_buf) {free(vertex_neis_buf); vertex_neis_buf = NULL;}
}

int remove_dups(int n, int v[]) {
  /* assume v is sorted, remove duplicates in v in place, returns the
     size of v after removal. */
  int i=0,j=0; /* j points to the slot to be filled next */
  if(n <= 0) return 0;
  /* the first is always kept */
  for(i=1,j=1; i<n; i++) {
    if(v[i]!=v[i-1]) { /* keep */
      if(i!=j) v[j] = v[i];
      j++;
    }
  }
  return j;
}

int ascending_int(const void* a, const void* b) {
  return (*(int*)a) - (*(int*)b);
}
int ascending_vertex_idx(const void* a, const void* b) {
  /* *a and *b contain vertex indices */
  int ia=0,ib=0, r=0;
  ia = (*(int*)a);
  ib = (*(int*)b);
  r = all_vertices[ia].idx - all_vertices[ib].idx;
  if(r != 0) return r;
  return ia - ib;
}
int ascending_edge_idx(const void* a, const void* b) {
  /* *a and *b contain edge indices */
  int ia=0,ib=0;
  ia = (*(int*)a);
  ib = (*(int*)b);
  return all_edges[ia].idx - all_edges[ib].idx;
}

int make_set(int n, int v[]) {
  /* v has n integers, sort it and remove duplicates, return the size
     of v afterwards. */
  qsort(v, n, sizeof(int), ascending_int);
  /* remove duplicates */
  return remove_dups(n,v);
}

static void index_vertices() {
  int i=0, *p=NULL, *q=NULL;
  vertex_neis* ve=NULL;

  clear_vertices();
  vertex_neis_buf = Malloc(2*n_edges, int, "index_vertices()");
  vertex_edges_buf = Malloc(2*n_edges, int, "index_vertices()");
  all_vertices = Malloc(n_vertices, vertex_neis, "index_vertices()");
  memset(all_vertices, 0, sizeof(vertex_neis)*n_vertices);

  /* count the edges for each vertex first */
  for(i=0; i<n_edges; i++) {
    all_vertices[all_edges[i].to].n_parents++;

    all_vertices[all_edges[i].from].n_edges++;
    all_vertices[all_edges[i].to].n_edges++;
  }
  /* point into the buffer accordingly */
  p = vertex_neis_buf;
  q = vertex_edges_buf;
  for(i=0; i<n_vertices; i++) {
    all_vertices[i].idx = 0;
    all_vertices[i].attribute = 0; /* to use as idx, but for parents */

    all_vertices[i].parents = p;
    p += all_vertices[i].n_parents;

    all_vertices[i].edges = q;
    q += all_vertices[i].n_edges;
  }
  /* fill in the indices of the edges and the neighbors */
  for(i=0; i<n_edges; i++) {
    ve = all_vertices + all_edges[i].from;
    ve->edges[ve->idx] = i;
    ve->idx++;

    ve = all_vertices + all_edges[i].to;
    ve->edges[ve->idx] = i;
    ve->parents[ve->attribute++] = all_edges[i].from;
    ve->idx++;
  }
  /* sort the neighbors ascendingly */
  for(i=0; i<n_vertices; i++) {
    all_vertices[i].n_parents = make_set(all_vertices[i].n_parents,
					 all_vertices[i].parents);
  }
}

void print_set(FILE* f, int n, int v[]) {
  int i=0;
  for(i=0; i<n; i++)
    fprintf(f, " %d", v[i]);
}

static void print_all_vertices(FILE* f) { /* mainly for debug */
  int i=0;
  fprintf(f, "=========================================\n");
  fprintf(f, "%d vertices:\n", n_vertices);
  for(i=0; i<n_vertices; i++) {
    fprintf(f, "%d:\t%d parents:", i, all_vertices[i].n_parents);
    print_set(f, all_vertices[i].n_parents, all_vertices[i].parents);
    fprintf(f,"\n");
  }
  fprintf(f, "=========================================\n");
}

/* some utilities to operate on sets, which are represented as sorted (ascendingly) array */

int is_member(int x, int n, int v[]) {
  /* return 1 if x is in v, 0 otherwise */
  int i=0;
  for(i=0; i<n; i++) {
    if(x == v[i]) return 1;
    if(x < v[i]) return 0;
  }
  return 0;
}

int skip_union(int n1, int v1[], int n2, int v2[], int n_to_skip, int to_skip[], int out[]) {
  /* out = (v1 + v2) - to_skip.
     return the size of the resulting union
   */
  int i=0, x=0;
  while(n1>0 && n2>0) {
    if(v1[0] < v2[0]) {
      x = v1[0];
      n1--; v1++;
    } else if(v1[0] == v2[0]) {
      x = v1[0];
      n1--; v1++;
      n2--; v2++;
    } else {
      x = v2[0];
      n2--; v2++;
    }

    if(!is_member(x, n_to_skip, to_skip)) {out[i++] = x;}
  }
  /* for anything that remains, only one of v1 or v2 may have remains */
  for(; n1>0; n1--,v1++) {
    if(!is_member(v1[0], n_to_skip, to_skip)) {out[i++] = v1[0];}
  }
  for(; n2>0; n2--,v2++) {
    if(!is_member(v2[0], n_to_skip, to_skip)) {out[i++] = v2[0];}
  }

  return i;
}

/***************************************************************/
/* calculate the vertex and edge betweenness using the algorithm as described in
   Newman, M. E. J. (2001) Phys. Rev. E 64, 016132.
 */

static void add_predecessor(int v, int e) {
  int p=0;
  p = all_vertices[v].predecessor;
  if(p == -1) { /* no predecessor yet */
    all_vertices[v].predecessor = e;
  } else if(p >= 0) { /* now has two predecessors */
    all_vertices[v].predecessor = -2;
  } else { /* <= -2, already more than 1 predecessor */
    all_vertices[v].predecessor--;
  }
}

int cal_distance_to_vertex(int v, int buf[]) {
  /* calculate the distance to v, modify the distance and predecessor
     in all_vertices.

     return the number of vertices with distance determined.

     buf is used as the queue in BFS, and also contains the indices of
     the vertices, in increasing order of distance, so buf should have
     length n_vertices.
   */
  int i=0, j=0, d=0, p=0, k=0;
  int w=0, z=0;
  int n=0, *es=NULL, e=0;

  for(i=0; i<n_vertices; i++) {
    all_vertices[i].distance = -1; /* -1 for unknown yet */
    all_vertices[i].predecessor = -1; /* -1 for unknown yet */
  }

  all_vertices[v].distance = 0;

  /* use simple breadth first search, because we consider each edge as
     having the same length. buf is used as the queue, and those
     already in the queue have their distance determined, but need to
     process their neighbors. */
  buf[0] = v;
  k = 1;
  p = 0;
  while(p < k) {
    w = buf[p++];
    d = all_vertices[w].distance; /* normally should be >= 0 */
    n = all_vertices[w].n_edges;
    es = all_vertices[w].edges;
    /* go through the edges instead of vertices, because we record as edge as predecessor */
    for(j=0; j<n; j++) {
      e = es[j];
      z = all_edges[e].from == w ? all_edges[e].to : all_edges[e].from;

      if(all_vertices[z].distance < 0) { /* newly encountered */
	all_vertices[z].distance = d+1;
	buf[k++] = z;
	add_predecessor(z, e);
      } else if(all_vertices[z].distance == d+1) { /* already encountered */
	add_predecessor(z, e);
      }
    }
  }
  return k;
}

int calc_betweeness(int nv, int vs[], int ne, int edges[], int buf[]) {
  /* vs contains nv indices to the vertices in one or more components,
     es contains ne indices to the edges of the component(s).

     calculate the vertex and edge betweeness, modifying all_vertices,
     and all_edges accordingly. buf is used as temporary use, also on
     return stores the indices of the edges in descreasing order of
     betweeness. Returns the number of components in the graph. */

  int i=0, j=0, k=0, d=0, v=0, w=0, z=0, n=0, np=0;
  int m=0, *es=NULL, e=0;
  int n_components=0;

  for(i=0; i<nv; i++) {
    v = vs[i];
    all_vertices[v].betweeness = 0;
    all_vertices[v].idx = 0; /* to indicate the 1-based component id,
				0 for unassigned yet */
  }
  for(i=0; i<ne; i++) all_edges[edges[i]].betweeness = 0;

  for(i=0; i<nv; i++) {
    v = vs[i];
    for(j=0; j<nv; j++) all_vertices[vs[j]].bk = 1;

    n = cal_distance_to_vertex(v, buf);

    for(j=n-1; j>=0; j--) {
      /* from the furthest vertex, propagate the bk */
      w = buf[j];
      if(all_vertices[w].predecessor >= 0) { /* only one */
	e = all_vertices[w].predecessor;
	z = all_edges[e].from == w ? all_edges[e].to : all_edges[e].from;
	all_vertices[z].bk += all_vertices[w].bk;
	all_edges[e].betweeness += all_vertices[w].bk;
      } else if(all_vertices[w].predecessor <= -2) { /* more than one */
	/* look through the neighbors to find the predecessors,
	   share the bk equally
	 */
	m = all_vertices[w].n_edges;
	es = all_vertices[w].edges;
	d = all_vertices[w].distance;
	/* the absolute value is the number of predecessors */
	np = -all_vertices[w].predecessor;
	for(k=0; k<m; k++) {
	  e = es[k];
	  z = all_edges[e].from == w ? all_edges[e].to : all_edges[e].from;
	  if(all_vertices[z].distance == d-1) { /* predecessor */
	    all_vertices[z].bk += all_vertices[w].bk/np;
	    all_edges[e].betweeness += all_vertices[w].bk/np;
	  }
	}
      }/* else no predecessor */
    }

    /* accumulate betweeness, only for the reachable vertices */
    for(j=n-1; j>=0; j--) {
      w = buf[j];
      all_vertices[w].betweeness += all_vertices[w].bk;
    }

    /* see if it is new component, now the order in buf can be modified */
    if(all_vertices[v].idx <= 0) { /* new component, give id */
      ++n_components;
      for(j=0; j<n; j++) {
	w = buf[j];
	all_vertices[w].idx = n_components;
      }
    }
  }

  /* the component of an edge is the component of either of its
     vertex. */
  for(i=0; i<ne; i++) {
    e = edges[i];
    all_edges[e].idx = all_vertices[all_edges[e].from].idx;
  }

  return n_components;
}

static void cut_edge_from_vertex(int e, int v) {
  int i=0, n=0, *es=NULL;
  n = all_vertices[v].n_edges;
  es = all_vertices[v].edges;
  for(i=n-1; i>=0; i--)
    if(es[i] == e) break;
  if(i >= 0) { /* found */
    es[i] = es[n-1];
    all_vertices[v].n_edges--;
  }
}
static void cut_edge(int e) {
  cut_edge_from_vertex(e, all_edges[e].from);
  cut_edge_from_vertex(e, all_edges[e].to);
}

int highest_edge_betweeness_idx(int ne, int es[]) {
  /* return the index into es with the highest edge betweeness */
  int i=0, idx=0, e=0;
  double mb=0;

  if(ne <= 0) return 0;

  mb = all_edges[es[0]].betweeness;
  for(i=1; i<ne; i++) {
    e = es[i];
    if(all_edges[e].betweeness > mb) {
      mb = all_edges[e].betweeness;
      idx = i;
    }
  }
  return idx;
}

int divide_one_component(int nv, int vs[], int ne, int es[], int buf[]) {
  /* vs contains nv indices to the vertices in one or more components,
     es contains ne indices to the edges of the component(s).

     successively remove the edge with highest betweeness until there
     are more than one component.  Modify vs and es in place, so that
     they are in increasing components.

     return the number of edges cut, they are placed at the end of es.

     buf is for temporary use, should be large enough to hold a copy
     of vs.
     all_vertices and all_edges will be modified in the process.

     The assigned component id are for identifying the componets
     within vs only.
   */
  int n_components=0, i=0, e=0, n_cut=0;

  while((n_components = calc_betweeness(nv,vs, ne,es, buf)) <= 1) {
    i = highest_edge_betweeness_idx(ne, es);
    /* cut the edge with highest betweeness, in es[i] */
    e = es[i];
    printf("*** Cut edge %d: %d -> %d\tbetweeness: %f\n", e,
	   all_edges[e].from, all_edges[e].to, all_edges[e].betweeness);
    cut_edge(e);

    es[i] = es[ne-1];
    n_cut++;
    ne--;
  }
  /* sort vs and es by components */
  qsort(vs, nv, sizeof(int), ascending_vertex_idx);
  /* an edge is associated with the component of either of its vertices */
  qsort(es, ne, sizeof(int), ascending_edge_idx);

  return n_cut;
}

void finalize_component(int id, int nv, int vs[], int buf1[], int buf2[]) {
  /* v contains the n vertex indices of the component with the id.
     print the component and its parents.
   */
  int i=0, v=0, *pbuf=NULL, n=0;

  /* the disjoint components */
  printf("= component %d, size %d:", id, nv);
  print_set(stdout, nv,vs);
  printf("\n");

  /* take also the parents as the partition */
  pbuf = buf1; /* use buf1 and buf2 alternatively */
  n = 0;
  for(i=0; i<nv; i++) {
    v = vs[i];
    all_vertices[v].idx = id; /* also finalize the component id */
    n = skip_union(all_vertices[v].n_parents, all_vertices[v].parents,
		   n, pbuf, nv,vs,  pbuf==buf1 ? buf2 : buf1);
    pbuf = pbuf==buf1 ? buf2 : buf1;
  }

  /* the component with its parents, in one partition */
  printf("= partition %d, size %d:", id, nv + n);
  print_set(stdout, nv,vs);
  print_set(stdout, n,pbuf);
  printf("\n");

  /* the component with its parents, in one partition, but the parents
     are marked */
  printf("= xpartition %d, size %d+%d=%d:", id, nv,n, nv + n);
  print_set(stdout, nv,vs);
  printf(" p");
  print_set(stdout, n,pbuf);

  printf("\n\n");
}

int divide_all_components(int nv, int vs[], int ne, int es[],
			  int stop_size, int id_component,
			  int buf1[], int buf2[]) {
  /* recursively divide the components until stopping criteria are
     matched.

     id_component is the next component id to be used for finalized
     component.

     stop if a component is not larger than stop_size.

     return the updated component id.
  */
  int n_cuts=0, new_ne=0, v_cs=0, v_ce=0, e_cs=0, e_ce=0;
  int idx=0;

  if(nv <= stop_size) { /* finalize a component */
    finalize_component(id_component, nv,vs, buf1,buf2);
    return id_component+1;
  }

  n_cuts = divide_one_component(nv, vs, ne, es, buf1);
  new_ne = ne - n_cuts;

  printf("===== Cut %d edges out of %d, ratio: %f\n", n_cuts, ne, n_cuts/(double)ne);

  /* go through the components, recursively divide if needed */
  for(v_cs=0, e_cs=0; v_cs<nv; v_cs=v_ce, e_cs=e_ce) {
    idx = all_vertices[vs[v_cs]].idx;
    for(v_ce=v_cs+1; v_ce<nv; v_ce++) {
      if(all_vertices[vs[v_ce]].idx != idx) break;
    }
    /* [v_cs, v_ce) contains the range of the vertices of one component */

    for(e_ce=e_cs+1; e_ce<new_ne; e_ce++) {
      if(all_edges[es[e_ce]].idx != idx) break;
    }
    /* [e_cs, e_ce) contains the range of the edges of one component */

    id_component = divide_all_components(v_ce-v_cs, vs+v_cs, e_ce-e_cs, es+e_cs,
					 stop_size, id_component, buf1,buf2);
  }
  return id_component;
}

void cut_clustering(int stop_size) {
  int i=0;
  int *buf1=NULL, *buf2=NULL, *vbuf=NULL, *ebuf=NULL;
  buf1 = Malloc(3*n_vertices, int, "cut_clustering()");
  buf2 = buf1 + n_vertices;
  vbuf = buf2 + n_vertices;
  ebuf = Malloc(n_edges, int, "cut_clustering()");

  for(i=0; i<n_vertices; i++) vbuf[i] = i;
  for(i=0; i<n_edges; i++) ebuf[i] = i;

  divide_all_components(n_vertices, vbuf, n_edges, ebuf,
			stop_size, 1, buf1,buf2);

  free(buf1);
  free(ebuf);
}

/***************************************************************/
#define DEFAULT_STOP_SIZE      60

#define xstr(s) str(s)
#define str(s) #s

/* The option structure of this program */
option all_options[] = {
  /* str,type(STRING),is_required(0),description,index(0),end_index(0), val_int(0),val_real(0),val_str(NULL) */
  {"-?",      FLAG,   0,  "Showing the usage.\n", 0,0,  0,  0.0,    NULL},

  {"-f",   STRING,    1,  "File name of the input edges. Each line in the file represents an edge, which consists of a 0-based 'from' index, a 0-based 'to' index, and a weight, separated by space..\n", 0,0,  -1,  0.0,    NULL},
  {"-s",   INT,       0,  "Threshold of component size. If a component is no larger than this, it is not further divided. Default " xstr(DEFAULT_STOP_SIZE) ".\n", 0,0,  DEFAULT_STOP_SIZE,  0.0,    NULL},

  {NULL,      STRING, 0,  NULL, 0,0,  0,  0.0,    NULL} /* end marker */
};
/***************************************************************/

int main(int argc, char *argv[])
{
  time_t starting_time, end_time;
  char* edge_file_name=NULL;
  int stop_size = DEFAULT_STOP_SIZE;

  /* do some initialization */
  if(parse_options(argc,argv,all_options)) {
    /* error parsing options */
    usage(stderr, argv[0], all_options);
    exit(0);
  }
  if(get_flag_option_value("-?",all_options,0))
    usage(stdout, argv[0], all_options);

  /*********************/
  time(&starting_time);
  printf("starting time: %ld\t%s", (long)starting_time, ctime(&starting_time));

  /*********************/
  stop_size = get_int_option_value("-s",all_options, stop_size);
  printf("Stop Component Size: %d\n", stop_size);

  edge_file_name = get_str_option_value("-f",all_options,NULL);
  printf("Input edges file name: %s\n", edge_file_name);
  read_edges_from_file(edge_file_name);
  printf("Read %d edges for %d vertices, ignored %d of them.\n", n_edges+n_ignored_edges, n_vertices, n_ignored_edges);
  print_all_edges(stdout);

  index_vertices();

  print_all_vertices(stdout);

  cut_clustering(stop_size);

  /*********************/

  time(&end_time);
  printf("ending time: %d\t%s", (long)end_time, ctime(&end_time));
  printf("Total seconds used: %d\n", end_time - starting_time);
  printf("====End\n\n");
  /* cleanup */
  clear_all_edges();
  clear_vertices();

  return 0;
}
