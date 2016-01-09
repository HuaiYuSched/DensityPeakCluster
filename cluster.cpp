#include <cstdio>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#define MAX_NODES 4096
#define MAX_LINES (MAX_NODES *MAX_NODES)

using namespace std;

struct Edge {
  int i, j;
  double distance;
};

struct Vertex {
  double x, y;
  double *dists;
};
Vertex vexs[MAX_NODES];

int maxnode;
double **dist;
int nlines;
int edgenum;
double *ro;
double *delta;
int *neigh;
int *roindex;

double dc;

int *cl;
int *icl;
int *halo;
int ncluster;
Edge edges[MAX_LINES];

inline bool cmpedge(const Edge &a, const Edge &b) {
  return a.distance < b.distance;
}
inline double caldist(Vertex a, Vertex b) {
  return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
void loaddata() {
  char file[64];

  cout << "Please Input the data file" << endl;
  cout << "The format of the data file is:" << endl;
  cout << "First column: x of node i" << endl;
  cout << "Second column: y of node j" << endl;

  cin >> file;
  FILE *filein = fopen(file, "r");
  while (filein == NULL) {
    cout << "Error filename" << endl;
    cin >> file;
    filein = fopen(file, "r");
  }

  nlines = 0;
  int meanless;
  while (fscanf(filein, "%lf %lf\t%d", &vexs[nlines].x, &vexs[nlines].y,
                &meanless) != EOF) {
    nlines++;
    if (nlines > MAX_NODES) {
      cout << "Error, too many nodes" << endl;
      exit(0);
    }
  }
  edgenum = 0;
  for (int i = 0; i < nlines; i++) {
    vexs[i].dists = new double[nlines];
    for (int j = 0; j < nlines; j++) {
      vexs[i].dists[j] = caldist(vexs[i], vexs[j]);
      if (j > i) {
        // fprintf(fp, "%d %d %6.2f\n", i+1, j+1, vexs[i].dists[j]);

        edges[edgenum].i = i;
        edges[edgenum].j = j;
        edges[edgenum].distance = vexs[i].dists[j];
        edgenum++;
      }
    }
  }

  fclose(filein);
  maxnode = nlines;
  nlines = edgenum;

  dist = new double *[maxnode];

  for (int i = 0; i < maxnode; i++) {
    dist[i] = new double[maxnode];
  }
  for (int i = 0; i < nlines; ++i) {
    dist[edges[i].i][edges[i].j] = dist[edges[i].j][edges[i].i] =
        edges[i].distance;
  }
}

void init() {
  double percent = 2;
  cout << "average percentage of neighbours (hard coded): " << percent << endl;

  int cutdown = round(edgenum * percent / 100);
  sort(edges, edges + edgenum, cmpedge);
  dc = edges[cutdown].distance;
  cout << "The dc is " << dc << endl;

  for (int i = 0; i < maxnode - 1; i++) {
    for (int j = i + 1; j < maxnode; j++) {
      double inc = exp(-pow(dist[i][j] / dc, 2));
      ro[i] += inc;
      ro[j] += inc;
    }
  }
}

inline bool cmpro(const int &a, const int &b) { return ro[a] < ro[b]; }
void finddist() {
  double maxdist = edges[nlines - 1].distance;
  roindex = new int[maxnode];
  for (int i = 0; i < maxnode; ++i) {
    roindex[i] = i;
  }
  sort(roindex, roindex + maxnode, cmpro);

  const int maxroindex = roindex[maxnode - 1];

  delta[maxroindex] = 0;
  neigh[maxroindex] = -1;
  for (int i = 0; i < maxnode - 1; ++i) {
    const int &id = roindex[i];
    delta[id] = maxdist;
    neigh[id] = -1;
    for (int j = i + 1; j < maxnode; ++j) {
      if (dist[id][roindex[j]] <= delta[id]) {
        delta[id] = dist[id][roindex[j]];
        neigh[id] = roindex[j];
      }
    }

    delta[maxroindex] = max(delta[id], delta[maxroindex]);
  }
  cout << "Generated file:Decision_graph" << endl;
  cout << "column 1:Density" << endl;
  cout << "column 2:Delta" << endl;
  FILE *fp = fopen("Decision_graph", "w");
  for (int i = 0; i < maxnode; ++i) {
    fprintf(fp, "%6.2f %6.2f\n", ro[i], delta[i]);
  }
  fclose(fp);
}

void cluster() {
  cout << "Please enter minimum of ro" << endl;
  double romin;
  cin >> romin;
  cout << "Please enter minimum of delta" << endl;
  double deltamin;
  cin >> deltamin;

  ncluster = 0;
  for (int i = 0; i < maxnode; ++i) {
    if ((ro[i] > romin) && (delta[i] > deltamin)) {
      cl[i] = ncluster++;
    } else {
      cl[i] = -1;
    }
  }
  cout << "NUMBER OF CLUSTERS: " << ncluster << endl;

  icl = new int[ncluster];
  for (int i = 0; i < maxnode; ++i) {
    if (cl[i] != -1) {
      icl[cl[i]] = i;
    }
  }

  cout << "Performing assignation" << endl;
  for (int i = maxnode - 1; i >= 0; --i) {
    if (cl[roindex[i]] == -1) {
      cl[roindex[i]] = cl[neigh[roindex[i]]];
    }
  }

  for (int i = 0; i < maxnode; ++i) {
    halo[i] = cl[i];
  }
  double *bordro = new double[ncluster];
  for (int i = 0; i < ncluster; ++i) {
    bordro[i] = 0;
  }

  for (int i = 0; i < maxnode - 1; ++i) {
    for (int j = i + 1; j < maxnode; ++j) {
      if ((cl[i] != cl[j]) && (dist[i][j] <= dc)) {
        double roAver = (ro[i] + ro[j]) / 2;
        bordro[cl[i]] = max(roAver, bordro[cl[i]]);
        bordro[cl[j]] = max(roAver, bordro[cl[j]]);
      }
    }
  }

  for (int i = 0; i < maxnode; ++i) {
    if (ro[i] < bordro[cl[i]]) {
      halo[i] = -1;
    }
  }
  if (ncluster != 0)
    // delete[] bordro;

    for (int i = 0; i < ncluster; ++i) {
      int nc = 0;
      int nh = 0;
      for (int j = 0; j < maxnode; ++j) {
        if (cl[j] == i) {
          ++nc;
        }
        if (halo[j] == i) {
          ++nh;
        }
      }
      cout << "CLUSTER: " << i + 1 << " CENTER: " << icl[i] + 1
           << " ELEMENTS: " << nc << " CORE: " << nh << " HALO: " << nc - nh
           << endl;
    }

  FILE *fp = fopen("Clusterfile", "w");
  cout << "Generated file:Clusterfile" << endl;
  cout << "column 1:element id" << endl;
  cout << "column 2:cluster assignation without halo control" << endl;
  cout << "column 3:cluster assignation with halo control" << endl;
  for (int i = 0; i < maxnode; i++) {
    fprintf(fp, "%lf %lf %d %d\n", vexs[i].x, vexs[i].y, cl[i] + 1,
            halo[i] + 1);
  }
  fclose(fp);
}

int main() {
  loaddata();
  ro = new double[maxnode];
  delta = new double[maxnode];
  neigh = new int[maxnode];

  cl = new int[maxnode];
  halo = new int[maxnode];

  init();
  finddist();
  cluster();
}
