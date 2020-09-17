#pragma once

#include <iostream>
#include <vector>
#include "include\Eigen\Sparse"
#include "include\Eigen\Dense"
#include "include\Eigen\Geometry"
#include <ctime>
#include "KDTree.hpp"

using namespace std;
using namespace Eigen;

struct Edge;
struct Triangle;

void dijkstra(MatrixXi G, int maxn, int n, int startnode, int endnode, vector<vector<int>>& paths);
struct Vertex
{
	float *coords, *normals; //3d coordinates etc
	int idx;				 //who am i; verts[idx]

	vector<int> vertList;		//adj vvertices;
	vector<Vertex *> vertPList; //adj vvertices;
	vector<int> triList;
	vector<Triangle *> triPList;
	vector<int> edgeList;
	vector<Edge *> edgePList;

	vector<Vertex *> collCandList;

	Vertex *collapse;
	float cost;

	Vertex(int i, float *c) : idx(i), coords(c), collapse(NULL), cost(0.){};
};

struct Edge
{
	int idx;	  //edges[idx]
	int v1i, v2i; //endpnts
	Vertex *v1p, *v2p;
	float length;
	Edge(int id, int v1, int v2) : idx(id), v1i(v1), v2i(v2), v1p(NULL), v2p(NULL) { computeLength(); };

	void computeLength()
	{
		length = 7;
	}
	int HasVertex(Vertex *v)
	{
		return (v->idx == v1i || v->idx == v2i || v == v1p || v == v2p);
	}
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
	Vertex *v1p, *v2p, *v3p;
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3), v1p(NULL), v2p(NULL), v3p(NULL){};
	float *normal;
	int HasVertex(Vertex *v)
	{
		return (v->idx == v1i || v->idx == v2i || v->idx == v3i || v == v1p || v == v2p || v == v3p);
	}
};

class Mesh
{
private:
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	bool makeVertsNeighbor(int v1i, int v2i);

	void compute_triangle_normal(Triangle *tri);
	void replace_vertex(Triangle *tri, Vertex *vold, Vertex *vnew);
	void compute_edge_cost_at_vertex(Vertex *v);
	void compute_all_edge_collapse_costs();
	float compute_edge_collapse_cost(Vertex *u, Vertex *v);

	float get_edge_cost(float *v1, float *v2, float *vi, float *pj);
	void remove_non_adjacent(Vertex *thisVert, Vertex *n);
	void collapse(Vertex *u, Vertex *v);
	Vertex *minimum_cost_edge();

	void delete_vertex(Vertex *vd);
	void delete_face(Triangle *td);
	void delete_edge(Edge *ed);
	void update_verts();
	void update_tris();
	void update_edges();
	void update();
	void check_faces_edges();
	float get_distance(float *v1, float *v2);
	float get_distance(point_t v1, point_t v2);
	void farthest_sampling_by_sphere(vector<vector<Vertex *>> &solution_set);
	void dijkstra_array_sol();
	float* find_bounding_sphere();
public:
	vector<Vertex *> verts;
	vector<Triangle *> tris;
	vector<Edge *> edges;

	vector<Edge *> skeletonedges;

	Mesh(){};

	void loadOff(char *name);
	void do_edge_collapse_process();
	void writeMesh(char *name);

	template <class Type>
	void Remove(Type t, vector<Type> vlist)
	{
		int i = 0;
		for (i = 0; i < vlist.size(); i++)
		{
			if (vlist[i] == t)
			{
				break;
			}
		}

		vlist.erase(vlist.begin() + i);
		for (i = 0; i < vlist.size(); i++)
		{
			if (vlist[i] == t)
				vlist.erase(vlist.begin() + i);
		}
	}

	template <class Type>
	int Contains(Type t, vector<Type> vlist)
	{
		int i = 0;
		int count = 0;
		for (i = 0; i < vlist.size(); i++)
		{
			if (vlist[i] == t)
				count++;
		}
		return count;
	}

	bool IsAFace(Vertex *t1, Vertex *t2)
	{
		if (Contains(t2, t1->vertPList) && Contains(t1, t2->vertPList))
		{
			vector<Vertex *> kList;

			for (int i = 0; i < t1->vertPList.size(); i++)
			{
				for (int j = 0; j < t2->vertPList.size(); j++)
				{
					if (t1->vertPList[i] == t2->vertPList[j])
						kList.push_back(t1->vertPList[i]);
				}
			}
			vector<Triangle *> ktriList;
			for (int m = 0; m < t1->triPList.size(); m++)
			{
				if (t1->triPList[m]->v1p == t2 || t1->triPList[m]->v2p == t2 || t1->triPList[m]->v3p == t2)
					ktriList.push_back(t1->triPList[m]);
			}

			for (int k = 0; k < kList.size(); k++)
			{
				for (int m = 0; m < ktriList.size(); m++)
				{
					if (ktriList[m]->v1p == kList[k] || ktriList[m]->v2p == kList[k] || ktriList[m]->v3p == kList[k])
						return true;
				}
			}
			return false;
		}
		else
			return false;
	}
};
