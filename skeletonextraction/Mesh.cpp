#include "Mesh.h"

void Mesh::loadOff(char *name)
{
	FILE *fPtr = fopen(name, "r");
	char str[334];

	fscanf(fPtr, "%s", str);

	int nVerts, nTris, n, i = 0;
	float x, y, z;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nTris, &n);
	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addVertex(x, y, z);
	}

	while (fscanf(fPtr, "%d", &i) != EOF)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addTriangle((int)x, (int)y, (int)z);
	}

	fclose(fPtr);
}

void Mesh::writeMesh(char *name)
{
	FILE *fPtr = fopen(name, "w");
	char str[334];

	fprintf(fPtr, "%s\n", str);

	int nVerts = verts.size(), nTris = tris.size(), n = 0, i = 0;
	float x, y, z;

	fprintf(fPtr, "%d %d %d\n", nVerts, nTris, n);
	while (i++ < nVerts)
	{
		x = verts[i]->coords[0];
		y = verts[i]->coords[1];
		z = verts[i]->coords[2];
		fprintf(fPtr, "%f %f %f\n", x, y, z);
	}
	int j = 0;
	while (j++ < nTris)
	{
		x = tris[j]->v1i;
		y = tris[j]->v2i;
		z = tris[j]->v3i;
		fprintf(fPtr, "%f %f %f\n", (int)x, (int)y, (int)z);
	}

	fclose(fPtr);
}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back(new Triangle(idx, v1, v2, v3));
	tris[idx]->v1p = verts[v1];
	tris[idx]->v2p = verts[v2];
	tris[idx]->v3p = verts[v3];

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	verts[v1]->triPList.push_back(tris[idx]);
	verts[v2]->triPList.push_back(tris[idx]);
	verts[v3]->triPList.push_back(tris[idx]);

	if (!makeVertsNeighbor(v1, v2))
		addEdge(v1, v2);

	if (!makeVertsNeighbor(v1, v3))
		addEdge(v1, v3);

	if (!makeVertsNeighbor(v2, v3))
		addEdge(v2, v3);
}

bool Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->vertList.size(); i++)
		if (verts[v1i]->vertList[i] == v2i)
			return true;

	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);

	verts[v1i]->vertPList.push_back(verts[v2i]);
	verts[v2i]->vertPList.push_back(verts[v1i]);
	return false;
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float *c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back(new Vertex(idx, c));
}

void Mesh::addEdge(int v1, int v2)
{
	int idx = edges.size();

	edges.push_back(new Edge(idx, v1, v2));

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);

	verts[v1]->edgePList.push_back(edges[idx]);
	verts[v2]->edgePList.push_back(edges[idx]);
}

void Mesh::do_edge_collapse_process()
{
	std::clock_t start;
	double duration;

	start = std::clock();
	compute_all_edge_collapse_costs();

	while (tris.size() > 0 && verts.size() > 0)
	{
		Vertex *mn = minimum_cost_edge();
		if (mn != mn->collapse)
		{
			collapse(mn, mn->collapse);
			update();
			compute_all_edge_collapse_costs();
			check_faces_edges();
		}
	}

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "edge collapses are completed in " << duration << " sec" << '\n';
	dijkstra_array_sol();
}
//-----------------------------------------------------------------------------
float Mesh::get_distance(float *v1, float *v2)
{
	MatrixXf Ki(3, 4);
	VectorXf edgev(3);
	edgev(0) = v2[0] - v1[0];
	edgev(1) = v2[1] - v1[1];
	edgev(2) = v2[2] - v1[2];

	float normp = std::sqrt(pow(edgev(0), 2) + pow(edgev(1), 2) + pow(edgev(2), 2));

	return normp;
}
//-----------------------------------------------------------------------------
float Mesh::get_distance(point_t v1, point_t v2)
{
	MatrixXf Ki(3, 4);
	VectorXf edgev(3);
	edgev(0) = v2[0] - v1[0];
	edgev(1) = v2[1] - v1[1];
	edgev(2) = v2[2] - v1[2];

	float normp = std::sqrt(pow(edgev(0), 2) + pow(edgev(1), 2) + pow(edgev(2), 2));

	return normp;
}
//-----------------------------------------------------------------------------
float Mesh::get_edge_cost(float *v1, float *v2, float *vi, float *pj)
{
	MatrixXf Ki(3, 4);
	VectorXf edgev(3);
	edgev(0) = v2[0] - v1[0];
	edgev(1) = v2[1] - v1[1];
	edgev(2) = v2[2] - v1[2];

	float normp = std::sqrt(pow(edgev(0), 2) + pow(edgev(1), 2) + pow(edgev(2), 2));

	//return normp;

	float *ai = new float[3];
	ai[0] = edgev(0) / normp;
	ai[1] = edgev(1) / normp;
	ai[2] = edgev(2) / normp;

	float *bi = new float[3];
	bi[0] = ai[1] * vi[2] - ai[2] * vi[1];
	bi[1] = -ai[0] * vi[2] + ai[2] * vi[0];
	bi[2] = ai[0] * vi[1] - ai[1] * vi[0];

	Ki(0, 0) = 0;
	Ki(0, 1) = -ai[2];
	Ki(0, 2) = ai[1];
	Ki(0, 3) = -bi[0];

	Ki(1, 0) = ai[2];
	Ki(1, 1) = 0;
	Ki(1, 2) = -ai[0];
	Ki(1, 3) = -bi[1];

	Ki(2, 0) = -ai[1];
	Ki(2, 1) = ai[0];
	Ki(2, 2) = 0;
	Ki(2, 3) = -bi[2];
	MatrixXf KiT(4, 3);
	KiT = Ki.transpose();
	MatrixXf qi = KiT * Ki;
	Vector3f pjm(pj[0], pj[1], pj[2]);
	RowVector4f pjmT = pjm.homogeneous().transpose();
	MatrixXf f1 = pjmT * qi;
	MatrixXf fi = f1 * pjm.homogeneous();
	return fi(0, 0);
}
//-----------------------------------------------------------------------------
void Mesh::compute_triangle_normal(Triangle *tri)
{
	Vector3f v0(verts[tri->v1i]->coords[0], verts[tri->v1i]->coords[1], verts[tri->v1i]->coords[2]);
	Vector3f v1(verts[tri->v2i]->coords[0], verts[tri->v2i]->coords[1], verts[tri->v2i]->coords[2]);
	Vector3f v2(verts[tri->v3i]->coords[0], verts[tri->v3i]->coords[1], verts[tri->v3i]->coords[2]);

	Vector3f normal = (v1 - v0).cross(v2 - v1);

	float mnorm = std::sqrt(pow(normal(0), 2) + pow(normal(1), 2) + pow(normal(2), 2));
	if (mnorm == 0)
		return;
	tri->normal = new float[3];
	tri->normal[0] = normal(0) / mnorm;
	tri->normal[1] = normal(1) / mnorm;
	tri->normal[2] = normal(2) / mnorm;
}
//-----------------------------------------------------------------------------

void Mesh::replace_vertex(Triangle *tri, Vertex *vold, Vertex *vnew)
{
	Vertex *v0 = tri->v1p;
	Vertex *v1 = tri->v2p;
	Vertex *v2 = tri->v3p;

	assert(vold && vnew);
	if (vold != v0 && vold != v1 && vold == v2)
		return;
	if (vnew == v0 || vnew == v1 || vnew == v2)
		return;
	if (vold == v0)
	{
		v0 = vnew;
		update();
	}
	else if (vold == v1)
	{
		v1 = vnew;
		update();
	}
	else
	{
		if (vold != v2)
			return;
		v2 = vnew;
		update();
	}

	int i = 0;

	Remove(tri->idx, vold->triList);
	Remove(tri, vold->triPList);
	update();

	if (!Contains(tri->idx, vnew->triList))
		return;

	vnew->triList.push_back(tri->idx);
	vnew->triPList.push_back(tri);
	update();

	for (i = 0; i < 3; i++)
	{
		Vertex *vi = i == 0 ? v0 : (i == 1 ? v1 : v2);
		remove_non_adjacent(vold, vi);
		update();
		remove_non_adjacent(vi, vold);
		update();
	}
	for (i = 0; i < 3; i++)
	{
		Vertex *vi = i == 0 ? v0 : (i == 1 ? v1 : v2);
		if (!Contains(tri->idx, vi->triList) == 1)
			return;
		for (int j = 0; j < 3; j++)
			if (i != j)
			{
				Vertex *vj = j == 0 ? v0 : (j == 1 ? v1 : v2);
				vi->vertList.push_back(vj->idx);
				vi->vertPList.push_back(vj);
				update();
			}
	}
	if (tri->v1p == vold)
	{
		tri->v1p = vnew;
		tri->v1i = vnew->idx;
	}
	else if (tri->v2p == vold)
	{
		tri->v2p = vnew;
		tri->v2i = vnew->idx;
	}
	else if (tri->v3p == vold)
	{
		tri->v3p = vnew;
		tri->v3i = vnew->idx;
	}
	compute_triangle_normal(tri);
}
//-----------------------------------------------------------------------------
void Mesh::remove_non_adjacent(Vertex *thisVert, Vertex *n)
{
	// removes n from neighbor list if n isn't a neighbor.
	if (!Contains(n->idx, thisVert->vertList))
		return;
	for (int i = 0; i < thisVert->triList.size(); i++)
	{
		Triangle *tri = tris[thisVert->triList[i]];
		if (tri->HasVertex(n))
			return;
	}
	Remove(n->idx, thisVert->vertList);
	Remove(n, thisVert->vertPList);
	update();
}
//-----------------------------------------------------------------------------
void Mesh::compute_edge_cost_at_vertex(Vertex *v)
{
	if (v->vertList.size() == 0)
	{
		v->collapse = NULL;
		return;
	}
	v->collapse = NULL;
	float fmincost = INT_MAX;

	for (int i = 0; i < v->vertPList.size(); i++)
	{
		if (v->vertPList[i] == v)
			continue;
		float fcost = compute_edge_collapse_cost(v, v->vertPList[i]);
		if (fcost < fmincost && v->collapse != v)
		{
			if (IsAFace(v, v->vertPList[i]))
			{
				v->collapse = v->vertPList[i]; // candidate for edge collapse
				v->cost = fcost;
				fmincost = fcost;
			}
		}
	}
}
//-----------------------------------------------------------------------------
void Mesh::compute_all_edge_collapse_costs()
{
	for (int i = 0; i < verts.size(); i++)
	{
		compute_edge_cost_at_vertex(verts[i]);
	}
}
//-----------------------------------------------------------------------------
float Mesh::compute_edge_collapse_cost(Vertex *u, Vertex *v)
{
	float fcost1 = fabs(get_edge_cost(u->coords, v->coords, u->coords, v->coords));
	return fcost1;
}
//-----------------------------------------------------------------------------
void Mesh::collapse(Vertex *u, Vertex *v)
{
	update();
	if (!v)
	{
		delete_vertex(u);
		return;
	}
	int i = 0;
	vector<Vertex *> tmp;

	for (i = 0; i < u->vertList.size(); i++)
	{
		tmp.push_back(verts[u->vertList[i]]);
	}

	for (i = u->triPList.size() - 1; i >= 0; i--)
	{
		if (u->triPList[i]->HasVertex(v))
		{
			delete_face(u->triPList[i]);
		}
	}

	for (i = u->edgePList.size() - 1; i >= 0; i--)
	{
		if (u->edgePList[i]->v1p == u)
		{
			u->edgePList[i]->v1p = v;
			u->edgePList[i]->v1i = v->idx;
		}
		else if (u->edgePList[i]->v2p == u)
		{
			u->edgePList[i]->v2p = v;
			u->edgePList[i]->v2i = v->idx;
		}
	}

	update_verts();

	for (i = u->triPList.size() - 1; i >= 0; i--)
	{
		replace_vertex(u->triPList[i], u, v);
	}

	delete_vertex(u);

	for (i = 0; i < tmp.size(); i++)
	{
		compute_edge_cost_at_vertex(tmp[i]);
	}
	v->cost = v->cost + u->cost;
}

Vertex *Mesh::minimum_cost_edge()
{
	Vertex *mn = verts[0];
	for (int i = 0; i < verts.size(); i++)
	{
		if (verts[i]->cost < mn->cost)
		{
			mn = verts[i];
		}
	}
	return mn;
}
//-----------------------------------------------------------------------------
void Mesh::delete_vertex(Vertex *vd)
{
	int vidx = vd->idx;
	verts.erase(verts.begin() + vidx);

	update();
}
//-----------------------------------------------------------------------------
void Mesh::update_verts()
{
	for (int vi = 0; vi < verts.size(); vi++)
		verts[vi]->idx = vi;

	for (int vi = 0; vi < verts.size(); vi++)
	{
		for (int i = 0; i < verts[vi]->edgePList.size(); i++)
		{
			if (verts[vi]->edgePList[i]->v1p && verts[vi]->edgePList[i]->v2p)
			{
				verts[vi]->edgePList[i]->v1i = verts[vi]->edgePList[i]->v1p->idx;
				verts[vi]->edgePList[i]->v2i = verts[vi]->edgePList[i]->v2p->idx;
			}
			else
			{
				verts[vi]->edgePList.erase(verts[vi]->edgePList.begin() + i);
			}
		}

		for (int j = 0; j < verts[vi]->triPList.size(); j++)
		{
			verts[vi]->triPList[j]->v1i = verts[vi]->triPList[j]->v1p->idx;
			verts[vi]->triPList[j]->v2i = verts[vi]->triPList[j]->v2p->idx;
			verts[vi]->triPList[j]->v3i = verts[vi]->triPList[j]->v3p->idx;
		}

		for (int k = 0; k < verts[vi]->vertPList.size(); k++)
		{
			verts[vi]->vertList[k] = verts[vi]->vertPList[k]->idx;
		}
	}
}
//-----------------------------------------------------------------------------
void Mesh::delete_face(Triangle *td)
{
	int fidx = td->idx;
	if (td->v1p)
	{
		if (Contains(td, td->v1p->triPList))
			Remove(td, td->v1p->triPList);
		if (Contains(td->idx, td->v1p->triList))
			Remove(td->idx, td->v1p->triList);
	}

	if (td->v2p)
	{
		if (Contains(td, td->v2p->triPList))
			Remove(td, td->v2p->triPList);
		if (Contains(td->idx, td->v2p->triList))
			Remove(td->idx, td->v2p->triList);
	}

	if (td->v3p)
	{
		if (Contains(td, td->v3p->triPList))
			Remove(td, td->v3p->triPList);
		if (Contains(td->idx, td->v3p->triList))
			Remove(td->idx, td->v3p->triList);
	}

	//Remove(td, tris);
	if (tris.size() > fidx)
		tris.erase(tris.begin() + fidx);

	update_tris();
}
//-----------------------------------------------------------------------------
void Mesh::check_faces_edges()
{
	for (int i = 0; i < tris.size(); ++i)
	{
		if (tris[i]->v1p == NULL || tris[i]->v2p == NULL || tris[i]->v3p == NULL ||
			tris[i]->v1i > verts.size() || tris[i]->v2i > verts.size() || tris[i]->v3i > verts.size())
		{
			delete_face(tris[i]);
			--i;
		}
	}
	for (int i = 0; i < verts.size(); ++i)
	{
		for (int j = 0; j < verts[i]->vertPList.size(); ++j)
		{
			if (verts[i]->vertPList[j]->idx > verts.size())
			{
				Remove(verts[i]->vertPList[j], verts[i]->vertPList);
			}
		}

		for (int j = 0; j < verts[i]->vertList.size(); ++j)
		{
			if (verts[i]->vertList[j] > verts.size())
			{
				Remove(verts[i]->vertList[j], verts[i]->vertList);
			}
		}
			
	} 

	update();
}
//-----------------------------------------------------------------------------
void Mesh::delete_edge(Edge *ed)
{
	int eidx = ed->idx;
	if (ed->v1p)
	{
		if (Contains(ed, ed->v1p->edgePList))
			Remove(ed, ed->v1p->edgePList);
		if (Contains(ed->idx, ed->v1p->edgeList))
			Remove(ed->idx, ed->v1p->edgeList);
	}

	if (ed->v2p)
	{
		if (Contains(ed, ed->v2p->edgePList))
			Remove(ed, ed->v2p->edgePList);
		if (Contains(ed->idx, ed->v2p->edgeList))
			Remove(ed->idx, ed->v2p->edgeList);
	}

	//Remove(ed, edges);
	if (edges.size() > eidx)
		edges.erase(edges.begin() + eidx);

	update_edges();
}
//-----------------------------------------------------------------------------
void Mesh::update_tris()
{
	for (vector<Triangle *>::iterator it = tris.begin(); it != tris.end(); it++)
		(*it)->idx = it - tris.begin();

	for (int vi = 0; vi < verts.size(); vi++)
	{
		for (int j = 0; j < verts[vi]->triPList.size(); j++)
		{
			verts[vi]->triList[j] = verts[vi]->triPList[j]->idx;
		}
	}
}
//-----------------------------------------------------------------------------
void Mesh::update_edges()
{
	for (vector<Edge *>::iterator it = edges.begin(); it != edges.end(); it++)
		(*it)->idx = it - edges.begin();

	for (int vi = 0; vi < verts.size(); vi++)
	{
		for (int j = 0; j < verts[vi]->edgePList.size(); j++)
		{
			verts[vi]->edgeList[j] = verts[vi]->edgePList[j]->idx;
		}
	}
}
//-----------------------------------------------------------------------------
void Mesh::update()
{
	update_verts();
	update_tris();
	update_edges();
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void Mesh::farthest_sampling_by_sphere(vector<vector<Vertex *>> &solution_set)
{
	int n = verts.size();

	vector<float> soldist;
	vector<vector<Vertex *>> solution_set2;
	vector<vector<int>> solIds2;
	
	for (int k = 0; k < n; k++)
	{
		vector<Vertex *> spls;
		vector<int> sids;

		vector<int> nIdxs;
		vector<float> nDsts;
		for (int pi = 0; pi < verts.size(); pi++)
		{
			if (pi != k)
			{
				nDsts.push_back(get_distance(verts[k]->coords, verts[pi]->coords));
				nIdxs.push_back(verts[pi]->idx);
			}
		}

		int maxIdx = 0;
		float maxVal = 0.;
		for (int mi = 0; mi < nDsts.size(); mi++)
		{
			if (nDsts[mi] > maxVal)
			{
				maxIdx = mi;
				maxVal = nDsts[mi];
			}
		}
		spls.push_back(verts[maxIdx]);
		spls.push_back(verts[k]);
		sids.push_back(k);
		sids.push_back(maxIdx);
		solution_set2.push_back(spls);
		solIds2.push_back(sids);
		soldist.push_back(maxVal);
	}
	int pntN = 0;
	while (pntN < 3)
	{
		float maxDist = 0.;
		int maxIndx = 0;
		for (int si = 0; si < soldist.size(); si++)
		{
			if (soldist[si] > maxDist)
			{
				maxDist = soldist[si];
				maxIndx = si;
			}
		}
		solution_set.push_back(solution_set2[maxIndx]);
		solution_set2.erase(solution_set2.begin() + maxIndx);
		++pntN;
	}
	


}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void Mesh::dijkstra_array_sol()
{
	vector<vector<Vertex *>> solution_set;
	farthest_sampling_by_sphere(solution_set);
	this->edges.clear();

	int vertsize = verts.size();

	MatrixXi graph(vertsize, vertsize);
	graph.setZero();

	for (int ssi = 0; ssi < solution_set.size(); ssi++)
	{
		for (int i = 0; i < vertsize; ++i)
		{
			for (int j = 0; j < verts[i]->vertPList.size(); ++j)
			{
				if (verts[i]->vertPList[j]->idx < vertsize)
					graph(i,verts[i]->vertPList[j]->idx) = (int)get_distance(verts[i]->coords, verts[i]->vertPList[j]->coords);
			}	
		}
		vector<vector<int>> paths;
		dijkstra(graph, vertsize, vertsize, solution_set[ssi][0]->idx, solution_set[ssi][1]->idx, paths);
		if (paths.size() > solution_set[ssi][1]->idx)
		{
			for (int k = 0; k < paths[solution_set[ssi][1]->idx].size() - 1; k++)
			{
				addEdge(paths[solution_set[ssi][1]->idx][k], paths[solution_set[ssi][1]->idx][k + 1]);
			}
			cout << ssi << endl;
		}
		else 
		{
			for (int i = 0; i < paths.size(); i++)
			{
				for (int k = 0; k < paths[i].size() - 1; k++)
				{
					addEdge(paths[i][k], paths[i][k + 1]);
				}
				cout << ssi << endl;
			}
		}
	}
}
//-----------------------------------------------------------------------------
void dijkstra(MatrixXi G, int maxn, int n, int startnode, int endnode, vector<vector<int>>& paths)
{
	int **cost = NULL;
	int *distance = new int[maxn];
	int *pred = new int[maxn];
	int *visited = new int[maxn];
	cost = new int*[maxn];
	for (int ci = 0; ci < maxn; ci++)
		cost[ci] = new int[maxn];

	int count = 0, mindistance, nextnode = 0, i = 0, j = 0;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			if (G(i,j) == 0)
				cost[i][j] = INFINITY;
			else
				cost[i][j] = G(i,j);
	for (i = 0; i < n; i++)
	{
		distance[i] = cost[startnode][i];
		pred[i] = startnode;
		visited[i] = 0;
	}
	distance[startnode] = 0;
	visited[startnode] = 1;
	count = 1;
	while (count < n - 1)
	{
		mindistance = INFINITY;
		for (i = 0; i < n; i++)
			if (distance[i] < mindistance && !visited[i])
			{
				mindistance = distance[i];
				nextnode = i;
			}
		visited[nextnode] = 1;
		for (i = 0; i < n; i++)
			if (!visited[i])
				if (mindistance + cost[nextnode][i] < distance[i])
				{
					distance[i] = mindistance + cost[nextnode][i];
					pred[i] = nextnode;
				}
		count++;
	}

	for (i = 0; i < n; i++)
	{
		if (i != startnode)
		{
			vector<int> nPath ;
			nPath.resize(500); nPath.clear();
			nPath.push_back(i);
			j = i;
			do
			{
				j = pred[j];
				nPath.push_back(j);
			} while (j != startnode);

			paths.push_back(nPath);
		}
	}
}
//-----------------------------------------------------------------------------
float* Mesh::find_bounding_sphere()
{
    float sphereRad = 0;
	register double dx, dy, dz;
    register double rad_sq, xspan, yspan, zspan, maxspan;
    double old_to_p, old_to_p_sq, old_to_new;
    float* xmin = new float[3], *xmax = new float[3], *ymin = new float[3], *ymax = new float[3], 
		   *zmin = new float[3], *zmax = new float[3], *dia1 = new float[3], *dia2 = new float[3];

    /* FIRST PASS: find 6 minima/maxima points */
    xmin[0] = ymin[1] = zmin[2] = INT_MAX; /* initialize for min/max compare */
    xmax[0] = ymax[1] = zmax[2] = -INT_MAX;

    for (int viter = 0; viter < verts.size(); ++viter)
    {
        float* caller_p = verts[viter]->coords;

        if (caller_p[0] < xmin[0])
            xmin = caller_p; /* New xminimum point */
        if (caller_p[0] > xmax[0])
            xmax = caller_p;
        if (caller_p[1] < ymin[1])
            ymin = caller_p;
        if (caller_p[1] > ymax[1])
            ymax = caller_p;
        if (caller_p[2] < zmin[2])
            zmin = caller_p;
        if (caller_p[2] > zmax[2])
            zmax = caller_p;
    }
    /* Set xspan = distance between the 2 points xmin & xmax (squared) */
    dx = xmax[0] - xmin[0];
    dy = xmax[1] - xmin[1];
    dz = xmax[2] - xmin[2];
    xspan = dx * dx + dy * dy + dz * dz;

    /* Same for y & z spans */
    dx = ymax[0] - ymin[0];
    dy = ymax[1] - ymin[1];
    dz = ymax[2] - ymin[2];
    yspan = dx * dx + dy * dy + dz * dz;

    dx = zmax[0] - zmin[0];
    dy = zmax[1] - zmin[1];
    dz = zmax[2] - zmin[2];
    zspan = dx * dx + dy * dy + dz * dz;

    /* Set points dia1 & dia2 to the maximally separated pair */
    dia1 = xmin;
    dia2 = xmax; /* assume xspan biggest */
    maxspan = xspan;
    if (yspan > maxspan)
    {
        maxspan = yspan;
        dia1 = ymin;
        dia2 = ymax;
    }
    if (zspan > maxspan)
    {
        dia1 = zmin;
        dia2 = zmax;
    }
	float* sphereCenterPoint = new float[3];
    sphereCenterPoint[0] = (dia1[0] + dia2[0]) / 2.0;
    sphereCenterPoint[1] = (dia1[1] + dia2[1]) / 2.0;
    sphereCenterPoint[2] = (dia1[2] + dia2[2]) / 2.0;
    /* calculate initial radius**2 and radius */
    dx = dia2[0] - sphereCenterPoint[0]; /* x component of radius vector */
    dy = dia2[1] - sphereCenterPoint[1]; /* y component of radius vector */
    dz = dia2[2] - sphereCenterPoint[2]; /* z component of radius vector */
    rad_sq = dx * dx + dy * dy + dz * dz;
    sphereRad = sqrt(rad_sq);

    /* SECOND PASS: increment current sphere */

    for (int viter = 0; viter < verts.size(); ++viter)
    {
        float* caller_p = verts[viter]->coords;
        dx = caller_p[0] - sphereCenterPoint[0];
        dy = caller_p[1] - sphereCenterPoint[1];
        dz = caller_p[2] - sphereCenterPoint[2];
        old_to_p_sq = dx * dx + dy * dy + dz * dz;
        if (old_to_p_sq > rad_sq) /* do r**2 test first */
        {                         /* this point is outside of current sphere */
            old_to_p = sqrt(old_to_p_sq);
            /* calc radius of new sphere */
            sphereRad = old_to_p;           //(sphereRad + old_to_p) / 2.0;
            rad_sq = sphereRad * sphereRad; /* for next r**2 compare */
            old_to_new = old_to_p - sphereRad;
            /* calc center of new sphere */
            sphereCenterPoint[0] = (sphereRad * sphereCenterPoint[0] + old_to_new * caller_p[0]) / old_to_p;
            sphereCenterPoint[1] = (sphereRad * sphereCenterPoint[1] + old_to_new * caller_p[1]) / old_to_p;
            sphereCenterPoint[2] = (sphereRad * sphereCenterPoint[2] + old_to_new * caller_p[2]) / old_to_p;
        }
    }
	//addVertex(sphereCenterPoint[0], sphereCenterPoint[1], sphereCenterPoint[2]);

	return sphereCenterPoint;
}
//-----------------------------------------------------------------------------
