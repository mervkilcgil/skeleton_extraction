#include "SkelatonExtractor.h"

SkelatonExtractor::SkelatonExtractor()
    : MeshViewer()
{
    init();
}

SkelatonExtractor::~SkelatonExtractor()
{
}

void SkelatonExtractor::init()
{
    std::string fname = "times_";
    fname = fname + ".txt";
    outfile1.open(fileName.c_str(), std::ofstream::out);
}

bool SkelatonExtractor::open_mesh(const char *_meshfilename)
{
    // load mesh
    if (MeshViewer::open_mesh(_meshfilename))
    {
        mesh_.add_property(vpos_);
        mesh_.add_property(ecotanweight_);
        mesh_.add_property(farea_);
        mesh_.add_property(v1ringarea_);

        // compute face & vertex normals
        mesh_.update_normals();
        // update face indices for faster rendering
        update_face_indices();
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void SkelatonExtractor::draw()
{
    MatrixXf vpos(mesh_.n_vertices(), 3);
    make_laplacian_smooth(5, vpos);
    //edge_collapses(vpos);
    outfile1.close();
    mesh_.update_normals();
    //update_face_indices();
    nfname = "new_" + fileName;
    OpenMesh::IO::write_mesh(mesh_, nfname.c_str());
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void SkelatonExtractor::get_laplacian_weight_matrix(std::map<int, float> oneringAreas,
                                                    std::map<int, float> edgeWeights,
                                                    std::vector<Triplet<float>> &L_coeff,
                                                    MatrixXf &WH0_diag, MatrixXf &vpos,
                                                    float &originalFaceAreaSum)
{
    MatrixXf vertexWeights(mesh_.n_vertices(), mesh_.n_vertices());
    vertexWeights.setZero();

    L_coeff.clear();

    for (Mesh::FaceIter fit = mesh_.faces_begin(); fit != mesh_.faces_end(); fit++)
    {
        Mesh::FaceVertexCWIter fvit = mesh_.fv_cwbegin(*fit);
        Mesh::VertexHandle v1 = *fvit;
        ++fvit;
        Mesh::VertexHandle v2 = *fvit;
        ++fvit;
        Mesh::VertexHandle v3 = *fvit;
        ++fvit;

        int edgeId1 = mesh_.edge_handle(mesh_.find_halfedge(v1, v2)).idx();
        int edgeId2 = mesh_.edge_handle(mesh_.find_halfedge(v2, v3)).idx();
        int edgeId3 = mesh_.edge_handle(mesh_.find_halfedge(v3, v1)).idx();

        float edgeWeight1 = 0.;
        float edgeWeight2 = 0.;
        float edgeWeight3 = 0.;

        std::map<int, float>::iterator it = edgeWeights.find(edgeId1);
        if (it != edgeWeights.end())
            edgeWeight1 = (*it).second;

        it = edgeWeights.find(edgeId2);
        if (it != edgeWeights.end())
            edgeWeight2 = (*it).second;

        it = edgeWeights.find(edgeId3);
        if (it != edgeWeights.end())
            edgeWeight3 = (*it).second;

        vertexWeights(v1.idx(), v1.idx()) = vertexWeights(v1.idx(), v1.idx()) + edgeWeight1 + edgeWeight3;

        vertexWeights(v2.idx(), v2.idx()) = vertexWeights(v2.idx(), v2.idx()) + edgeWeight1 + edgeWeight2;

        vertexWeights(v3.idx(), v3.idx()) = vertexWeights(v3.idx(), v3.idx()) + edgeWeight2 + edgeWeight3;
    }

    for (Mesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter)
    {
        int vid = (*viter).idx();
        for (Mesh::VertexVertexIter witer = mesh_.vv_begin(*viter); witer != mesh_.vv_end(*viter); ++witer)
        {
            int wid = (*witer).idx();
            WH0_diag(vid, wid) = 1;
            int edgeId = mesh_.edge_handle(mesh_.find_halfedge((*viter), (*witer))).idx();

            std::map<int, float>::iterator it = edgeWeights.find(edgeId);
            if (it == edgeWeights.end())
                continue;
            double edgeWeight = (*it).second;
            
            L_coeff.push_back(Triplet<float>(vid, wid, edgeWeight));
        }

        if (isinf(vertexWeights(vid, vid)) || isnan(vertexWeights(vid, vid)))
        {
            if (isinf(vertexWeights(vid, vid)))
                cout << "infinite edge weight" << endl;
            if (isnan(vertexWeights(vid, vid)))
                cout << "nan edge weight" << endl;
            L_coeff.push_back(Triplet<float>(vid, vid, 0));
        }
        else
            L_coeff.push_back(Triplet<float>(vid, vid, -vertexWeights(vid, vid)));

        vpos(vid, 0) = mesh_.point(*viter)[0];
        vpos(vid, 1) = mesh_.point(*viter)[1];
        vpos(vid, 2) = mesh_.point(*viter)[2];
        for (Mesh::VertexIter viter2 = mesh_.vertices_begin(); viter2 != mesh_.vertices_end(); ++viter2)
        {
            int vid2 = (*viter2).idx();
            WH0_diag(vid, vid2) = 1;
        }

        std::map<int, float>::iterator nit = oneringAreas.find(vid);
        if (nit == oneringAreas.end())
            continue;
        originalFaceAreaSum += (*nit).second;
    }
}
//-----------------------------------------------------------------------------
void SkelatonExtractor::make_laplacian_smooth(int numit, MatrixXf &vpos)
{
    float SL = 2.;
    std::map<int, float> edgeWeights = compute_weights();

    std::map<int, float> oneringAreas = calculate_one_ring_area();

    std::vector<Triplet<float>> WL_diag_coeff;
    std::vector<Triplet<float>> L_coeff;

    int n = mesh_.n_vertices();

    MatrixXf zeros(n, 3);
    zeros.setZero();

    MatrixXf WH0_diag(n, n);
    MatrixXf WL_diag(n, n);

    float originalFaceAreaSum = 0.;
    float initialFaceWeight = 1 / 1000 * sqrt(get_average_face_area());

    get_laplacian_weight_matrix(oneringAreas, edgeWeights, L_coeff, WH0_diag, vpos, originalFaceAreaSum);

    SparseMatrix<float> Lsparse(n, n);
    Lsparse.reserve(n);
    Lsparse.setZero();
    Lsparse.setFromTriplets(L_coeff.begin(), L_coeff.end());

    VectorXi diag(1);
    diag(0) = 0;

    SparseMatrix<float> WH0 = spdiags(WH0_diag, diag, n, n);

    SparseMatrix<float> WH = spdiags(WH0_diag, diag, n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j){
            WL_diag(i, j) = initialFaceWeight;
            WH0_diag(i,j) = 1;
        }
    }

    SparseMatrix<float> WL = spdiags(WL_diag, diag, n, n);

    for (int i = 0; i < numit; ++i)
    {
        std::clock_t start;
        double duration;

        start = std::clock();
        MatrixXf WL_dense = WL.toDense();
        MatrixXf WH_dense = WH.toDense();
        MatrixXf lwl = Lsparse * WL;

        MatrixXf A(lwl.rows() + WH_dense.rows(), lwl.cols()); // <-- D(A.rows() + B.rows(), ...)
        A << lwl, WH_dense;

        MatrixXf hwl = WH_dense * vpos;

        MatrixXf b2(zeros.rows() + hwl.rows(), zeros.cols());
        b2 << zeros, hwl;

        MatrixXf cpts(n, 3);
        cpts.setZero();

        for (int j = 0; j < 3; j++)
        {
            VectorXf soln = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b2.col(j));
            for (int k = 0; k < n; k++)
                cpts(k, j) = soln(k);
        }

        for (Mesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter)
        {
            int vid = (*viter).idx();
            mesh_.set_point(*viter, Mesh::Point(cpts(vid, 0), cpts(vid, 1), cpts(vid, 2)));
        }

        std::map<int, float> noneringAreas = calculate_one_ring_area();

        float newringareas = 0.;
        for (Mesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter)
        {
            std::map<int, float>::iterator it = noneringAreas.find((*viter).idx());
            if (it == noneringAreas.end())
                continue;
            newringareas += (*it).second;
        }

        float changeinarea = std::pow(newringareas, -0.5);
        //area_ratios.append(np.sum(newringareas) / originalFaceAreaSum)

        WL_dense = WL_dense.array() * SL;
        WH_dense = WH0.toDense().array() * changeinarea;

        std::map<int, float> nedgeWeights = compute_weights();

        vpos.setZero();

        float nFaceAreaSum = 0;
        MatrixXf nWH0_diag(n, n);

        L_coeff.clear();
        get_laplacian_weight_matrix(noneringAreas, nedgeWeights, L_coeff, nWH0_diag, vpos, nFaceAreaSum);

        Lsparse.setZero();
        Lsparse.setFromTriplets(L_coeff.begin(), L_coeff.end());

        mesh_.update_vertex_normals();
        mesh_.update_face_normals();
        duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        cout << i + 1 << " iteration is completed in " << duration << " sec" << '\n';
        outfile1 << std::to_string(duration) << std::endl;
    }
}
//-----------------------------------------------------------------------------
std::map<int, float> SkelatonExtractor::compute_weights()
{
    std::map<int, float> edgeWeights;
    Mesh::EdgeIter eiter, e_end(mesh_.edges_end());
    for (eiter = mesh_.edges_begin(); eiter != e_end; ++eiter)
    {
        float wt = 0.;

        Mesh::HalfedgeHandle he = mesh_.halfedge_handle(*eiter, 0);
        float sangle = mesh_.calc_sector_angle(he);
        if (sangle != 0 && !isinf(1 / tan(sangle)))
            wt = 1 / tan(sangle);
        Mesh::HalfedgeHandle sh = mesh_.halfedge_handle(*eiter, 1);
        float sangle2 = 0.;
        if (sh.is_valid())
        {
            sangle2 = mesh_.calc_sector_angle(sh);
            if (sangle2 != 0 && !isinf(1 / tan(sangle2)))
                wt += 1 / tan(sangle2);
        }
        //mesh_.property(ecotanweight_, *eiter) = wt;
        bool isEdgeValid = sangle != 0;
        if (sh.is_valid())
            isEdgeValid = isEdgeValid && sangle2 != 0;
        if (isEdgeValid)
        {
            std::pair<int, float> p1;
            p1.first = (*eiter).idx();
            p1.second = wt;
            edgeWeights.insert(p1);
        }
    }
    return edgeWeights;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
std::map<int, float> SkelatonExtractor::compute_face_area()
{
    std::map<int, float> faceAreas;
    for (Mesh::FaceIter fit = mesh_.faces_begin(); fit != mesh_.faces_end(); fit++)
    {
        Mesh::FaceVertexCWIter fvit = mesh_.fv_cwbegin(*fit);
        Mesh::Point v1 = mesh_.point(*fvit);
        ++fvit;
        Mesh::Point v2 = mesh_.point(*fvit);
        ++fvit;
        Mesh::Point v3 = mesh_.point(*fvit);
        ++fvit;

        float a = get_distance_bw_2points(v1, v2);
        float b = get_distance_bw_2points(v2, v3);
        float c = get_distance_bw_2points(v3, v1);

        float s = (a + b + c) / 2.;

        float farea = std::sqrt(s * (s - a) * (s - b) * (s - c));

        //mesh_.property(farea_, *fit) = farea;
        std::pair<int, float> p1;
        p1.first = (*fit).idx();
        p1.second = farea;
        faceAreas.insert(p1);
    }
    return faceAreas;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
float SkelatonExtractor::get_average_face_area()
{
    std::map<int, float> faceareas = compute_face_area();
    float fsum = 0.;
    int nfaces = 0;

    for (Mesh::FaceIter fit = mesh_.faces_begin(); fit != mesh_.faces_end(); fit++)
    {
        //fsum += mesh_.property(farea_, *fit);
        std::map<int, float>::iterator it = faceareas.find((*fit).idx());
        if (it == faceareas.end())
            continue;
        fsum += (*it).second;
        ++nfaces;
    }
    return fsum / nfaces;
}
//-----------------------------------------------------------------------------
float SkelatonExtractor::get_distance_bw_2points(Mesh::Point p1, Mesh::Point p2)
{
    return std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1], 2) + std::pow(p1[2] - p2[2], 2));
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
std::map<int, float> SkelatonExtractor::calculate_one_ring_area()
{
    std::map<int, float> oneRingAreas;
    std::map<int, float> faceareas = compute_face_area();
    for (Mesh::VertexIter vit = mesh_.vertices_begin(); vit != mesh_.vertices_end(); ++vit)
    {
        float oneringarea = 0;
        for (Mesh::VertexFaceCWIter vfit = mesh_.vf_cwbegin(*vit); vfit != mesh_.vf_cwend(*vit); ++vfit)
        {
            std::map<int, float>::iterator it = faceareas.find((*vfit).idx());
            if (it == faceareas.end())
                continue;
            oneringarea += (*it).second;
            //oneringarea += mesh_.property(farea_, *vfit);
        }

        //mesh_.property(v1ringarea_, *vit) = oneringarea;
        std::pair<int, float> p1;
        p1.first = (*vit).idx();
        p1.second = oneringarea;
        oneRingAreas.insert(p1);
    }
    return oneRingAreas;
}

//-----------------------------------------------------------------------------
void SkelatonExtractor::edge_collapses(MatrixXf vpos)
{
    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    int n = mesh_.n_vertices();

    int collapsedVertexId = -1;
    
    float collapsedCost = 0.;
    for (int i = 0; i < n; ++i)
        collapse(collapsedVertexId, collapsedCost);

    cout << "edge collapse process is finished" << endl;
}
//-----------------------------------------------------------------------------
void SkelatonExtractor::collapse(int& vid, float& collapsedCost)
{
    
    for (Mesh::VertexIter vit = mesh_.vertices_begin(); vit != mesh_.vertices_end(); ++vit)
    {
        float fmincost = INT_MAX;
        int maxHID = -1;
        
        for (Mesh::VertexEdgeIter veit = mesh_.ve_begin(*vit); veit != mesh_.ve_end(*vit); ++veit)
        {
            if (vid != -1)
            {
                if (mesh_.to_vertex_handle(mesh_.halfedge_handle(*veit, 0)).idx() != vid &&
                    mesh_.from_vertex_handle(mesh_.halfedge_handle(*veit, 0)).idx() != vid)
                    {
                        continue;
                    }   
            }
            Mesh::Point vj = mesh_.point(mesh_.to_vertex_handle(mesh_.halfedge_handle(*veit, 0)));
            Mesh::Point vi = mesh_.point(mesh_.from_vertex_handle(mesh_.halfedge_handle(*veit, 0)));

            float fcost1 = fabs(get_edge_cost(vi, vj, vi, vj));
            if (fcost1 <= fmincost)
            {
                fmincost = fcost1;
                maxHID = (mesh_.halfedge_handle(*veit, 0)).idx();
            }
        }
        for (Mesh::VertexIHalfedgeIter vhit = mesh_.vih_begin(*vit); vhit != mesh_.vih_end(*vit); ++vhit)
        {
            if ((*vhit).idx() == maxHID)
            {
                //collapsedhalfedges.push_back(*vhit);
                mesh_.collapse(*vhit);
                collapsedCost += fmincost;
                break;
            }
        }
        vid = (*vit).idx();
        mesh_.update_vertex_normals();
        mesh_.update_face_normals();
        break;
    }
}
//-----------------------------------------------------------------------------
float SkelatonExtractor::get_edge_cost(Mesh::Point v1, Mesh::Point v2, Mesh::Point vi, Mesh::Point pj)
{
    MatrixXf Ki(3, 4);
    VectorXf edgev(3);
    edgev(0) = v2[0] - v1[0];
    edgev(1) = v2[1] - v1[1];
    edgev(2) = v2[2] - v1[2];

    float normp = std::sqrt(pow(edgev(0), 2) + pow(edgev(1), 2) + pow(edgev(2), 2));
    return normp;

    /* Mesh::Point ai = Mesh::Point(edgev(0) / normp, edgev(1) / normp, edgev(2) / normp);

    Mesh::Point bi = Mesh::Point(ai[1] * vi[2] - ai[2] * vi[1], -ai[0] * vi[2] + ai[2] * vi[0], ai[0] * vi[1] - ai[1] * vi[0]);

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
    MatrixXf qi = Ki.transpose() * Ki;
    Vector3f pjm(pj[0], pj[1], pj[2]);
    Vector4f pjh = pjm.homogeneous();
    MatrixXf f1 = pjh.transpose() * qi;
    MatrixXf fi = f1 * pjh;
    return fi(0, 0); */
}

//=============================================================================

//=============================================================================