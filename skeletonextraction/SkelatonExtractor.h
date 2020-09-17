#ifndef SKELETON_EXTRACTOR_HH
#define SKELETON_EXTRACTOR_HH

#include "MeshViewer.h"
#include <ctime>

class SkelatonExtractor : public MeshViewer
{
public:
    SkelatonExtractor();
    ~SkelatonExtractor();

    virtual bool            open_mesh(const char *_meshfilename);
    void                    draw();
    void                    init();
    string                  getNewFileName() { return nfname;}

protected:

    void                    make_laplacian_smooth(int numit, MatrixXf& vpos);
    void                    edge_collapses(MatrixXf vpos);

    std::map<int, float>    compute_weights();
    std::map<int, float>    compute_face_area();
    std::map<int, float>    calculate_one_ring_area();
    
    void                    get_laplacian_weight_matrix(std::map<int, float> oneringAreas,
                                                        std::map<int, float> edgeWeights,
                                                        std::vector<Triplet<float> >& L_coeff,
                                                        MatrixXf& WH0_diag, MatrixXf& vpos, 
                                                        float& originalFaceAreaSum);

    float                   get_edge_cost(Mesh::Point v1, Mesh::Point v2, Mesh::Point vi, Mesh::Point pj);

    float                   get_distance_bw_2points(Mesh::Point p1, Mesh::Point p2);
    float                   get_average_face_area();
    void                    collapse(int& vid, float& collapsedCost);

    
    OpenMesh::VPropHandleT<Mesh::Point>     vpos_;
    OpenMesh::EPropHandleT<float>           ecotanweight_;
    OpenMesh::FPropHandleT<float>           farea_;
    OpenMesh::VPropHandleT<float>           v1ringarea_;

    std::ofstream                           outfile1;

    string                                  nfname;

    
};

#endif // SKELETON_EXTRACTOR_HH
