#ifndef MESH_VIEWER_WIDGET_HH
#define MESH_VIEWER_WIDGET_HH
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMeshT.hh>
#include <OpenMesh/Core/Utils/BaseProperty.hh>


#include <fstream>
#include <iostream>
#include <ostream>
#include <numeric>
#include "include\Eigen\Sparse"
#include "include\Eigen\Dense"
#include "spdiags.hpp"
#include "include\Eigen\Geometry"

#if _ITERATOR_DEBUG_LEVEL == 0 && _SECURE_SCL != 0 
#error _SECURE_SCL != 0 when _ITERATOR_DEBUG_LEVEL == 0 
#endif 

#ifndef PI
#define PI 3.141592653589793238462643383279
#endif

using namespace std;
using namespace Eigen;
class MeshViewer
{
public:
    MeshViewer();
    ~MeshViewer() {}

    virtual bool    open_mesh(const char *_filename);

protected:
    typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

protected:

    void                        update_face_indices();

    Mesh                        mesh_;

    int                         ninterior;
    int                         nboundary;
    bool                        youWantToPaintEachVertexDifferently;
    string                      fileName;
    std::vector<unsigned int>   indices_;
};

#endif // MESH_VIEWER_WIDGET_HH
