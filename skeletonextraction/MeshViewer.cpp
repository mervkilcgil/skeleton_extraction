#include "MeshViewer.h"

MeshViewer::MeshViewer()
{
    ninterior = 0;
    nboundary = 0;
    youWantToPaintEachVertexDifferently = false;
    fileName = ".off";
}

//-----------------------------------------------------------------------------
bool MeshViewer::open_mesh(const char *_filename)
{
    // load mesh
    fileName = string(_filename);
    if (OpenMesh::IO::read_mesh(mesh_, _filename))
    {
        mesh_.request_face_normals();
        mesh_.request_vertex_normals();
        mesh_.request_vertex_status();
        mesh_.request_vertex_texcoords2D();
        mesh_.request_edge_status();
        mesh_.request_edge_colors();
        mesh_.request_face_colors();
        mesh_.request_face_status();
        mesh_.request_face_texture_index();
        mesh_.request_halfedge_status();
        mesh_.request_halfedge_normals();
        mesh_.request_halfedge_texcoords1D();
        mesh_.request_halfedge_texcoords2D();
        mesh_.request_halfedge_texcoords3D();
        mesh_.request_vertex_colors();
        mesh_.request_vertex_texcoords1D();
        mesh_.request_vertex_texcoords3D();
        mesh_.update_normals();
        update_face_indices();

        Mesh::VertexIter viter;
        for (viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter)
        {
            if (mesh_.is_boundary(*viter))
                ++nboundary;
            else
                ++ninterior;
        }

    }

    return false;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void MeshViewer::update_face_indices()
{
    Mesh::ConstFaceIter f_it,
        f_end(mesh_.faces_end());
    Mesh::ConstFaceVertexIter fv_it;

    indices_.clear();
    indices_.reserve(mesh_.n_faces() * 3);
    std::cout << "mesh indices updated" << std::endl;

    for (f_it = mesh_.faces_begin(); f_it != f_end; ++f_it)
    {
        fv_it = mesh_.cfv_iter(*f_it);
        Mesh::VertexHandle id1 = *fv_it;
        Mesh::VertexHandle id2 = *(++fv_it);
        Mesh::VertexHandle id3 = *(++fv_it);

        if (id1.is_valid() && id2.is_valid() && id3.is_valid()){
            indices_.push_back(id1.idx());
            indices_.push_back(id2.idx());
            indices_.push_back(id3.idx());
        }
        
    }
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//=============================================================================
