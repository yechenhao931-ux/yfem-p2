#include <graph/dual_graph.hpp>
#include <iostream>
#include <cassert>
typedef std::vector<size_t> Face;


std::vector<Face> select_main_element(std::vector<size_t>& vertices,int dim){
    std::vector<Face> faces;
    assert(0);
    switch (dim)
    {
    case 2: //2D
        
    case 3: //3D
    default:
        break;
    }

}

std::vector<std::vector<size_t>> build_dual_graph(const std::vector<std::vector<size_t>>& elements, int dim){
    size_t ne = elements.size();
    std::cout << ne << "\n";
    std::vector<std::vector<size_t>> dual_adj(ne);

    std::unordered_map<Face, size_t, FaceHash> face2elem;

    for(size_t ei = 0; ei < ne; ei++){
        const auto& vertices = elements[ei];
        std::vector<Face> faces;

        // 处理不同类型
        if(dim == 2){
            // triangle
            if(vertices.size() == 3){
                faces = {
                    {vertices[0],vertices[1]},
                    {vertices[1],vertices[2]},
                    {vertices[2],vertices[0]}
                };
            }
            // square
            else if (vertices.size() == 4){
                faces = {
                    {vertices[0], vertices[1]},
                    {vertices[1], vertices[2]},
                    {vertices[2], vertices[3]},
                    {vertices[3], vertices[0]}
                };
            } else{
                assert(0);
            }
        }else if(dim == 3){
            // tetrahedron (4 个三角面)
            if (vertices.size() == 4){
                faces = {
                    {vertices[0], vertices[1], vertices[2]},
                    {vertices[0], vertices[1], vertices[3]},
                    {vertices[1], vertices[2], vertices[3]},
                    {vertices[0], vertices[2], vertices[3]}
                };
            }
            // hexahedron / cube (6 个四边形面，MFEM 顶点序)
            else if (vertices.size() == 8){
                faces = {
                    {vertices[0], vertices[1], vertices[2], vertices[3]}, // bottom
                    {vertices[4], vertices[5], vertices[6], vertices[7]}, // top
                    {vertices[0], vertices[1], vertices[5], vertices[4]},
                    {vertices[1], vertices[2], vertices[6], vertices[5]},
                    {vertices[2], vertices[3], vertices[7], vertices[6]},
                    {vertices[3], vertices[0], vertices[4], vertices[7]}
                };
            } else{
                assert(0);
            }
        }
        

        for(auto& f: faces){
            std::sort(f.begin(),f.end());
            auto it = face2elem.find(f);
            if(it == face2elem.end()){
                face2elem[f] = ei;
            } else {
                size_t ej = it->second;
                dual_adj[ei].push_back(ej);
                dual_adj[ej].push_back(ei);
            }
        }

    }

    return dual_adj;
}