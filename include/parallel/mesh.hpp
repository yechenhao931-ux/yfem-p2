#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include <iostream>
// 假设元素只含有三角形
class Element{
public:
    int vertice[3]; //三角形顶点
    Element() = default;
    Element(int v0, int v1, int v2){
        vertice[0] = v0;
        vertice[1] = v1;
        vertice[2] = v2;
    }
    ~Element() = default;
    
    void Print(){
        printf("Tri:{%d, %d, %d}",vertice[0],vertice[1],vertice[2]);
    }
};


class Mesh
{
public:
    int num_vertices;
    std::vector<Element> elements;
public:
    Mesh() = default;
    Mesh(int num_vertices_){
        num_vertices = num_vertices_;
        elements.resize(num_vertices_);
    }
    void Print(){
        printf("Mesh\n");
        for(auto &e:elements){
            e.Print();
        }
    }
    ~Mesh() = default;
};


#endif