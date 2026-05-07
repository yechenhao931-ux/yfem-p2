#ifndef MESH_HPP
#define MESH_HPP

enum class GeometryType : int {
    POINT       = 0,
    SEGMENT     = 1,
    TRIANGLE    = 2,
    SQUARE      = 3,
    TETRAHEDRON = 4,
    CUBE        = 5,
    PRISM       = 6
};
#endif