#include <io/read_msh.hpp>
#ifdef USE_GMSH
#include <gmsh.h>
#endif
#include <vector>
#include <map>
#include <string>
#include <fstream>
#ifdef USE_GMSH

int num_nodes_per_element(ElementType elementType)
{
    switch (elementType)
    {
    case ElementType::LINE:
        return 2; // line
    case ElementType::TRIANGLE:
        return 3; // triangle
    case ElementType::QUADRANGLE:
        return 4; // quad
    case ElementType::TRTRAHEDRON:
        return 4; // tetra
    case ElementType::HEXAHEDRON:
        return 8; // hex
    case ElementType::PRISM:
        return 6; // prism
    case ElementType::PYRAMID:
        return 5; // pyramid
    case ElementType::POINT:
        return 1; // node
    default:
        std::cerr << "Unsupported element type: " << elementType << std::endl;
        exit(1);
    }
}

void Gmsh41::read_gmsh(const char *str)
{

    gmsh::initialize();
    gmsh::open(str);

    std::vector<std::size_t> nodeTags;
    std::vector<double> coords, params;
    gmsh::model::mesh::getNodes(nodeTags, coords, params);

    int hdim = getHighestElementDimension();

    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);
    nodes.push_back({0,0,0,0});

    for (auto &e : entities)
    {
        int entityDim = e.first;
        int entityTag = e.second;

        std::vector<std::size_t> nodeTags;
        std::vector<double> coords, params;

        gmsh::model::mesh::getNodes(nodeTags, coords, params,
                                    entityDim, entityTag);

        GNodeBlock block;
        block.entityDim = entityDim;
        block.entityTag = entityTag;
        block.parametri = 0;
        block.numNodesInBlock = nodeTags.size();
        for (size_t i = 0; i < nodeTags.size(); i++)
        {
            double x = coords[3 * i];
            double y = coords[3 * i + 1];
            double z = coords[3 * i + 2];
            // printf("entity(%d,%d)  node[%zu] (%f,%f,%f)\n", entityDim, entityTag, nodeTags[i], x, y, z);
            block.nodes.push_back({nodeTags[i], x, y, z});
            nodes.push_back({nodeTags[i], x, y, z});
        }

        node_blocks.push_back(std::move(block));
    }

    for (auto &ent : entities)
    {
        int32_t dim = ent.first;
        int32_t tag = ent.second;

        GEntityBlock block;

        std::vector<int> elemTypes;
        std::vector<std::vector<std::size_t>> elemTags, elemNodes;

        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodes,
                                       dim, tag);

        if (elemTypes.empty())
            continue;
        block.entityDim = dim;
        block.entityTag = tag;
        block.elementType = elemTypes.back();

        // std::cout << "Entity (dim=" << dim
        //           << ", tag=" << tag << ",type=" << elemTypes.front() << ")\n";

        for (size_t t = 0; t < elemTypes.size(); t++)
        {

            ElementType type =  (ElementType)elemTypes[t];
            int nper = num_nodes_per_element(type);

            const auto &tags = elemTags[t];
            const auto &nodes = elemNodes[t];
            block.numElementsInBlock = tags.size();
            for (size_t e = 0; e < tags.size(); e++)
            {
                GElement el;
                std::vector<size_t> be;
                el.elementTag = tags[e];
                el.numNodesinElement = nper;
                for (int j = 0; j < nper; j++)
                {
                    // std::cout << nodes[e * nper + j] << " ";
                    el.nodes.push_back(nodes[e * nper + j]);
                    be.push_back(nodes[e * nper + j]);
                }
                if (dim == hdim)
                {
                    basic_elements.push_back(be);
                    vmap.push_back(tags[e]);
                }

                block.elements.push_back(std::move(el));
            }
            assert(block.elements.size() == tags.size());
        }

        entity_blocks.push_back(std::move(block));
    }

    gmsh::finalize();

}
#endif
