#include <io/save.hpp>


void SavePartitionResult(const Graph& g, size_t n_part, std::string path){
    std::ofstream file(path, std::ios::binary);
    
    if (!file) {
        std::cerr << "无法打开文件！" << std::endl;
        return;
    }
    
    
    size_t n = g.nvtxs;
    
    // 以二进制格式写入数据
    printf("res: %ld \n",g.where.size());
    file.write(reinterpret_cast<char*>(&n_part), sizeof(n));
    file.write(reinterpret_cast<char*>(&n), sizeof(n));

    for(const int p: g.where){
        file.write((const char *)(&p), sizeof(p));
    }

    file.close();
}