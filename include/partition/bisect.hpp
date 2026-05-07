#ifndef BISECT_HPP
#define BISECT_HPP

#include "../graph/graph.hpp"
#include <vector>
#include <limits>
#include <functional>

// ============================================================
//  BisectionOptions  —— 二分划分控制参数
// ============================================================
struct BisectionOptions
{
    int     nTrials     = 10;       //多次随机尝试次数
    int     nFMPasses   = 8;        //每次尝试的FM最大迭代轮数
    double  ubFactor    = 1.03;     //负载不平衡容忍系数（1.03）
    int     seed        = 42;       // 随机种子
    bool    verbose     = false;    // 打印调试信息
};

// ============================================================
//  BisectionResult  —— 返回统计结果
// ============================================================
struct BisectionResult{
    int mincut  = std::numeric_limits<int>::max();
};

#endif