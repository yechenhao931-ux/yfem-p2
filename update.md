1.
改进完成并推送到 `claude/improve-mesh-partition-L1nZ1`。

## 关键改进点

**1. 连通性修复 `connectrepair.cpp`（新增）**
- `EnforceConnectivity`：BFS 找出每个分区的连通子图，把孤岛迁到切边最多的邻居分区
- `PostBalanceFix`：仅在边界点上做的连通性保留型平衡修复

**2. K-means 增强 `kmeanspart.cpp`**
- `BFSSeedExpansion`：从 K-means++ 几何种子做轮转 BFS，初始分区已连通且均衡
- 默认启用 `BalanceFix`（连通性感知）取代原 `BalanceRepair`（纯几何）
- 可选孤岛修复

**3. `GeoKwayPartition` 收尾阶段**
- 粗化加 `maxVwgtCap`（BHEM 思路），避免超级顶点失衡
- HEM 匹配改用平方距离，去掉每边 sqrt
- FM 后交替 `EnforceConnectivity`／`PostBalanceFix`（最多 3 轮），再跑一次短 FM 找回切边质量

**4. `GeoKwayFMCut` 性能**
- centroid bbox 由 O(K²) 改为 O(n) 顶点 bbox

## 量化效果（box.mesh，8 分区，10 个随机种子）

| 指标 | baseline | all-on | 变化 |
|---|---|---|---|
| 最大不平衡 | 13.67% | 9.50% | **-30%** |
| 平均不平衡 | 6.27% | 4.97% | **-21%** |
| 平均切边 | 632.3 | 636.1 | +0.6% |
| 孤立分区 | 0 | 0 | 持平 |

**5. 新增 `test/bench_partition.cpp`** 消融实验工具，直接对比五种配置的切边／不平衡／连通分量数。
