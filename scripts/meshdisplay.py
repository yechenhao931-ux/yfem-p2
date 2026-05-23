import re
import sys
from pathlib import Path
import numpy as np
import struct
import plotly.graph_objects as go


from enum import IntEnum

class GeometryType(IntEnum):
    POINT       = 0
    SEGMENT     = 1
    TRIANGLE    = 2
    SQUARE      = 3
    TETRAHEDRON = 4
    CUBE        = 5
    PRISM       = 6


class MFEMMesh:
    def __init__(self, file_path):
        self.dimension = 0
        self.num_vertices = 0
        self.elements = []
        self.boundary = []
        self.vertices = []
        self._load(file_path)

    def _load(self, file_path):
        with open(file_path, 'r') as f:
            lines = f.readlines()

        current_section = None
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # 跳过空行和注释
            if not line or line.startswith('#'):
                i += 1
                continue

            # 解析关键标识符
            if line == "MFEM mesh v1.0":
                pass
            elif line == "dimension":
                self.dimension = int(lines[i+1].strip())
                i += 1
            elif line == "elements":
                num_elements = int(lines[i+1].strip())
                i += 2
                for _ in range(num_elements):
                    # 格式: attribute geometry_type v1 v2 ...
                    parts = list(map(int, lines[i].split()))
                    self.elements.append({
                        'attr': parts[0],
                        'geom': parts[1],
                        'connectivity': parts[2:]
                    })
                    i += 1
                continue
            elif line == "boundary":
                num_bdr = int(lines[i+1].strip())
                i += 2
                for _ in range(num_bdr):
                    parts = list(map(int, lines[i].split()))
                    self.boundary.append({
                        'attr': parts[0],
                        'geom': parts[1],
                        'connectivity': parts[2:]
                    })
                    i += 1
                continue
            elif line == "vertices":
                num_vertices = int(lines[i+1].strip())
                self.num_vertices = num_vertices
                # 找到下一处非空行：数字 => 普通坐标；"nodes" => 高阶网格
                j = i + 2
                while j < len(lines) and lines[j].strip() == "":
                    j += 1
                nxt = lines[j].strip() if j < len(lines) else ""
                if nxt == "nodes":
                    # 坐标在后面的 nodes GridFunction 中给出，交由 nodes 分支解析
                    i = j
                    continue
                v_dim = int(nxt)
                i = j + 1
                for _ in range(num_vertices):
                    coords = list(map(float, lines[i].split()))
                    # 补齐到三维，便于统一绘制
                    while len(coords) < 3:
                        coords.append(0.0)
                    self.vertices.append(coords)
                    i += 1
                continue
            elif line == "nodes":
                # 高阶（曲边）网格：节点坐标以 GridFunction 形式给出
                #   FiniteElementSpace
                #   FiniteElementCollection: <name>
                #   VDim: <vd>
                #   Ordering: <ord>      (0 = byNODES, 1 = byVDIM)
                #   <values...>
                # MFEM 约定前 num_vertices 个自由度即为网格顶点，据此恢复顶点坐标。
                vd, ordering = 3, 0
                i += 1
                while i < len(lines):
                    l = lines[i].strip()
                    if l.startswith("VDim"):
                        vd = int(l.split(":")[1])
                    elif l.startswith("Ordering"):
                        ordering = int(l.split(":")[1])
                    elif l == "" or l == "FiniteElementSpace" or l.startswith("FiniteElementCollection"):
                        pass
                    else:
                        break  # 到达第一个数值行
                    i += 1
                vals = []
                while i < len(lines):
                    l = lines[i].strip()
                    if l:
                        vals.extend(map(float, l.split()))
                    i += 1
                ndof = len(vals) // vd if vd else 0
                self.vertices = []
                for vi in range(self.num_vertices):
                    if ordering == 0:  # byNODES
                        coord = [vals[d * ndof + vi] for d in range(vd)]
                    else:              # byVDIM
                        coord = [vals[vi * vd + d] for d in range(vd)]
                    while len(coord) < 3:
                        coord.append(0.0)
                    self.vertices.append(coord)
                continue
            i += 1

    def info(self):
        print(f"维度: {self.dimension}")
        print(f"单元数量: {len(self.elements)}")
        print(f"边界元素数量: {len(self.boundary)}")
        print(f"顶点数量: {len(self.vertices)}")


def read_parts(file: str)->tuple[int, list[int]]:
    with open(file, 'rb') as file:
        # 读取整数
        n = struct.unpack('Q', file.read(8))[0]
        c = struct.unpack('Q', file.read(8))[0]
        
        parts = []
        format_string = f'{c}i'   # 'i' 表示有符号整数
        data = struct.unpack(format_string, file.read(c * 4))
        parts.extend(data)
    return (n,parts)


    
COLORS = [
    '#e6194b', '#3cb44b', '#4363d8', '#f58231',
    '#911eb4', '#42d4f4', '#f032e6', '#bfef45',
    '#fabed4', '#469990', '#dcbeff', '#9A6324',
]

def draw_vertices(fig, mesh: MFEMMesh):

    # --- 節点（頂点）のオーバーレイ表示 ---
    points = mesh.vertices
    x =  [sub[0] for sub in mesh.vertices]
    y =  [sub[1] for sub in mesh.vertices]
    z =  [sub[2] for sub in mesh.vertices]

    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(size=2, color='black', opacity=0.5),
        name='Nodes'
    ))

    return

# 各几何类型的棱边（局部顶点索引对）
EDGES_BY_GEOM = {
    GeometryType.TRIANGLE:    [(0, 1), (1, 2), (2, 0)],
    GeometryType.SQUARE:      [(0, 1), (1, 2), (2, 3), (3, 0)],
    GeometryType.TETRAHEDRON: [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    GeometryType.CUBE:        [(0, 1), (1, 2), (2, 3), (3, 0),
                               (4, 5), (5, 6), (6, 7), (7, 4),
                               (0, 4), (1, 5), (2, 6), (3, 7)],
}


def draw_elements(fig, mesh: MFEMMesh, part_file):

    (N, parts) = read_parts(part_file)
    vertices = mesh.vertices

    all_lines = {p: [[], [], []] for p in range(N)}

    for idx, e in enumerate(mesh.elements):
        part = parts[idx]
        edges = EDGES_BY_GEOM.get(e['geom'])
        if edges is None:
            continue  # 暂不支持的几何类型（如 PRISM）
        v = e['connectivity']
        for a, b in edges:
            va, vb = v[a], v[b]
            all_lines[part][0] += [vertices[va][0], vertices[vb][0], None]
            all_lines[part][1] += [vertices[va][1], vertices[vb][1], None]
            all_lines[part][2] += [vertices[va][2], vertices[vb][2], None]

    # draw
    for p in range(N):
        fig.add_trace(go.Scatter3d(
            x=all_lines[p][0], y=all_lines[p][1], z=all_lines[p][2],
            mode='lines', line=dict(color=COLORS[p % len(COLORS)]),
        ))

    return


if __name__ == "__main__":
    # 替换为你实际的 mesh 文件路径
    mesh_file = "../resource/box.mesh" 


    if len(sys.argv) >= 2:
        mesh_file =  sys.argv[1]
    


    path_str = mesh_file
    path_obj = Path(path_str)

    # stem 属性可以直接获取“文件名去掉后缀”的部分
    part_file  = "../build/"+ path_obj.stem + ".part"


    try:
        mesh = MFEMMesh(mesh_file)
        mesh.info()


        fig = go.Figure()

        draw_vertices(fig, mesh)
        draw_elements(fig, mesh, part_file)

        fig.update_layout(
            title=f"mesh:{mesh_file}",
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Z',
                aspectmode='data'
            ),
            legend=dict(
                title="Physical Groups",
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='gray',
                borderwidth=1
            ),
            margin=dict(l=0, r=0, b=0, t=50)
        )

        fig.show()
            
    except FileNotFoundError:
        print("未找到文件，请检查路径。")

    import re


