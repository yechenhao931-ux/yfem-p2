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
                v_dim = int(lines[i+2].strip())
                i += 3
                for _ in range(num_vertices):
                    coords = list(map(float, lines[i].split()))
                    self.vertices.append(coords)
                    i += 1
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

def draw_elements(fig, mesh: MFEMMesh, part_file):

   #(N,parts) = read_parts(part_file)
    
    (N,parts) = read_parts(part_file)
    vertices = mesh.vertices

    all_lines = {}
    for i in range(N):
        all_lines[i] = [[],[],[]]

    for i, e in enumerate(mesh.elements):
        part = parts[i]
        if e['geom'] == 4:
            v = e['connectivity'] 
            i, j, k, l = v[0], v[1], v[2], v[3]
            for a, b in [(i, j), (j, k), (k, i), (i,l), (j,l), (k,l)]:
                all_lines[part][0] += [vertices[a][0], vertices[b][0], None]
                all_lines[part][1] += [vertices[a][1], vertices[b][1], None]
                all_lines[part][2] += [vertices[a][2], vertices[b][2], None]


    # draw
    for i in range(N):
        fig.add_trace(go.Scatter3d(
            x=all_lines[i][0], y=all_lines[i][1], z=all_lines[i][2],
            mode='lines',line=dict(color=COLORS[i]),
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


