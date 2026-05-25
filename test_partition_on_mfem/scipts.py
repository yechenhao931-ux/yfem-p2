import subprocess
import re
import statistics

N_RUNS = 20
NP = 8
EXEC = "./build/test_mfem"
MESH = "../resource/box.mesh"

times = []

pattern = re.compile(r"FEM solve time:\s*([0-9.]+)")

for i in range(N_RUNS):
    print(f"Run {i+1}/{N_RUNS}")

    cmd = [
        "mpirun",
        "-np", str(NP),
        EXEC,
        "-m", MESH
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True
    )

    match = pattern.search(result.stdout)

    if match:
        t = float(match.group(1))
        times.append(t)
    else:
        print("Warning: no time found")
        print(result.stdout)

# 写文件
with open("result.txt", "w") as f:
    for t in times:
        f.write(f"{t}\n")

# 统计
avg = statistics.mean(times)
min_v = min(times)
max_v = max(times)

print("\n=== Statistics ===")
print(f"avg = {avg:.6f}")
print(f"min = {min_v:.6f}")
print(f"max = {max_v:.6f}")