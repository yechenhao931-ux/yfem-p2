# --- 路径定义 ---
MFEM_DIR   := $(HOME)/mfem/mfem-4.9p
HYPRE_DIR  := $(HOME)/mfem/hypre/src/hypre
METIS_DIR  := $(HOME)/mfem/metis-4.0
YFEM_DIR   := $(HOME)/yfem-op1

# --- 编译选项 ---
CXX      := mpic++
CXXFLAGS := -std=c++17 -O3
INCLUDES := -I$(YFEM_DIR)/include -I$(HYPRE_DIR)/include -I$(MFEM_DIR)
LDFLAGS  := -L$(MFEM_DIR) -L$(HYPRE_DIR)/lib -L$(METIS_DIR)
LIBS     := -lmfem -lHYPRE -lmetis -lrt

# --- YFEM 静态库 ---
YFEM_LIBS := $(YFEM_DIR)/build/src/common/libcommon.a \
             $(YFEM_DIR)/build/src/graph/libgraph.a \
             $(YFEM_DIR)/build/src/io/libio.a \
             $(YFEM_DIR)/build/src/partition/libpartition.a

# --- 目标文件 ---
BUILD_DIR := build

.PHONY: all dir origin op

all: dir origin op

dir:
	mkdir -p $(BUILD_DIR)

origin:
	$(CXX) $(CXXFLAGS) test_mfem.cpp -o $(BUILD_DIR)/test_mfem \
		-I$(HYPRE_DIR)/include -I$(MFEM_DIR) \
		-L$(MFEM_DIR) -L$(HYPRE_DIR)/lib -L$(METIS_DIR) \
		-lmfem -lHYPRE -lmetis -lrt

op:
	$(CXX) $(CXXFLAGS) test_mfemm.cpp -o $(BUILD_DIR)/test_mfemm \
		$(INCLUDES) $(LDFLAGS) $(LIBS) $(YFEM_LIBS)

clean:
	rm -rf $(BUILD_DIR)
