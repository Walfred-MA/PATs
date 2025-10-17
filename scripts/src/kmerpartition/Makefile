# =========================================================
# Makefile for kmerpartition
# Author: Wangfei Ma
# =========================================================

# Compiler and flags
CXX      := g++
CXXFLAGS := -O3 -std=c++17 -fopenmp -Wall -Wextra -mtune=native
LDFLAGS  := -fopenmp

# Directories
SRC_DIR  := src
OBJ_DIR  := obj
BIN_DIR  := .

# Target executable
TARGET   := $(BIN_DIR)/kmerpartition

# Source and object files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Include paths
INCLUDES := -I$(SRC_DIR)

# =========================================================
# Build rules
# =========================================================

all: $(TARGET)

$(TARGET): $(OBJ_DIR) $(OBJS)
	@echo "🔗 Linking $(TARGET)..."
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)
	@echo "✅ Build complete."

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "⚙️  Compiling $< ..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# =========================================================
# Housekeeping
# =========================================================

clean:
	@echo "🧹 Cleaning build artifacts..."
	rm -rf $(OBJ_DIR) $(TARGET)

rebuild: clean all

.PHONY: all clean rebuild

