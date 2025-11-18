# === Compiler and Flags ===
CXX := g++
CXXFLAGS := -Wall -std=c++17 -O3 -flto
LDFLAGS := -lz -pthread -lhts

# === Sources ===
SRC_DIR := src
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

# === Targets ===
TARGETS := kmer_searcher

all: $(TARGETS)

kmer_searcher: $(SRCS)
	$(CXX) $(CXXFLAGS) -DCOUNTINT=uint8_t $^ -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: all clean

