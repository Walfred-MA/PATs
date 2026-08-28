CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -Wall -std=c++17 -O3
LDFLAGS ?=
LDLIBS ?= -lhts -lz -pthread

# === Sources ===
SRC_DIR := src
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

# === Targets ===
TARGETS := kmer_convertor kmer_convertor8

all: $(TARGETS)

kmer_convertor: $(SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOUNTINT=uint16_t $^ -o $@ $(LDFLAGS) $(LDLIBS)

kmer_convertor8: $(SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DCOUNTINT=uint8_t $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: all clean
