# Compiler and flags
CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -MMD -MP
LDFLAGS ?=
LDLIBS ?= -lz -pthread

# Sources and target
SRCS := src/fasta.cpp src/gzfile.cpp src/KmerCounter.cpp src/KmerFilter.cpp src/KmerHash.cpp src/KmerStruct.cpp src/main.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)
TARGET := kmer_selector

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
	@echo "Build complete: $(TARGET)"

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)
	@echo "Cleaned up."
