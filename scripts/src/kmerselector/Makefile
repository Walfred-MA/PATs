# Compiler and flags
CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra 
LDFLAGS = -lz -pthread

# Sources and target
SRCS = src/fasta.cpp src/gzfile.cpp src/KmerCounter.cpp src/KmerFilter.cpp src/KmerHash.cpp src/KmerStruct.cpp src/main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = kmerselector

# Default rule
all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✅ Build complete: $(TARGET)"

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
	@echo "🧹 Cleaned up."

.PHONY: all clean

