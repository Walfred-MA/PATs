CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -Wall -std=c++17 -O3 -DCOUNTINT=uint8_t
LDFLAGS ?=
LDLIBS ?= -lhts -lcurl -lz -pthread

# === Source Files ===
SRC_DIR  := src
SRCS     := $(wildcard $(SRC_DIR)/*.cpp)
TARGET   := kmer_searcher

# === Rules ===
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
