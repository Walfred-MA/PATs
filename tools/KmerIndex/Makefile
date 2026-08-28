CXX ?= c++
CPPFLAGS ?=
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDFLAGS ?=
LDLIBS ?= -pthread

TARGET := kmerindex
SRC := src/main.cpp src/kmer_index.cpp
OBJ := $(SRC:.cpp=.o)
DEP := $(OBJ:.o=.d)

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

src/%.o: src/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

test: $(TARGET)
	python3 tests/run_tests.py ./$(TARGET)

clean:
	rm -f $(TARGET) $(OBJ) $(DEP)

-include $(DEP)
