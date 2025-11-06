CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra

TARGET := KmerStrd
SRCS := src/main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f $(TARGET)

