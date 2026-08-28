CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -std=c++17 -O3 -Wall -Wextra
LDFLAGS ?=
LDLIBS ?=

TARGET := kmerstrd
SRCS := src/main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET)
