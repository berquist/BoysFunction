CC       = clang
CXX      = clang++
FLAGS    = -O3 -g -Wall
CFLAGS   = $(FLAGS) -std=c99
CXXFLAGS = $(FLAGS) -std=c++11

all: boys.x

%.x: %.C
	$(CXX) $(CXXFLAGS) $< -o $@

%.x: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -f *.o *.x
