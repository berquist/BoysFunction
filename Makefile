CC       = clang
CXX      = clang++
FLAGS    = -O3 -g -Wall -DDEBUG
CFLAGS   = $(FLAGS) -std=gnu99 -lm
CXXFLAGS = $(FLAGS) -std=c++11

all: boys.x newboys.x

%.x: %.C
	$(CXX) $(CXXFLAGS) $< -o $@

%.x: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -f *.o *.x
	$(RM) -rf *.x.dSYM
