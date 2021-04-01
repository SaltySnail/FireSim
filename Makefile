# A simple Makefile for compiling small projects

# set the compiler
CC := g++

# set the compiler flags
CPPFLAGS := -std=c++17 -Wall -Wpedantic `sdl2-config --cflags --libs` -lSDL2_image -lm -I include/

# add header files here
HDRS := FastNoiseLite.h

# add source files here
SRCS := main.cpp

# generate names of object files
OBJS := $(SRCS:.cpp=.o)

# name of executable
EXEC := fastnoise

# default recipe
all: clean $(EXEC)

# recipe for building the final executable
$(EXEC): $(OBJS) $(HDRS) Makefile
	$(CC) -o $@ $(OBJS) $(CPPFLAGS)

# recipe for building object files
#$(OBJS): $(@:.o=.c) $(HDRS) Makefile
#	$(CC) -o $@ $(@:.o=.c) -c $(CFLAGS)

# recipe to clean the workspace
clean:
	rm -f $(wildcard *.o) $(EXEC)
