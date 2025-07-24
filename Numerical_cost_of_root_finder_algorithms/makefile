# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g -std=c++11 $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


# Source directory
SRCDIR = src
# Source files (inside src/)
SRCS = NumericalCost.cc core_statmech_functions.cc fermi_gas.cc bose_gas.cc boltzmann_gas.cc chemical_potential.cc
SRCS := $(addprefix $(SRCDIR)/, $(SRCS))

# Build directory
OBJDIR = build
# Object files (inside build/)
OBJS = $(SRCS:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

# Target executable
TARGET = NumericalCost


# Default rule to build and run the target executable
all: $(TARGET)

# Rule to link object files into the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Rule to compile .cc files into .o files inside /build
$(OBJDIR)/%.o: $(SRCDIR)/%.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Make sure the directory exists, i.e., create it if it doesn't
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Rule to run the target executable
run: $(TARGET)
	./$(TARGET)

# Clean rule to remove generated files
clean:
	rm -f $(OBJS)
