# Compiler
CXX = g++

#Compiler flags
CXXFLAGS = -Wall -g
#CPPFLAGS = -I/home/path/gsl/include
#LDFLAGS =  -L/home/path/gsl/lib
LDLIBS = -lgsl -lgslcblas -lm
#Target executable
TARGET = NumericalCost

#Source files
SRCS = NumericalCost.cc core_statmech_functions.cc fermi_gas.cc bose_gas.cc boltzmann_gas.cc chemical_potential.cc

#Object files
OBJS = $(SRCS:.cc=.o)

#Default rule to build and run the executable
all: $(TARGET)

# Rule to link object files into the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Rule to compile .cc files into .o files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to run the executable
run: $(TARGET)
	./$(TARGET)

# Clean rule to remove generated files
clean:
	rm -f $(TARGET) $(OBJS)