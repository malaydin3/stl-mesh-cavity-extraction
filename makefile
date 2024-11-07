# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -O2

# Target executable
TARGET = problem_1.out

# Source file
SRC = problem_1.cpp

# Rule to build the target
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Clean rule to remove the target
clean:
	rm -f $(TARGET)
