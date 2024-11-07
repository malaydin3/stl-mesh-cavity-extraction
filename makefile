# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -O2

# Target executable
TARGET = main.out

# Source file
SRC = main.cpp

# Rule to build the target
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

# Clean rule to remove the target
clean:
	rm -f $(TARGET)
