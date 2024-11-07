# Cavity Extraction Tool for STL Meshes

## Overview

This application identifies and extracts internal cavities from a 3D STL mesh file. It processes the mesh to find disconnected regions (or cavities) that are completely enclosed within the surface. The program computes adjacency relationships between triangles, identifies cavities through traversal, and outputs these regions in a new STL file.

This tool was developed as part of a computational geometry challenge, using C++ and the `stl_reader` library.

## Features

- **Identify Cavities**: Detects internal disconnected regions (cavities) in a 3D mesh.
- **Save Results**: Outputs the extracted (individual) cavities as an STL file.
- **Efficient Traversal**: Uses BFS and DFS for cavity detection.
- **Automatic Orientation Correction**: Ensures triangle normals are consistently oriented.

## Prerequisites

- **stl_reader**: Ensure the `stl_reader.h` library is available in the project directory. This library allows easy access to mesh triangle vertices and connectivity from STL files.
- **C++11** or later for standard library support.

## Getting Started

1. **Compile the Program**:
   ```bash
   g++ main.cpp -o main.out

2. **Run the Program**:
   ```bash
   ./main.out
