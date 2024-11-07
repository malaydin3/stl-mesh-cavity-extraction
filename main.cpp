// Mert Alaydin
// October 2024
// Computational Geomerty Challenge
// compile with :
// g++ main.cpp -o main.out

#include <iostream>
#include <map>
#include <queue>
#include <unordered_set>
#include "stl_reader.h"


using namespace stl_reader;
using namespace std;


// This function checks if a given normal is aligned with one the Cartesian Axis
bool isAlignedWithAxis(const vector<float>& normal, float epsilon = 1e-12) {
    int nonZeroCount = 0;
    
    // Check if each component is (approximately) zero
    if (std::abs(normal[0]) > epsilon) nonZeroCount++;
    if (std::abs(normal[1]) > epsilon) nonZeroCount++;
    if (std::abs(normal[2]) > epsilon) nonZeroCount++;
    
    // If only one component is non-zero, the vector is aligned with an axis
    return (nonZeroCount == 1);
}

// A function to calculate triangle surface normal in 3D for given XYZ coordinates of triangle vertices, v0, v1, v2
vector<float> computeTriangleNormal(const float* v0, const float* v1, const float* v2){
    // compute the triangle normal
    vector<float> normal(3,0), e1(3,0), e2(3,0); 
    // compute the first edge (v2v1)
    e1[0] = v1[0] - v0[0];
    e1[1] = v1[1] - v0[1];
    e1[2] = v1[2] - v0[2];
    // compute the second edge (v3v1)
    e2[0] = v2[0] - v0[0];
    e2[1] = v2[1] - v0[1];
    e2[2] = v2[2] - v0[2];
    // take the cross product between e1 and e2
    // Set Normal[0] to (multiply U[1] by V[2]) minus (multiply U[2] by V[1])
	// Set Normal[1] to (multiply U[2] by V[0]) minus (multiply U[0] by V[2])
	// Set Normal[2] to (multiply U[0] by V[1]) minus (multiply U[1] by V[0])
    normal[0] = (e1[1]*e2[2]) - (e1[2]*e2[1]);
    normal[1] = (e1[2]*e2[0]) - (e1[0]*e2[2]);
    normal[2] = (e1[0]*e2[1]) - (e1[1]*e2[0]);
    return normal;
}

// Function to detect all the disconnected bodies to a given subset of triangles identified previously
vector<vector<unsigned int>> detectDisconnectedBodiesInSubset(const vector<vector<unsigned int>>& adjList, const vector<unsigned int>& cavityTriangles) {
    // Disconnected bodies
    vector<vector<unsigned int>> bodies;
    // set for triangle indices
    unordered_set<unsigned int> cavitySet(cavityTriangles.begin(), cavityTriangles.end()); 
    // keep track of all visited triangles 
    unordered_set<unsigned int> visited; 

    // Traverse only the triangles in cavityTriangles
    for (unsigned int triangle : cavityTriangles) {
         // If this triangle hasn't been visited - this also includes previously exhausted bodie(s)
        if (visited.find(triangle) == visited.end()) { 
            vector<unsigned int> body;  // A new disconnected body
            queue<unsigned int> toVisit;
            // push it to our level queue
            toVisit.push(triangle);
            // mark it as visited
            visited.insert(triangle);

            // Perform BFS to find all triangles connected to this one, within the subset
            while (!toVisit.empty()) {
                unsigned int current = toVisit.front();
                toVisit.pop();
                body.push_back(current);  // Add the current triangle to this body

                // Visit all adjacent triangles
                for (unsigned int neighbor : adjList[current]) {
                    // check if this triangle belongs to overall cavity set and has not been visited before 
                    if (cavitySet.find(neighbor) != cavitySet.end() && visited.find(neighbor) == visited.end()) {
                        visited.insert(neighbor);
                        toVisit.push(neighbor);
                    }
                }
            }

            // at this point we have exhausted all the triangles that are connected to each other in the first disjoint set 
            // Store the connected component (body)
            bodies.push_back(body);
        }
    }

    std::cout << "Total number of individual cavities found: " << bodies.size() << endl;
    return bodies; 
}



// Conduct a DFS to find all the connected triangles to a given triangle
void traverseMesh(int triangle, vector<vector<unsigned int>>& adjList, vector<bool>& visited) {
    // Mark the current triangle as visited
    visited[triangle] = true;
    
    // Recursively visit all neighbor triangles
    for (unsigned int neighbor : adjList[triangle]) {
        if (!visited[neighbor]) {
            traverseMesh(neighbor, adjList, visited);
        }
    }
}

// Store all unvisited triangles after DFS
vector<unsigned int> findCavityTriangles(vector<vector<unsigned int>>& adjList, unsigned int startTriangle, const size_t& numTriangles) {
    vector<bool> visited(numTriangles, false);
    vector<unsigned int> cavityTriangles;

    // Traverse (DFS) all the neighbors starting from an outer surface triangle
    traverseMesh(startTriangle, adjList, visited);

    // After traversal, visited array will mark all triangles that are connected to the startTriagle and its neighbors
    for (unsigned int i = 0; i < numTriangles; i++) {
        // all the unvisited triangles are part of a cavity because they are not connected to the exterior!
        if (!visited[i]) {
            // Triangles that are not visited are part of cavities or disconnected components
            cavityTriangles.push_back(i);
        }
    }
    return cavityTriangles;
}

// Adjancecy list for conducting a DFS
vector<vector<unsigned int>> buildAdjList(const StlMesh <float, unsigned int>& mesh ) {
    size_t numTriangles = mesh.num_tris();
    vector<vector<unsigned int>> adjList(numTriangles);


    // Map to store edges and the triangle that owns that edge
    // TODO:
    // map is probably not the most efficient choice here due to red-black balance tree properites
    // one can use unordered_map with a sutiable hash function to achieve the same behaviour
    map<pair<unsigned int, unsigned int>, unsigned int> edgeToTriangle;

    // Iterate over each triangle
    for (unsigned int i = 0; i < numTriangles; i++) {
        // Extract the vertices of the triangle
        unsigned int v0 = mesh.tri_corner_ind(i,0);
        unsigned int v1 = mesh.tri_corner_ind(i,1);
        unsigned int v2 = mesh.tri_corner_ind(i,2);

        // Get the edges for this triangle, sorted to avoid direction issues
        pair<unsigned int, unsigned int> edges[3] = {
            make_pair(min(v0, v1), max(v0, v1)),
            make_pair(min(v1, v2), max(v1, v2)),
            make_pair(min(v2, v0), max(v2, v0))
        };

        // For each edge, check if it already exists in the map
        for (unsigned int e = 0; e < 3; e++) {
            auto edge = edges[e];
            // if it doesn't exist then add this triangle
            if (edgeToTriangle.find(edge) == edgeToTriangle.end()) {
                // If the edge is not yet mapped, add it with the current triangle index
                edgeToTriangle[edge] = i;
            } else {
                // If the edge is already in the map, we found a neighboring triangle
                int neighbor = edgeToTriangle[edge];

                // Add both triangles to each other's adjacency list
                // add to my neighbor's list
                adjList[i].push_back(neighbor);
                // add me to my neighbor's list 
                adjList[neighbor].push_back(i);
            }
        }
    }

    return adjList;
}

// A helper function to find a triangle that is aligned with one of the cartesian axis 
unsigned int findExteriorTriangle(const StlMesh <float, unsigned int>& mesh){
    // exterior triangle index
    unsigned int result = 0;
    for(size_t triangle = 0; triangle < mesh.num_tris(); ++triangle) {
        vector<float> normal = computeTriangleNormal(mesh.tri_corner_coords(triangle, 0), mesh.tri_corner_coords(triangle, 1), mesh.tri_corner_coords(triangle, 2));
        if(isAlignedWithAxis(normal)){
            result = triangle;
            std::cout << "Exterior triangle is " << triangle << " with normal " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
            break;
        }
    }
    // if not found then return 1
    return result;
}

// Function to write STL file
void writeSTL(const std::string& filename, const vector<unsigned int>& triangles, const StlMesh <float, unsigned int>& mesh) {
    std::ofstream file;
    file.open(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open the file!" << std::endl;
        return;
    }
    
    // Start the STL file with the header
    file << "solid I_like_cg" << std::endl;
    
    // Write each triangle to the STL file
    for (size_t i = 0; i < triangles.size(); i++) {
        // normals arrived corrupted in the original STL file so let's compute them 
        vector<float> normal = computeTriangleNormal(mesh.tri_corner_coords(triangles[i], 0), mesh.tri_corner_coords(triangles[i], 1), mesh.tri_corner_coords(triangles[i], 2));
        file << "facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
        file << "    outer loop" << std::endl;
        file << "        vertex " << mesh.tri_corner_coords(triangles[i], 0)[0] << " " << mesh.tri_corner_coords(triangles[i], 0)[1] << " " << mesh.tri_corner_coords(triangles[i], 0)[2] << std::endl;
        file << "        vertex " << mesh.tri_corner_coords(triangles[i], 1)[0] << " " << mesh.tri_corner_coords(triangles[i], 1)[1] << " " << mesh.tri_corner_coords(triangles[i], 1)[2] << std::endl;
        file << "        vertex " << mesh.tri_corner_coords(triangles[i], 2)[0] << " " << mesh.tri_corner_coords(triangles[i], 2)[1] << " " << mesh.tri_corner_coords(triangles[i], 2)[2] << std::endl;
        file << "    endloop" << std::endl;
        file << "endfacet" << std::endl;

    }
    
    // End the STL file
    file << "endsolid I_like_cg" << std::endl;
    
    file.close();
    std::cout << "STL file '" << filename << "' written successfully." << std::endl;
}



// a function to fix winding order of a given triangle 
bool orientTriangle(const unsigned int* targetTriangle, const vector<unsigned int>& commonVertices){
    int v1Index = -1, v2Index = -1;

    // Find indices of common vertices in the target triangle
    for (int i = 0; i < 3; ++i) {
        if (targetTriangle[i] == commonVertices[0]) v1Index = i;
        if (targetTriangle[i] == commonVertices[1]) v2Index = i;
    }

    // If the common vertices in the target triangle are in reverse order, we need to reorient
    if (v2Index > v1Index) {
        return true; // Reorientation is required
    }

    return false; // No re-orientation needed
}


// BFS through the seed triangle and check the orientation consistency in all connected triangles to the seed
vector<unsigned int>  orientMesh(const unsigned int& seedTriangle, const vector<vector<unsigned int>>& adjList, StlMesh<float, unsigned int>& mesh) {
    // Keep track of visited triangles 
    vector<bool> visitList(mesh.num_tris(), false);
    // Queue for BFS
    queue<unsigned int> levelSet;

    // Re-oriented triangles
    vector<unsigned int> reorientedTriangles;

    // Push the seed triangle to the queue if not visited
    if (!visitList[seedTriangle]) {
        visitList[seedTriangle] = true;
        levelSet.push(seedTriangle);
    }

    // Also let's add the seed triangle to the reoriented list
    reorientedTriangles.push_back(seedTriangle);

    // Iterate through the levels of neighborhood
    while (!levelSet.empty()) {
        // Get the triangle in the front of the queue
        unsigned int currentTriangle = levelSet.front();
        levelSet.pop();

        // Iterate through adjacent triangles of the current triangle
        for (const unsigned int& adjacentTriangle : adjList[currentTriangle]) {
            // Proceed if this triangle has not been checked before for consistent orientation
            if (!visitList[adjacentTriangle]) {
                // Find out the two shared vertices between the triangles (the common edge)
                vector<unsigned int> commonVertices;
                for (unsigned int i = 0; i < 3 && commonVertices.size() < 2; i++) {
                    for (unsigned int j = 0; j < 3; j++) {
                        if (mesh.tri_corner_inds(currentTriangle)[i] == mesh.tri_corner_inds(adjacentTriangle)[j]) {
                            commonVertices.push_back(mesh.tri_corner_inds(currentTriangle)[i]);
                            if (commonVertices.size() == 2) break; // Found both vertices
                        }
                    }
                }

                // Check the size of common vertices. If it is 2, proceed with orientation check
                if (commonVertices.size() == 2) {
                    // Mark the triangle as visited before checking orientation
                    visitList[adjacentTriangle] = true;

                    // Reorient if needed
                    bool requiredOrientation = orientTriangle(mesh.tri_corner_inds(adjacentTriangle), commonVertices);
                    if (requiredOrientation) {
                        // Swap the first and third vertex of the adjacent triangle
                        unsigned int* indices = const_cast<unsigned int*>(mesh.tri_corner_inds(adjacentTriangle));
                        // Fix the winding order by swapping first and last indices
                        std::swap(indices[0], indices[2]);
                        // record the triangle
                        reorientedTriangles.push_back(adjacentTriangle);
                    }


                    // Push the adjacent triangle onto the queue
                    levelSet.push(adjacentTriangle);
                }
            }
        }
    }

    // return the list of triangles
    return reorientedTriangles;
}


// A helper function to load the STL file
StlMesh <float, unsigned int> readSTL(const std::string& filename){

    // Create an instance of the mesh class 
    StlMesh <float, unsigned int> mesh (filename);

    return mesh;
}


// Driver function 
int main(){
    // the name of the file provided in the problem statement 
    StlMesh <float, unsigned int> mesh  = readSTL("geometry_with_voids.stl");

    // Build an adjancecy list
    vector<vector<unsigned int>> adjList = buildAdjList(mesh);

    // Find a surface triangle
    unsigned int surfaceTriangle = findExteriorTriangle(mesh);

    // Find the triangles associated with cavities  
    vector<unsigned int> cavityTriangles = findCavityTriangles(adjList, surfaceTriangle, mesh.num_tris());

    // Identify individual bodies in the cavities in case we have multiple bodies
    vector<vector<unsigned int>> bodies = detectDisconnectedBodiesInSubset(adjList, cavityTriangles);

    // Write the cavities to an STL file (i.e. bodies that are not connected to the exterior surface)
    int bodyCounter = 1;
    for(const vector<unsigned int>& body: bodies){
        string filename = "cavity_";
        filename += std::to_string(bodyCounter);
        filename.append(".stl");
        writeSTL(filename, body, mesh);
        bodyCounter++;
    }
    
    // Orient the mesh based on given seed triangle
    cout << "Re-orienting the triangles based on seed triangle ID: " << surfaceTriangle << "\n";
    vector<unsigned int> reorientedTriangles =  orientMesh(surfaceTriangle, adjList, mesh);

    // write the oriented triangles to an STL file
    writeSTL("reoriented.stl", reorientedTriangles, mesh);

    return 0;
}

