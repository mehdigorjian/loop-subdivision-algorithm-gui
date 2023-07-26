#include "loop_subdiv.h"
#include <functional>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

struct hash_pair
{
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2> &p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

void loop_subdivision(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const int num_iters,
    Eigen::MatrixXd &SV,
    Eigen::MatrixXi &SF)
{
    if (num_iters <= 0)
        return;

    const int V_rows = V.rows();
    const int F_rows = F.rows();

    // build edge to (one or two) faces
    std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> edge_to_faces;
    int current_vertex_index, next_vertex_index, previous_vertex_index;
    std::pair<int, int> pair;

    for (int i = 0; i < F_rows; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            current_vertex_index = F(i, j);
            next_vertex_index = F(i, (j + 1) % 3);
            pair = (current_vertex_index <= next_vertex_index) ? std::make_pair(current_vertex_index, next_vertex_index) : std::make_pair(next_vertex_index, current_vertex_index);
            if (edge_to_faces.find(pair) == edge_to_faces.end())
            {
                std::vector<int> entry;
                edge_to_faces[pair] = entry;
            }
            edge_to_faces[pair].push_back(i);
        }
    }

    SV.resize(V_rows + edge_to_faces.size(), 3);
    SF.resize(4 * F_rows, 3);

    // build vertex to neightboring vertices
    std::unordered_map<int, std::vector<int>> vertex_to_vertices;
    int a, b;

    for (auto entry : edge_to_faces)
    {
        a = entry.first.first;
        b = entry.first.second;

        if (vertex_to_vertices.find(a) == vertex_to_vertices.end())
        {
            std::vector<int> entry_vec;
            vertex_to_vertices[a] = entry_vec;
        }
        vertex_to_vertices[a].push_back(b);

        if (vertex_to_vertices.find(b) == vertex_to_vertices.end())
        {
            std::vector<int> entry_vec;
            vertex_to_vertices[b] = entry_vec;
        }
        vertex_to_vertices[b].push_back(a);
    }

    std::vector<int> faces;
    Eigen::RowVector3d position, A, B, C, D;
    std::vector<int> neighbour_vertices;
    int size;
    float BETA;

    // Add (modified) old vertices into SV
    int SV_index = 0;
    for (int i = 0; i < V_rows; i++)
    {
        neighbour_vertices = vertex_to_vertices[i];
        size = neighbour_vertices.size();
        if (size == 2)
        {
            // Boundary vertex
            A = V.row(neighbour_vertices[0]);
            B = V.row(neighbour_vertices[1]);
            position = 1.0 / 8.0 * (A + B) + 3.0 / 4.0 * V.row(i);
        }
        else
        {
            // Interior vertex
            BETA = 1.0 / size * (5.0 / 8.0 - (3.0 / 8.0 + 1.0 / 4.0 * cos(2.0 * M_PI / size)) * (3.0 / 8.0 + 1.0 / 4.0 * cos(2.0 * M_PI / size)));
            position = V.row(i) * (1 - size * BETA);
            for (int w = 0; w < size; w++)
            {
                position += BETA * V.row(neighbour_vertices[w]);
            }
        }
        SV.row(SV_index++) = position;
    }

    // Add new (edge) points to SV and populate edge_to_SV
    std::unordered_map<std::pair<int, int>, int, hash_pair> edge_to_SV; // Key represents an edge. Value Represents the row number of the edge point in SV.
    for (auto entry : edge_to_faces)
    {
        pair = entry.first;
        faces = entry.second;
        if (faces.size() == 2)
        {
            A = V.row(pair.first);
            B = V.row(pair.second);
            for (int j = 0; j < 3; j++)
            {
                if (F(faces[0], j) != pair.first && F(faces[0], j) != pair.second)
                {
                    C = V.row(F(faces[0], j));
                    break;
                }
            }
            for (int j = 0; j < 3; j++)
            {
                if (F(faces[1], j) != pair.first && F(faces[1], j) != pair.second)
                {
                    D = V.row(F(faces[1], j));
                    break;
                }
            }
            position = 3.0 / 8.0 * (A + B) + 1.0 / 8.0 * (C + D);
        }
        else if (faces.size() == 1)
        {
            A = V.row(pair.first);
            B = V.row(pair.second);
            position = 1.0 / 2.0 * (A + B);
        }
        else
        {
            std::cerr << "Edge connects to " << faces.size() << " faces\n";
        }
        SV.row(SV_index) = position;
        edge_to_SV[pair] = SV_index;
        SV_index++;
    }

    // Populate SF
    for (int i = 0; i < F_rows; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            SF(i * 4 + j, 0) = F(i, j);
            pair = (F(i, j) <= F(i, (j + 1) % 3)) ? std::make_pair(F(i, j), F(i, (j + 1) % 3)) : std::make_pair(F(i, (j + 1) % 3), F(i, j));
            SF(i * 4 + j, 1) = edge_to_SV[pair];
            SF(i * 4 + 3, j) = SF(i * 4 + j, 1);
            pair = (F(i, j) <= F(i, (j + 2) % 3)) ? std::make_pair(F(i, j), F(i, (j + 2) % 3)) : std::make_pair(F(i, (j + 2) % 3), F(i, j));
            SF(i * 4 + j, 2) = edge_to_SV[pair];
        }
    }
    loop_subdivision(Eigen::MatrixXd(SV), Eigen::MatrixXi(SF), num_iters - 1, SV, SF);
}
