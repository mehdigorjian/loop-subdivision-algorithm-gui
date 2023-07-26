#include "Viewer.h"
#include "loop_subdiv.h"
#include <cstddef>
#include <igl/readOBJ.h>
#include <igl/readSTL.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/writeOBJ.h>
#include <string>

int main(int argc, char *argv[])
{
    // std::string stl_file_path = "../data/hw1.stl";
    std::cout << "Enter '.stl' file path: ";
    std::string stl_file_path;
    std::cin >> stl_file_path;

    std::size_t file_split = stl_file_path.find_last_of("/");
    std::string path = stl_file_path.substr(0, file_split);
    std::string file = stl_file_path.substr(file_split + 1);

    std::size_t ext_split = file.find_last_of(".");
    std::string ext_name = file.substr(0, ext_split);
    // std::cout << file << "---" << path << "---" << ext_name << std::endl;
    std::string obj_file_path = path + "/obj_" + ext_name + ".obj";
    // std::cout << obj_file_path << std::endl;
    std::string subdiv_obj_file_path = path + "/subdiv_obj_" + ext_name + ".obj";

    Eigen::MatrixXd V, temp_V;
    Eigen::MatrixXi F, N, SVI, SVJ;

    igl::readSTL(stl_file_path, temp_V, F, N);
    igl::remove_duplicate_vertices(temp_V, F, 0, V, SVI, SVJ, F);
    igl::writeOBJ(obj_file_path, V, F);
    igl::readOBJ(obj_file_path, V, F);

    if (F.cols() != 3 || F.minCoeff() < 0)
    {
        std::cerr << "ERROR: NOT A TRINAGULAR MESH!" << std::endl;
        return EXIT_FAILURE;
    }
    // reset the subdivision
    Eigen::MatrixXd OV = V;
    Eigen::MatrixXi OF = F;
    bool show_lines = true;

    Viewer v;
    std::cout << R"(GUI shortcut:
[spacebar]    run subdivision
  2           run subdivision 2 times
  R,r         reset mesh
  S,s         save subdivided obj

)";
    v.set_mesh(V, F);
    v.callback_key_pressed = [&v, &V, &F, &OV, &OF, &show_lines, &subdiv_obj_file_path](igl::opengl::glfw::Viewer &, unsigned int key, int) -> bool
    {
        switch (key)
        {
        default:
            return false;
        case 'R':
        case 'r':
            V = OV;
            F = OF;
            v.set_mesh(V, F);
            break;
        case '2':
            loop_subdivision(Eigen::MatrixXd(V), Eigen::MatrixXi(F), 2, V, F);
            v.set_mesh(V, F);
            break;
        case ' ':
            loop_subdivision(Eigen::MatrixXd(V), Eigen::MatrixXi(F), 1, V, F);
            v.set_mesh(V, F);
            break;
        case 'S':
        case 's':
            igl::writeOBJ(subdiv_obj_file_path, V, F);
            break;
        }
        return true;
    };
    v.launch();
}
