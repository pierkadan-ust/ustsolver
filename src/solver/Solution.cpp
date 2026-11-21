#include "solver/Solution.hpp"

void Solution::write_tecplot(Variable* variable, string path, label if_print, string type) {

    if (if_print) {
        std::cout << "\nSolution    : ";
    }

    if (variable->ZONE->type_mesh == 1) {  // 1 - Quad
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"" << variable->NAME << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([3]=cellcentered)" << std::endl;
        outfile << "zonetype=fequadrilateral" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        if (type == "magnitude") {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << sqrt(pow(variable->X(i, 0), 2) + pow(variable->X(i, 1), 2)) << std::endl; }
        }
        else {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i) << std::endl; }
        }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << " "
                << variable->ZONE->cells(i, 3) + 1 << std::endl;
        }
        outfile.close();
    }

    else if (variable->ZONE->type_mesh == 2) {  // 2 - Tri
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"" << variable->NAME << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([3]=cellcentered)" << std::endl;
        outfile << "zonetype=fetriangle" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        if (type == "magnitude") {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << sqrt(pow(variable->X(i, 0), 2) + pow(variable->X(i, 1), 2)) << std::endl; }
        }
        else {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i) << std::endl; }
        }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << std::endl;
        }
        outfile.close();
    }

    else if (variable->ZONE->type_mesh == 3) {  // 3 - Hex
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"z\",\"" << variable->NAME << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([4]=cellcentered)" << std::endl;
        outfile << "zonetype=febrick" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 2) << std::endl; }
        if (type == "magnitude") {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << sqrt(pow(variable->X(i, 0), 2) + pow(variable->X(i, 1), 2) + pow(variable->X(i, 2), 2)) << std::endl; }
        }
        else {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i) << std::endl; }
        }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << " "
                << variable->ZONE->cells(i, 3) + 1 << " "
                << variable->ZONE->cells(i, 4) + 1 << " "
                << variable->ZONE->cells(i, 5) + 1 << " "
                << variable->ZONE->cells(i, 6) + 1 << " "
                << variable->ZONE->cells(i, 7) + 1 << std::endl;
        }
        outfile.close();
    }

    else if (variable->ZONE->type_mesh == 4) {  // 4 - Tetra
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"z\",\"" << variable->NAME << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([4]=cellcentered)" << std::endl;
        outfile << "zonetype=fetetrahedron" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 2) << std::endl; }
        if (type == "magnitude") {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << sqrt(pow(variable->X(i, 0), 2) + pow(variable->X(i, 1), 2) + pow(variable->X(i, 2), 2)) << std::endl; }
        }
        else {
            for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i) << std::endl; }
        }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << " "
                << variable->ZONE->cells(i, 3) + 1 << std::endl;
        }
        outfile.close();
    }

    if (if_print) {
        std::cout << "The solution of [" << variable->NAME << "] is saved at \"" << path << "\"\n";
    }
}

void Solution::write_tecplot(Variable* variable, label dim, string name, string path, label if_print) {

    if (if_print) {
        std::cout << "\nSolution    : ";
    }

    if (variable->ZONE->type_mesh == 1) {  // 1 - Quad
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"" << name << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([3]=cellcentered)" << std::endl;
        outfile << "zonetype=fequadrilateral" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i, dim) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << " "
                << variable->ZONE->cells(i, 3) + 1 << std::endl;
        }
        outfile.close();
    }

    else if (variable->ZONE->type_mesh == 2) {  // 2 - Tri
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"" << name << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([3]=cellcentered)" << std::endl;
        outfile << "zonetype=fetriangle" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i, dim) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << std::endl;
        }
        outfile.close();
    }

    else if (variable->ZONE->type_mesh == 3) {  // 3 - Hex
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"z\",\"" << name << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([4]=cellcentered)" << std::endl;
        outfile << "zonetype=febrick" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 2) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i, dim) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << " "
                << variable->ZONE->cells(i, 3) + 1 << " "
                << variable->ZONE->cells(i, 4) + 1 << " "
                << variable->ZONE->cells(i, 5) + 1 << " "
                << variable->ZONE->cells(i, 6) + 1 << " "
                << variable->ZONE->cells(i, 7) + 1 << std::endl;
        }
        outfile.close();
    }

    else if (variable->ZONE->type_mesh == 4) {  // 4 - Tetra
        std::ofstream outfile(path);
        outfile << "variables=\"x\",\"y\",\"z\",\"" << name << "\"" << std::endl;
        outfile << "zone N=" << variable->ZONE->num_nodes << ", E=" << variable->ZONE->num_cells << ", datapacking=block" << std::endl;
        outfile << "varlocation=([4]=cellcentered)" << std::endl;
        outfile << "zonetype=fetetrahedron" << std::endl;
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 0) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 1) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_nodes; i++) { outfile << variable->ZONE->nodes(i, 2) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) { outfile << variable->X(i, dim) << std::endl; }
        for (label i = 0; i < variable->ZONE->num_cells; i++) {
            outfile << variable->ZONE->cells(i, 0) + 1 << " "
                << variable->ZONE->cells(i, 1) + 1 << " "
                << variable->ZONE->cells(i, 2) + 1 << " "
                << variable->ZONE->cells(i, 3) + 1 << std::endl;
        }
        outfile.close();
    }

    if (if_print) {
        std::cout << "The solution of [" << variable->NAME << "(" << name << ")" << "] is saved at \"" << path << "\"\n";
    }
}
