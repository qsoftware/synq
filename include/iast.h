//
// Created by alex.
//

#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <variant>
#include <complex>

class nodeVisitor;

// Adicione std::vector<std::complex<double>> na lista
using return_type = std::variant<double, std::vector<double>, std::vector<std::complex<double>>, Eigen::Matrix2cf>;

/**
 * IASTnode class is the base class of the Abstract Syntax Tree
 * used to parse quantum operations into quantum assembler.
 */
class IASTnode
{
public:
    inline static int index = 0;
    int name;
    int num_qubits = 0;
    bool reverse_gate = false;
    virtual ~IASTnode() = default;
    virtual void accept(nodeVisitor &visitor) = 0;
    virtual return_type get_data() = 0;
    virtual int get_num_qubits();
private:
    std::string mName;
};

class uGate:public IASTnode {

};
