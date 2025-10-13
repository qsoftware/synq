//
// Created by adenilton on 10/13/25.
//

#pragma once

#include <variant>
#include <vector>
#include<string>

class nodeVisitor;

using return_type = std::variant<double, std::vector<double> >;

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


