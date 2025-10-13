//
// Created by adenilton on 10/4/25.
//

#include "../include/ucrzNode.h"

#include <algorithm>
#include<cmath>
#include <iostream>
#include <ostream>
#include "../include/nodeVisitor.h"

ucrzNode::ucrzNode(const std::vector<double>* angles, const bool _first = false, const bool _reverse=false){
    num_qubits = static_cast<int>(log2(static_cast<double>(angles->size()))) + 1;
    name = index++;
    const int half_size = static_cast<int>(angles->size()) / 2;
    first = _first;
    reverse_gate = _reverse;

    if (num_qubits > 2) {
        std::vector<double> angles_left(half_size);
        std::vector<double> angles_right(half_size);

        for (int i = 0; i < angles->size() / 2; i++) {
            angles_left[i] = (angles->at(i) + angles->at(i+half_size)) / 2;
            angles_right[i] = (angles->at(i) - angles->at(i+half_size)) / 2;
        }
        gate1 = std::make_unique<ucrzNode>(&angles_left, false, false);
        gate2 = std::make_unique<ucrzNode>(&angles_right, false, true);
    }
    else {
        gate1 = std::make_unique<rzNode>((angles->at(0) + angles->at(1)) / 2.0);
        gate2 = std::make_unique<rzNode>((angles->at(0) - angles->at(1)) / 2.0);
    }
}

void ucrzNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type ucrzNode::get_data() {
    return angles;
}

firstUcrzNode::firstUcrzNode(const std::vector<double>* angles) : ucrzNode(angles, true){
}

void firstUcrzNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}



rzNode::rzNode(const double theta) {
    num_qubits = 1;
    angle = theta;
}

rzNode::rzNode() {
    num_qubits = 1;
    angle = 0.0;
}

void rzNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type rzNode::get_data() {
    return angle;
}
