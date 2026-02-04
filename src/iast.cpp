//
// Created by alex.
//

#include "../include/iast.h"
#include<cmath>
#include <ostream>
#include "../include/nodeVisitor.h"

UCRotationNode::UCRotationNode(const std::vector<double>* angles, const bool _first, const bool _reverse, const bool _inverse) {
    num_qubits = static_cast<int>(log2(static_cast<double>(angles->size()))) + 1;
    name = index++;
    this->angles = *angles;
    first = _first;
    reverse_gate = _reverse;
    inverse = _inverse;
}

void UCRotationNode::init() {
    if (num_qubits > 2) {
        std::vector<double> angles_left(angles.size() / 2);
        std::vector<double> angles_right(angles.size() / 2);

        for (size_t i = 0; i < angles.size() / 2; i++) {
            angles_left[i] = (angles[i] + angles[i + angles.size() / 2]) / 2;
            angles_right[i] = (angles[i] - angles[i + angles.size() / 2]) / 2;
        }
        gate1 = createSubNode(&angles_left, false);
        gate2 = createSubNode(&angles_right, true);
    } else {
        gate1 = createBaseRotation((angles[0] + angles[1]) / 2.0);
        gate2 = createBaseRotation((angles[0] - angles[1]) / 2.0);
    }
}

void UCRotationNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type UCRotationNode::get_data() {
    return angles;
}
