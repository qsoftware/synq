#include "../include/hNode.h"
#include "../include/nodeVisitor.h"

hNode::hNode(const int position) : position(position) {
    num_qubits = 1;
}

void hNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type hNode::get_data() {
    return std::monostate{};
}

std::unique_ptr<IASTnode> hNode::createHNode(const int position) {
    return std::make_unique<hNode>(position);
}
