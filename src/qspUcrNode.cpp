#include "../include/qspUcrNode.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include "../include/nodeVisitor.h"

qspUcrNode::qspUcrNode(const std::vector<std::complex<double>>* state) {
    num_qubits = static_cast<int>(log2(static_cast<double>(state->size())));
    name = index++;
    this->state = *state;
    first = false;
    reverse_gate = false;
    this->next_state.resize(state->size()/2);
    this->angles_ry.resize(state->size()/2);
    this->angles_rz.resize(state->size()/2);
    init();
}

void qspUcrNode::init(){
  
    if(num_qubits > 1){
        for (int i = 0; i < this->state.size()/2; i++){

            this->angles_rz[i] = std::remainder(std::arg(this->state[2*i+1]) - std::arg(this->state[2*i]), 2*std::numbers::pi);

            double norm = std::sqrt((std::pow(std::abs(this->state[2*i]),2) + std::pow(std::abs(this->state[2*i+1]),2)));
            double phase = (std::arg(this->state[2*i]) + std::arg(this->state[2*i+1])) / 2.0;
            this->next_state[i] = std::polar(norm, phase);

            if (norm < 1e-12){
                this->angles_ry[i] = 0.0;
            } else {
                this->angles_ry[i] = 2 * std::asin(std::abs(this->state[2*i+1]) / std::abs(this->next_state[i]));
            }

            
        }

        next_qsp = createSubQsp(&next_state);
        ucry = createUcry(&angles_ry);
        ucrz = createUcrz(&angles_rz);  
    }
    else {

        this->global_phase = (std::arg(this->state[0]) + std::arg(this->state[1])) / 2.0;

        double base_norm = std::sqrt((std::pow(std::abs(this->state[0]),2) + std::pow(std::abs(this->state[1]),2)));
        double angle_rz = std::remainder(std::arg(this->state[1]) - std::arg(this->state[0]), 2*std::numbers::pi);
        double angle_ry = 0.0;
        if (base_norm > 1e-12) {
            angle_ry = std::abs(this->state[1]) / base_norm;
        }

        this->ucry = createBaseRy(2 * std::asin(angle_ry));
        this->ucrz = createBaseRz(angle_rz);
    }

}

void qspUcrNode::accept(nodeVisitor &visitor) {
    visitor.visit(*this);
}

return_type qspUcrNode::get_data() {
    return state;
}

std::unique_ptr<ryNode> qspUcrNode::createBaseRy(double angle_ry) {
    return std::make_unique<ryNode>(angle_ry);
}

std::unique_ptr<rzNode> qspUcrNode::createBaseRz(double angle_rz) {
    return std::make_unique<rzNode>(angle_rz);
}

std::unique_ptr<qspUcrNode> qspUcrNode::createSubQsp(const std::vector<std::complex<double>>* subState) {
    return std::make_unique<qspUcrNode>(subState);
}

std::unique_ptr<ucryNode> qspUcrNode::createUcry(const std::vector<double>* angles_ry) {
    return std::make_unique<firstUcryNode>(angles_ry, true);
}

std::unique_ptr<ucrzNode> qspUcrNode::createUcrz(const std::vector<double>* angles_rz) {
    return std::make_unique<firstUcrzNode>(angles_rz, true);
}