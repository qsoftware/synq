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
            // The average phase must be consistent with the *wrapped* rz angle:
            // Rz(phi) maps it to arg(a) and arg(b) only if phase = arg(a) + phi/2.
            // Using (arg(a) + arg(b))/2 is off by pi whenever arg(b)-arg(a) is wrapped.
            double phase = std::arg(this->state[2*i]) + this->angles_rz[i] / 2.0;
            this->next_state[i] = std::polar(norm, phase);

            if (norm < 1e-12){
                this->angles_ry[i] = 0.0;
            } else {
                // atan2 is equivalent to 2*asin(|b|/norm) but never returns NaN due to rounding.
                this->angles_ry[i] = 2 * std::atan2(std::abs(this->state[2*i+1]), std::abs(this->state[2*i]));
            }

            
        }

        next_qsp = createSubQsp(&next_state);
        ucry = createUcry(&angles_ry);
        ucrz = createUcrz(&angles_rz);  
    }
    else {

        double angle_rz = std::remainder(std::arg(this->state[1]) - std::arg(this->state[0]), 2*std::numbers::pi);
        // Consistent with the wrapped angle_rz (see the comment in the loop above).
        this->global_phase = std::arg(this->state[0]) + angle_rz / 2.0;

        double base_norm = std::sqrt((std::pow(std::abs(this->state[0]),2) + std::pow(std::abs(this->state[1]),2)));
        double angle_ry = 0.0;
        if (base_norm > 1e-12) {
            angle_ry = 2 * std::atan2(std::abs(this->state[1]), std::abs(this->state[0]));
        }

        this->ucry = createBaseRy(angle_ry);
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