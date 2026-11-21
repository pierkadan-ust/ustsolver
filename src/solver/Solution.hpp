// Write Solution

#pragma once

#include "base/Global.hpp"
#include "field/Variable.hpp"

namespace Solution {

    void write_tecplot(Variable* variable, string path, label if_print = 1, string type = "none");

    void write_tecplot(Variable* variable, label dim, string name, string path, label if_print = 1);
}
