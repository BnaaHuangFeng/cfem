#pragma once
#include "InputSystem/EnumDataType.h"
#include <map>
#include <vector>
using namespace std;
namespace NodeVarInfo{
    extern const map<NodeVariableType,int>             nodeVarCpntNum;   /**< component num of elemental variable*/
    extern const map<NodeVariableType,vector<string>>  nodeVarCpntName;
    extern const map<NodeVariableType,string>          nodeVarName;
}