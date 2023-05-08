#include "MeshSystem/NodeVarInfo.h"
const map<NodeVariableType,int> NodeVarInfo::nodeVarCpntNum={
    {NodeVariableType::NONE,0},
    {NodeVariableType::COORD,3},
    {NodeVariableType::U,3},
    {NodeVariableType::UINC,3},
    {NodeVariableType::RF,3},
    {NodeVariableType::RESIDUAL,3},
    {NodeVariableType::LOAD,3}
};
const map<NodeVariableType,vector<string>>  NodeVarInfo::nodeVarCpntName={
    {NodeVariableType::NONE,    {}},
    {NodeVariableType::COORD,   {"X1","X2","X3"}},
    {NodeVariableType::U,       {"U1","U2","U3"}},
    {NodeVariableType::UINC,    {"DU1","DU2","DU3"}},
    {NodeVariableType::RF,      {"RF1","RF2","RF3"}},
    {NodeVariableType::RESIDUAL,{"RES1","RES2","RES3"}},
    {NodeVariableType::LOAD,    {"F1","F2","F3"}} 
};

const map<NodeVariableType,string>  NodeVarInfo::nodeVarName={
    {NodeVariableType::NONE,    ""},
    {NodeVariableType::COORD,   "X"},
    {NodeVariableType::U,       "U"},
    {NodeVariableType::UINC,    "DU"},
    {NodeVariableType::RF,      "RF"},
    {NodeVariableType::RESIDUAL,"RES"},
    {NodeVariableType::LOAD,    "F"} 
};