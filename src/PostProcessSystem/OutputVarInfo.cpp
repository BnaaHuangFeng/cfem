#include "PostProcessSystem/OutputVarInfo.h"
const map<FieldVariableType,VarPosition> FieldVarInfo::varPosition={
    {FieldVariableType::VONMISES,  VarPosition::ELEMENT},
    {FieldVariableType::STRESS,    VarPosition::ELEMENT},
    {FieldVariableType::LOGSTRAIN, VarPosition::ELEMENT},
    {FieldVariableType::U,         VarPosition::NODE}
};

const map<FieldVariableType,VarMathType> FieldVarInfo::varMathType={
    {FieldVariableType::VONMISES,  VarMathType::SCALAR},
    {FieldVariableType::STRESS,    VarMathType::TENSORRANK2},
    {FieldVariableType::LOGSTRAIN, VarMathType::TENSORRANK2},
    {FieldVariableType::U,         VarMathType::VECTOR}
};

const map<FieldVariableType,int> FieldVarInfo::varType={
    {FieldVariableType::VONMISES,  int(ElementVariableType::VONMISES)},
    {FieldVariableType::STRESS,    int(ElementVariableType::CAUCHYSTRESS)},
    {FieldVariableType::LOGSTRAIN, int(ElementVariableType::LOGSTRAIN)},
    {FieldVariableType::U,         int(NodeVariableType::U)}
};

const map<HistoryVariableType,string> HisVarInfo::varName={
    {HistoryVariableType::U,        "U"},
    {HistoryVariableType::U1,       "U1"},
    {HistoryVariableType::U2,       "U2"},
    {HistoryVariableType::U3,       "U3"},
    {HistoryVariableType::RF,       "RF"},
    {HistoryVariableType::RF1,      "RF1"},
    {HistoryVariableType::RF2,      "RF2"},
    {HistoryVariableType::RF3,      "RF3"}
};

const map<HistoryVariableType,VarPosition> HisVarInfo::varPosition={
    {HistoryVariableType::U,        VarPosition::NODE},
    {HistoryVariableType::U1,       VarPosition::NODE},
    {HistoryVariableType::U2,       VarPosition::NODE},
    {HistoryVariableType::U3,       VarPosition::NODE},
    {HistoryVariableType::RF,       VarPosition::NODE},
    {HistoryVariableType::RF1,      VarPosition::NODE},
    {HistoryVariableType::RF2,      VarPosition::NODE},
    {HistoryVariableType::RF3,      VarPosition::NODE}
};

const map<HistoryVariableType,VarMathType> HisVarInfo::varMathType={
    {HistoryVariableType::U,        VarMathType::VECTOR},
    {HistoryVariableType::U1,       VarMathType::SCALAR},
    {HistoryVariableType::U2,       VarMathType::SCALAR},
    {HistoryVariableType::U3,       VarMathType::SCALAR},
    {HistoryVariableType::RF,       VarMathType::VECTOR},
    {HistoryVariableType::RF1,      VarMathType::SCALAR},
    {HistoryVariableType::RF2,      VarMathType::SCALAR},
    {HistoryVariableType::RF3,      VarMathType::SCALAR}
};

const map<HistoryVariableType,int> HisVarInfo::varType={
    {HistoryVariableType::U,        int(NodeVariableType::U)},
    {HistoryVariableType::U1,       int(NodeVariableType::U)},
    {HistoryVariableType::U2,       int(NodeVariableType::U)},
    {HistoryVariableType::U3,       int(NodeVariableType::U)},
    {HistoryVariableType::RF,       int(NodeVariableType::RF)},
    {HistoryVariableType::RF1,      int(NodeVariableType::RF)},
    {HistoryVariableType::RF2,      int(NodeVariableType::RF)},
    {HistoryVariableType::RF3,      int(NodeVariableType::RF)}
};

const map<HistoryVariableType,int> HisVarInfo::cpntInd={
    {HistoryVariableType::U,        -1},
    {HistoryVariableType::U1,       0},
    {HistoryVariableType::U2,       1},
    {HistoryVariableType::U3,       2},
    {HistoryVariableType::RF,       -1},
    {HistoryVariableType::RF1,      0},
    {HistoryVariableType::RF2,      1},
    {HistoryVariableType::RF3,      2}
};
const map<HistoryVariableType,int> HisVarInfo::cpntNum={
    {HistoryVariableType::U,        3},
    {HistoryVariableType::U1,       1},
    {HistoryVariableType::U2,       1},
    {HistoryVariableType::U3,       1},
    {HistoryVariableType::RF,       3},
    {HistoryVariableType::RF1,      1},
    {HistoryVariableType::RF2,      1},
    {HistoryVariableType::RF3,      1}
};
const map<HistoryVariableType,VarOutputForm> HisVarInfo::varOutputForm={
    {HistoryVariableType::U,        VarOutputForm::ANY},
    {HistoryVariableType::U1,       VarOutputForm::ANY},
    {HistoryVariableType::U2,       VarOutputForm::ANY},
    {HistoryVariableType::U3,       VarOutputForm::ANY},
    {HistoryVariableType::RF,       VarOutputForm::ANY},
    {HistoryVariableType::RF1,      VarOutputForm::SUM},
    {HistoryVariableType::RF2,      VarOutputForm::SUM},
    {HistoryVariableType::RF3,      VarOutputForm::SUM}
};