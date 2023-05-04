#include "MaterialSystem/ElmtVarInfo.h"
const map<ElementVariableType,int> ElmtVarInfo::elmtVarCpntNum={
    {ElementVariableType::NONE,0},
    {ElementVariableType::VONMISES,1},
    {ElementVariableType::CAUCHYSTRESS,6},
    {ElementVariableType::LOGSTRAIN,6},
    {ElementVariableType::PRESSURE,1},
    {ElementVariableType::KIRCHOFFSTRESS,6},
    {ElementVariableType::JACOBIAN,1}
};
const map<ElementVariableType,vector<string>>  ElmtVarInfo::elmtVarCpntName={
    {ElementVariableType::NONE,{}},
    {ElementVariableType::VONMISES,{"mises"}},
    {ElementVariableType::CAUCHYSTRESS,{"S11","S22","S33","S12","S13","S23"}},
    {ElementVariableType::LOGSTRAIN,{"LE11","LE22","LE33","LE12","LE13","LE23"}},
    {ElementVariableType::PRESSURE,{"P"}},
    {ElementVariableType::KIRCHOFFSTRESS,{"T11","T22","T33","T12","T13","T23"}},
    {ElementVariableType::JACOBIAN,{"J"}} 
};