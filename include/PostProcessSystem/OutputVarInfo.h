#pragma once
#include "InputSystem/EnumDataType.h"
#include <map>
#include <vector>
using namespace std;
enum class VarPosition{
    NODE,
    ELEMENT
};
enum class VarMathType{
    SCALAR,
    VECTOR,
    TENSORRANK2,
    TENSORRANK4
};
enum class VarOutputForm{
    ANY,
    ALL,
    SUM
};
namespace FieldVarInfo{
    extern const map<FieldVariableType,VarPosition>     varPosition;    /**< position where variable is defined*/
    extern const map<FieldVariableType,VarMathType>     varMathType;    /**< varible's form in math*/
    extern const map<FieldVariableType,int>             varType;        /**< corresponding NodeVariableType/ElementVariableType*/
};
namespace HisVarInfo{
    extern const map<HistoryVariableType,string>        varName;
    extern const map<HistoryVariableType,VarPosition>   varPosition;    /**< component num of elemental variable*/
    extern const map<HistoryVariableType,VarMathType>   varMathType;    /**< varible's form in math*/
    extern const map<HistoryVariableType,int>           varType;        /**< corresponding NodeVariableType/ElementVariableType*/  
    extern const map<HistoryVariableType,int>           cpntInd;        /**< corresponding component index*/ 
    extern const map<HistoryVariableType,int>           cpntNum;        /**< variable's component num*/
    extern const map<HistoryVariableType,VarOutputForm> varOutputForm;  /**< variable's output form*/ 
};