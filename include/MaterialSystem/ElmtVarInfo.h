#include "InputSystem/EnumDataType.h"
#include <map>
#include <vector>
using namespace std;
namespace ElmtVarInfo{
    extern const map<ElementVariableType,int>             elmtVarCpntNum;   /**< component num of elemental variable*/
    extern const map<ElementVariableType,vector<string>>  elmtVarCpntName;
}