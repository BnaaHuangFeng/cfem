#include "MaterialSystem/ElasticConst.h"
#include "Utils/MessagePrinter.h"
void ElasticConst::getLame_GByE_Nu(double E,double nu,double *lame,double *G){
    *lame=getLameByE_Nu(E,nu);
    *G=getGByE_Nu(E,nu);
}
double ElasticConst::getLameByE_Nu(double E,double nu){
    return E*nu/((1.0+nu)*(1.0-2.0*nu));
}
double ElasticConst::getKByE_Nu(double E,double nu){
    return E/(3.0*(1-2.0*nu));
}
double ElasticConst::getGByE_Nu(double E,double nu){
    return 0.5*E/(1.0+nu);
}