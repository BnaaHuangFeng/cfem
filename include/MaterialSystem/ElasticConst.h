#pragma once
namespace ElasticConst{
    void getLame_GByE_Nu(double E,double nu,double *lame,double *G);
    double getLameByE_Nu(double E,double nu);
    double getKByE_Nu(double E,double nu);
    double getGByE_Nu(double E,double nu);
}