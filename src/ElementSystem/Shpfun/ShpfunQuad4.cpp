#include "ElementSystem/Shpfun/ShpfunQuad4.h"
/**
 * COMPUTES SHAPE FUNCTIONS AND SHAPE FUNCTION DERIVATIVES FOR
 * ELEMENT 'QUAD_4':
 *                         4         3
 *                          o-------o
 *                          |       |     STANDARD ISOPARAMETRIC
 *                          |       |     BI-LINEAR 4-NODE QUADRILATERAL 
 *                          |       |
 *                          o-------o
 *                         1         2
 *
 * REFERENCE: Expression (4.42)
*/
const double ShpfunQuad4::m_hgModal[4]={1.0, -1.0, 1.0, -1.0};
ShpfunQuad4::ShpfunQuad4(){
    m_x0[0](0)=-1.0;
    m_x0[0](1)=-1.0;

    m_x0[1](0)=1.0;
    m_x0[1](1)=-1.0;

    m_x0[2](0)=1.0;
    m_x0[2](1)=1.0;

    m_x0[3](0)=-1.0;
    m_x0[3](1)=1.0;
}
ShpfunQuad4::ShpfunQuad4(Vector2d t_r, Vector2d t_x0[]){
    Shpfun2D::setNatCoords(t_r);
    for(int i=0;i<m_funs;i++){
        m_x0[i]=t_x0[i];
    }
}
void ShpfunQuad4::init(Vector2d t_r,Vector2d t_x0[]){
    Shpfun2D::setNatCoords(t_r);
    ShpfunQuad4::setRefCoords(t_x0);
}
void ShpfunQuad4::setRefCoords(Vector2d t_x0[]){
    for(int i=0;i<m_funs;i++){
        m_x0[i]=t_x0[i];
    }   
}
void ShpfunQuad4::getShpfunVal(double t_N[]){
    // Shape functions and derivatives on element DOMAIN
    // -------------------------------------------------
    double S=m_r(0);
    double T=m_r(1);
    double ST=S*T;
    // Shape functions
    t_N[0]=(1.0-T-S+ST)*0.25;
    t_N[1]=(1.0-T+S-ST)*0.25;
    t_N[2]=(1.0+T+S+ST)*0.25;
    t_N[3]=(1.0+T-S-ST)*0.25;
}
void ShpfunQuad4::getCoords0(Vector2d t_x0[]){
    for(int i=0;i<m_funs;i++){
        t_x0[i]=m_x0[i];
    }
}
void ShpfunQuad4::getDer2Nat(Vector2d t_dNdr[]){
    // Shape function derivatives
    double xi=m_r(0);
    double eta=m_r(1);
    t_dNdr[0](0)=(eta-1.0)/4.0;
    t_dNdr[0](1)=(xi -1.0)/4.0;

    t_dNdr[1](0)= (1.0-eta)/4.0;
    t_dNdr[1](1)=-(1.0+xi )/4.0;

    t_dNdr[2](0)= (1.0+eta)/4.0;
    t_dNdr[2](1)= (1.0+xi )/4.0;

    t_dNdr[3](0)=-(1.0+eta)/4.0;
    t_dNdr[3](1)= (1.0-xi )/4.0;

}
void ShpfunQuad4::getDer2Ref(Vector2d t_dNdx0[]){
    getDer2Nat(m_dNdr);
    Rank2Tensor2d dx0dr;
    for(int i=0;i<m_funs;i++){
        for(int rowI=0;rowI<m_dim;++rowI){
            for(int colI=0;colI<m_dim;++colI){
                dx0dr(rowI,colI)+=m_dNdr[i](colI)*m_x0[i](rowI);                   
            }
        }
    }
    Rank2Tensor2d drdx0=dx0dr.inverse();
    for(int i=0;i<m_funs;i++){
        t_dNdx0[i]=m_dNdr[i]*drdx0;
    }
}
void ShpfunQuad4::getHGShpVec(Vector2d t_dNdx2[4],Vector2d t_x2[],double t_gamma[]){
    for(int nI=0;nI<m_funs;++nI){
        t_gamma[nI]=m_hgModal[nI];
        for(int nJ=0;nJ<m_funs;++nJ){
            for(int di=0;di<m_dim;++di){
                t_gamma[nI]-=t_dNdx2[nI](di)*t_x2[nJ](di)*m_hgModal[nJ];
            }
        }
    }
}