#include "ElementSystem/Shpfun/ShpfunQuad4.h"
#include "InputSystem/EnumDataType.h"
#include "MathUtils/Rank2Tensor2d.h"
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
ShpfunQuad4::ShpfunQuad4(){
    m_x0[0](1)=-1.0;
    m_x0[0](2)=-1.0;

    m_x0[1](1)=1.0;
    m_x0[1](2)=-1.0;

    m_x0[2](1)=1.0;
    m_x0[2](2)=1.0;

    m_x0[3](1)=-1.0;
    m_x0[3](2)=1.0;
}
ShpfunQuad4::ShpfunQuad4(Vector2d t_x0[]){
    for(int i=0;i<m_funs;i++){
        m_x0[i]=t_x0[i];
    }
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
    double S=m_r(1);
    double T=m_r(2);
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
    double xi=m_r(1);
    double eta=m_r(2);
    t_dNdr[1-1](1)=(eta-1.0)/4.0;
    t_dNdr[1-1](2)=(xi -1.0)/4.0;

    t_dNdr[2-1](1)= (1.0-eta)/4.0;
    t_dNdr[2-1](2)=-(1.0+xi )/4.0;

    t_dNdr[3-1](1)= (1.0+eta)/4.0;
    t_dNdr[3-1](2)= (1.0+xi )/4.0;

    t_dNdr[4-1](1)=-(1.0+eta)/4.0;
    t_dNdr[4-1](2)= (1.0-xi )/4.0;

}
void ShpfunQuad4::getDer2Ref(Vector2d t_dNdx0[]){
    getDer2Nat(m_dNdr);
    Rank2Tensor2d m_dx0dr(0.0);
    for(int i=0;i<m_funs;i++){
        m_dx0dr.getIthRow(1)+=m_dNdr[i]*m_x0[i](1);
        m_dx0dr.getIthRow(2)+=m_dNdr[i]*m_x0[i](2);
    }
    Rank2Tensor2d drdx0=m_dx0dr.inverse();
    for(int i=0;i<m_funs;i++){
        t_dNdx0[i]=m_dNdr[i]*drdx0;
    }
}