#include "MaterialSystem/LinearElasticMat2D.h"
#include "Utils/MessagePrinter.h"
#include "MaterialSystem/ElasticConst.h"
#include "MathUtils/ViogtRank2Tensor2D.h"
#include <cmath>
LinearElasticMat2D::LinearElasticMat2D(nlohmann::json *t_propPtr):m_strain(0.0){
    initProperty(t_propPtr);
}
void LinearElasticMat2D::initProperty(nlohmann::json *t_propPtr){
    if(t_propPtr->contains("lame")&&t_propPtr->contains("G")){
        if(t_propPtr->at("lame").is_number_float()&&t_propPtr->at("G").is_number_float()){
            m_lame=t_propPtr->at("lame"); m_G=t_propPtr->at("G");
        }
        else{
            MessagePrinter::printErrorTxt("properties K or G is not D float-point number");
            MessagePrinter::exitcfem();
        }
    }
    else if(t_propPtr->contains("E")&&t_propPtr->contains("nu")){
        if(t_propPtr->at("E").is_number_float()&&t_propPtr->at("nu").is_number_float()){
            ElasticConst::getLame_GByE_Nu(t_propPtr->at("E"),t_propPtr->at("nu"),&m_lame,&m_G);
        }
        else{
            MessagePrinter::printErrorTxt("properties E or nu is not D float-point number");
            MessagePrinter::exitcfem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("properties are not paired");
        MessagePrinter::exitcfem();        
    }
}
double exp2x(double x){return exp(2*x);}
void LinearElasticMat2D::updateMaterialBydudx(void *t_incStrainPtr,bool *t_converged){
    Rank2Tensor2d * incStrainPtr=(Rank2Tensor2d *)t_incStrainPtr;

    if(m_nLarge){
        // cal last converged B
        ViogtRank2Tensor2D B=m_strain0.isotropicFunc(&exp2x);
        // update B
        B=*incStrainPtr*B*incStrainPtr->transpose();
        // cal current strain e=0.5*lnB
        m_strain=B.isotropicFunc(&log)*0.5;
        m_F=(*incStrainPtr)*m_F0;
        m_J=m_F.det();
    }
    else{
        // cal strain_inc
        ViogtRank2Tensor2D strain_inc=(*incStrainPtr+incStrainPtr->transpose())*0.5;
        m_strain=m_strain0+strain_inc;
    }
    m_S=(m_strain*(2*m_G)+TensorConst2D::I*(m_lame*m_strain.trace()))/m_J;
    *t_converged=true;
}
void LinearElasticMat2D::getTangentModulus(void *t_incStrainPtr,void *t_D){
    // D_ijkl=lame*del_ij*del_kl+G*(del_ik*del_jl+del_il*del_jk)
    Rank2Tensor2d * incStrainPtr=(Rank2Tensor2d *)t_incStrainPtr;
    ViogtRank4Tensor2D *D=(ViogtRank4Tensor2D *)t_D;
    if(incStrainPtr){}
    if(!D){
        MessagePrinter::printErrorTxt("tangent modulus ptr D need preallocation");
        MessagePrinter::exitcfem();
    }
    *D=TensorConst2D::IXI*m_lame+TensorConst2D::IISym*(2*m_G);
}
double oneDivideX(double x){return 1.0/x;}
void LinearElasticMat2D::getSpatialTangentModulus(void *t_incStrainPtr,MatrixXd *t_a){
    // D_ijkl=lame*del_ij*del_kl+G*(del_ik*del_jl+del_il*del_jk)
    // a_ijkl=(1/2J)*[D:L:BMatrix]_ijkl-S_il*del_jk
    Rank2Tensor2d * incStrainPtr=(Rank2Tensor2d *)t_incStrainPtr;
    if(incStrainPtr){}
    if(!t_a){
        MessagePrinter::printErrorTxt("spatial tangent modulus ptr D need preallocation");
        MessagePrinter::exitcfem();
    }
    /****************************          e trial */
    ViogtRank4Tensor2D D;   //D = dtau / dE        */    
    D=TensorConst2D::IXI*m_lame+TensorConst2D::IISym*(2*m_G);
    // cal last converged B
    ViogtRank2Tensor2D B=m_strain0.isotropicFunc(&exp2x);
    // update B
    B=*incStrainPtr*B*incStrainPtr->transpose();
    ViogtRank4Tensor2D L=B.iostropicFuncDeriv(&log,&oneDivideX);
    MatrixXd BMatrix =TensorConst2D::I.ikjl(B)+B.iljk(TensorConst2D::I);
    *t_a=((D*L).toFullMatrix()*BMatrix)/(2*m_J)-m_S.iljk(TensorConst2D::I);
}
void LinearElasticMat2D::getElementVariable(ElementVariableType elmtVarType,void *elmtVarPtr){
    switch (elmtVarType)
    {
    case ElementVariableType::CAUCHYSTRESS:
        *(ViogtRank2Tensor2D *)elmtVarPtr=m_S;
        break;
    case ElementVariableType::JACOBIAN:
        *(double *)elmtVarPtr=m_J;
        break;
    case ElementVariableType::KIRCHOFFSTRESS:
        *(ViogtRank2Tensor2D *)elmtVarPtr=m_S/m_J;
        break;
    case ElementVariableType::LOGSTRAIN:
        *(ViogtRank2Tensor2D *)elmtVarPtr=m_strain;
        break;
    case ElementVariableType::PRESSURE:{
        double S33=m_lame*m_strain.trace();
        *(double *)elmtVarPtr=(m_S.trace()+S33)/(-3.0);
        break;
    }
    case ElementVariableType::VONMISES:{
        double S33=m_lame*m_strain.trace();
        double Sm=(m_S.trace()+S33)/3.0;
        *(double *)elmtVarPtr=sqrt(1.5*((m_S(0)-Sm)*(m_S(0)-Sm)+(m_S(1)-Sm)*(m_S(1)-Sm)+(S33-Sm)*(S33-Sm)+2*(m_S(2)-Sm)*(m_S(2)-Sm)));
        break;
    }
    default:
        MessagePrinter::printErrorTxt("Required material variale is not suuported by linearElastic material lib now");
        MessagePrinter::exitcfem();
        break;
    }
}