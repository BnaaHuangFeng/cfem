#include "MaterialSystem/NeoHookeanAbq2d.h"
#include "Utils/MessagePrinter.h"
#include "MaterialSystem/ElasticConst.h"
NeoHookeanAbq2d::NeoHookeanAbq2d(bool ifLarge,double t_det_dx0dr,nlohmann::json *t_propPtr):
                        Material2D(ifLarge,t_det_dx0dr),m_ifPropInit(false),
                        m_B(ViogtRank2Tensor2D::InitMethod::IDENTITY),
                        m_B0(ViogtRank2Tensor2D::InitMethod::IDENTITY),
                        m_K(0.0),m_G(0.0),m_planeState(false){
    checkIfLargeStrain();
    initProperty(t_propPtr);
}
void NeoHookeanAbq2d::initProperty(nlohmann::json *t_propPtr){
    if(m_ifPropInit)return;
    if(t_propPtr->contains("K")&&t_propPtr->contains("G")){
        if(t_propPtr->at("K").is_number_float()&&t_propPtr->at("G").is_number_float()){
            m_K=t_propPtr->at("K"); m_G=t_propPtr->at("G");
        }
        else{
            MessagePrinter::printErrorTxt("properties K or G is not D float-point number");
            MessagePrinter::exitcfem();
        }
    }
    else if(t_propPtr->contains("E")&&t_propPtr->contains("nu")){
        if(t_propPtr->at("E").is_number_float()&&t_propPtr->at("nu").is_number_float()){
            double E=t_propPtr->at("E"), nu=t_propPtr->at("nu");
            m_K=ElasticConst::getKByE_Nu(E,nu);
            m_G=ElasticConst::getGByE_Nu(E,nu);
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
    m_ifPropInit=true;
}
void NeoHookeanAbq2d::updateMaterialBydudx(void *t_incStrainPtr,bool *t_converged){
    Rank2Tensor2d * incStrainPtr=(Rank2Tensor2d *)t_incStrainPtr;
    // update B
    m_B.setFromRank2Tensor2D(*incStrainPtr*m_B0*incStrainPtr->transpose());
    m_F=(*incStrainPtr)*m_F0;
    m_J=m_F.det();
    double m_J_pow=pow(m_J,-2.0/3.0);
    ViogtRank2Tensor2D Biso=m_B*m_J_pow;
    ViogtRank2Tensor2D Biso_dev=Biso-TensorConst2D::I*((Biso.trace()+m_J_pow)/3.0);
    m_S=(Biso_dev*m_G+TensorConst2D::I*(m_K*m_J*(m_J-1.0)))/m_J;
    m_T33=m_G*(m_J_pow-(Biso.trace()+m_J_pow)/3.0)+m_K*m_J*(m_J-1.0);
    *t_converged=true;
}
void NeoHookeanAbq2d::updateConvergence(){
    m_F0=m_F;
    m_B0=m_B;
}
void NeoHookeanAbq2d::getTangentModulus(void *t_incStrainPtr,void *t_D){
    int a=*(int *)t_incStrainPtr;
    int b=*(int *)t_D;
    if(a||b){}
    return;
}
void NeoHookeanAbq2d::getSpatialTangentModulus(void *t_incStrainPtr,MatrixXd *t_a){
    Rank2Tensor2d * incStrainPtr=(Rank2Tensor2d *)t_incStrainPtr;
    // update B
    m_B.setFromRank2Tensor2D(*incStrainPtr*m_B0*incStrainPtr->transpose());
    m_F=(*incStrainPtr)*m_F0;
    m_J=m_F.det();
    double m_J_pow=pow(m_J,-2.0/3.0);
    ViogtRank2Tensor2D Biso=m_B*m_J_pow;
    double Tr_Biso=Biso.trace()+m_J_pow;
    ViogtRank2Tensor2D Biso_dev=Biso-TensorConst2D::I*(Tr_Biso/3.0);
    m_S=(Biso_dev*m_G+TensorConst2D::I*(m_K*m_J*(m_J-1.0)))/m_J;
    m_T33=m_G*(m_J_pow-(Biso.trace()+m_J_pow)/3.0)+m_K*m_J*(m_J-1.0);
    double S33=m_T33/m_J;
    double p = (m_S(0,0)+m_S(1,1)+S33)/3.0;     // p=Sii/3
    ViogtRank2Tensor2D S_dev=m_S-TensorConst2D::I*p;
    const int n=4;
    MatrixXd item1(n,n,0.0), item2(n,n,0.0), item3(n,n,0.0), item4(n,n,0.0), item5(n,n,0.0);
    item1 = (TensorConst2D::IIDev*(2.0*m_G*Tr_Biso/(3.0*m_J))).toFullMatrix();
    item2 = (TensorConst2D::IISym*(-2.0*p)).toFullMatrix();
    item3 = (TensorConst2D::I.ijkl(S_dev)+S_dev.ijkl(TensorConst2D::I))*(-2.0/3.0);
    item4 = (TensorConst2D::IXI*(m_K*(2.0*m_J-1.0))).toFullMatrix();
    item5 = TensorConst2D::I.ikjl(m_S);
    *t_a=item1+item2+item3+item4+item5;
}
void NeoHookeanAbq2d::getMatVariable(ElementVariableType elmtVarType,void *elmtVarPtr){
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
    case ElementVariableType::LOGSTRAIN:{
        ViogtRank2Tensor2D strain=m_B.isotropicFunc(&log)*0.5;
        *(ViogtRank2Tensor2D *)elmtVarPtr=strain;
        break;
    }

    case ElementVariableType::PRESSURE:{
        double S33=m_T33*m_J;
        *(double *)elmtVarPtr=(m_S.trace()+S33)/(-3.0);
        break;
    }
    case ElementVariableType::VONMISES:{
        double S33=m_T33*m_J;
        double Sm=(m_S.trace()+S33)/3.0;
        *(double *)elmtVarPtr=sqrt(1.5*((m_S(0)-Sm)*(m_S(0)-Sm)+(m_S(1)-Sm)*(m_S(1)-Sm)+(S33-Sm)*(S33-Sm)+2*m_S(2)*m_S(2)));
        break;
    }
    default:
        MessagePrinter::printErrorTxt("Required material variale is not suuported by linearElastic material lib now");
        MessagePrinter::exitcfem();
        break;
    }    
}
void NeoHookeanAbq2d::getMatVariableArray(ElementVariableType elmtVarType,PetscScalar *elmtVarPtr){
    switch (elmtVarType)
    {
    case ElementVariableType::CAUCHYSTRESS:{
        double S33=m_T33*m_J;
        elmtVarPtr[0]=m_S(0);      elmtVarPtr[1]=m_S(1);    elmtVarPtr[2]=S33;
        elmtVarPtr[3]=m_S(2);      elmtVarPtr[4]=0.0;       elmtVarPtr[5]=0.0;
        break;
    }
    case ElementVariableType::JACOBIAN:
        *elmtVarPtr=m_J;
        break;
    case ElementVariableType::KIRCHOFFSTRESS:{
        double S33=m_T33*m_J;
        elmtVarPtr[0]=m_S(0)/m_J;      elmtVarPtr[1]=m_S(1)/m_J;    elmtVarPtr[2]=S33/m_J;
        elmtVarPtr[3]=m_S(2)/m_J;      elmtVarPtr[4]=0.0;           elmtVarPtr[5]=0.0;
        break;
    }
    case ElementVariableType::LOGSTRAIN:{
        ViogtRank2Tensor2D strain=m_B.isotropicFunc(&log)*0.5;
        elmtVarPtr[0]=strain(0);     elmtVarPtr[1]=strain(1);       elmtVarPtr[2]=0.0;
        elmtVarPtr[3]=strain(2)*2.0; elmtVarPtr[4]=0.0;             elmtVarPtr[5]=0.0;
        break;
    }
    case ElementVariableType::PRESSURE:{
        double S33=m_T33*m_J;
        *elmtVarPtr=(m_S.trace()+S33)/(-3.0);
        break;
    }
    case ElementVariableType::VONMISES:{
        double S33=m_T33*m_J;
        double Sm=(m_S.trace()+S33)/3.0;
        *elmtVarPtr=sqrt(1.5*((m_S(0)-Sm)*(m_S(0)-Sm)+(m_S(1)-Sm)*(m_S(1)-Sm)+(S33-Sm)*(S33-Sm)+2*m_S(2)*m_S(2)));
        break;
    }
    default:
        MessagePrinter::printErrorTxt("Required material variale is not suuported by linearElastic material lib now");
        MessagePrinter::exitcfem();
        break;
    }      
}
PetscErrorCode NeoHookeanAbq2d::checkIfLargeStrain(){
    if(!m_nLarge){
        MessagePrinter::printRankError("step->nLarge need to be set to true when use material Neo-Hookean (Abaqus version)");
        MessagePrinter::exitcfem();
    }
    return 0;
}