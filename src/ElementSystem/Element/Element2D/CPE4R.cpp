#include "ElementSystem/Element/Element2D/CPE4R.h"
#include "MathUtils/VectorXd.h"
const int CPE4R::m_dim=2;                     /**< element dimension*/
const int CPE4R::m_mDof_node=2;               /**< dof num per node*/
const int CPE4R::m_mNode=4;                   /**< a element's nodes number*/
const int CPE4R::m_mQPoint=1;
const int CPE4R::m_QPW=4.0;
const double CPE4R::m_HG_coeff=0.01;
ShpfunQuad4 CPE4R::m_shpfun=ShpfunQuad4(Vector2d(0.0,0.0)); /**< shape function relative computer*/
PetscErrorCode CPE4R::initElement(PetscInt t_elmt_rId, bool nLarge,MeshSystem *t_meshSysPtr, PetscScalar *elmtParamPtr){
    m_elmt_rId=t_elmt_rId;
    m_nLarge=nLarge;
    if(elmtParamPtr){}
    Vector2d elmt_coord0[m_mNode];
    t_meshSysPtr->getElmtNodeCoord(t_elmt_rId,0,elmt_coord0);
    Vector2d elmt_dNdr[m_mNode];
    m_shpfun.getDer2Nat(elmt_dNdr);
    Rank2Tensor2d dx0dr(Rank2Tensor2d::InitMethod::ZERO);
    for(int nodeI=0;nodeI<m_mNode;nodeI++){
        for(int m=0;m<m_mDof_node;m++){
            for(int n=0;n<m_mDof_node;n++){
                dx0dr(m,n)+=elmt_dNdr[nodeI](n)*elmt_coord0[nodeI](m);
            }
        }
    }    
    m_det_dx0dr=dx0dr.det();
    m_Q1[0]=0.; m_Q1[1]=0.; m_Q2[0]=0.; m_Q2[1]=0.;
    return 0;
}

PetscErrorCode CPE4R::getElmtInnerForce(void *t_elmtCoord2, void *t_elmtDofInc, VectorXd *t_elmtInnerForce, bool *t_converged){
    Vector2d *elmtCoord2=(Vector2d *)t_elmtCoord2;
    Vector2d *elmtDofInc=(Vector2d *)t_elmtDofInc;
    const int mBMatrix=3;   /** row number of B-matrix*/
    Vector2d elmtCoord1[m_mNode];
    Vector2d qPCoord;
    Vector2d dNdx[m_mNode];
    Vector2d dNdx2[m_mNode];    /**< derivate to last converged coord*/
    ViogtRank2Tensor2D stress;                                  /**< cauchy stress*/
    MatrixXd BMatrix(mBMatrix,m_mNode*m_mDof_node);             /**< the discrete symmetric gradient operation, B-matrix*/
    for (int i=0;i<m_mNode;i++){
        elmtCoord1[i]=elmtCoord2[i]+elmtDofInc[i];
    }
    this->m_shpfun.setRefCoords(elmtCoord1);
    this->m_shpfun.getDer2Ref(dNdx);
    if(!m_ifHGUpdated) updateHourglass(elmtCoord2,dNdx);
    this->getBMatrix(dNdx,m_mNode,m_mDof_node,&BMatrix);
    /**
     * update material
    */
    if(!m_nLarge){  // for small strain
        Rank2Tensor2d duIncdx;
        for(int nodeI=0;nodeI<m_mNode;nodeI++){
            for(int m=0;m<m_mDof_node;m++){
                for(int n=0;n<m_mDof_node;n++){
                    duIncdx(m,n)+=dNdx[nodeI](n)*elmtDofInc[nodeI](m);
                }
            }
        }
        m_matPtr->updateMaterialBydudx(&duIncdx,t_converged);
    }
    else{   // for large strain
        Rank2Tensor2d FInc(Rank2Tensor2d::InitMethod::ZERO);
        this->m_shpfun.setRefCoords(elmtCoord2);
        this->m_shpfun.getDer2Ref(dNdx2);
        for(int nodeI=0;nodeI<m_mNode;nodeI++){
            for(int m=0;m<m_mDof_node;m++){
                for(int n=0;n<m_mDof_node;n++){
                    FInc(m,n)+=dNdx2[nodeI](n)*elmtCoord1[nodeI](m);
                }
            }
        }
        m_matPtr->updateMaterialBydudx(&FInc,t_converged); 
    }
    // require cauchy stress
    m_matPtr->getMatVariable(ElementVariableType::CAUCHYSTRESS,&stress);
    if(!t_converged) return 1;
    // evaluate elemental volume
    double J=0;
    m_matPtr->getMatVariable(ElementVariableType::JACOBIAN,&J);
    double volume=m_QPW*J*m_det_dx0dr;
    *t_elmtInnerForce=BMatrix.transpose()*(stress*volume);
    // update current hourglass general force of hourglass
    for(int di=0;di<m_mDof_node;++di){
        m_Q1[di]=m_Q2[di];
        for(int nJ=0;nJ<m_mNode;++nJ){
            m_Q1[di]+=0.5*m_HG_coeff*m_Kmax2*m_gamma2[nJ]*elmtDofInc[nJ](di);
        }
    }
    // hourglass control force*/
    VectorXd f_HG(m_mNode*m_mDof_node,0.0);
    for(int nI=0;nI<m_mNode;++nI){
        for(int di=0;di<m_mDof_node;++di){
            int dof=nI*m_mDof_node+di;
            f_HG(dof)=0.5*m_gamma2[nI]*m_Q2[di];
            for(int nJ=0;nJ<m_mNode;++nJ){
                f_HG(dof)+=0.25*m_HG_coeff*m_Kmax2*m_gamma2[nI]*m_gamma2[nJ]*elmtDofInc[nJ](di);
                (*t_elmtInnerForce)(dof)+=f_HG(dof);
            }
        }
    }
    // for dubug
    // VectorXd f_HG_(m_mNode*m_mDof_node,0.0);
    // for(int nI=0;nI<m_mNode;++nI){
    //     for(int di=0;di<m_mDof_node;++di){
    //         int dof1=nI*m_mDof_node+di;
    //         f_HG_(dof1)=0.5*m_gamma2[nI]*m_Q2[di];
    //         for(int nJ=0;nJ<m_mNode;++nJ){
    //             for(int dj=0;dj<m_mDof_node;++dj){
    //                 int dof2=nJ*m_mDof_node+dj;
    //                 f_HG_(dof1)+=m_KHG2(dof1,dof2)*elmtDofInc[nJ](dj);
    //             }
    //         }
    //     }
    // }
    return 0;
}

PetscErrorCode CPE4R::getElmtStfMatrix(void *t_elmtCoord2, void *t_elmtDofInc, MatrixXd *t_stfMatrix){
    Vector2d *elmtCoord2=(Vector2d *)t_elmtCoord2;
    Vector2d *elmtDofInc=(Vector2d *)t_elmtDofInc;
    Vector2d elmtCoord1[m_mNode];
    /** cal current elmt coords*/
    for (int i=0;i<m_mNode;i++){
        elmtCoord1[i]=elmtCoord2[i]+elmtDofInc[i];
    }
    /** cal dNi/dx*/
    m_shpfun.setRefCoords(elmtCoord1);
    Vector2d dNdx[4];
    m_shpfun.getDer2Ref(dNdx);
    if(!m_ifHGUpdated) updateHourglass(elmtCoord2,dNdx);
    // evaluate elemental volume
    double J=0;
    m_matPtr->getMatVariable(ElementVariableType::JACOBIAN,&J);
    double volume=m_QPW*J*m_det_dx0dr;
    if(!m_nLarge){  // for small strain
        const int mBMatrix=3;   /** B-matrix's row num*/
        MatrixXd B(mBMatrix,m_mNode*m_mDof_node,0.0);
        getBMatrix(dNdx,m_mNode,m_mDof_node,&B);
        ViogtRank4Tensor2D D(ViogtRank4Tensor2D::InitMethod::ZERO);
        /** cal duInc/dx*/
        Rank2Tensor2d duIncdx;
        for(int nodeI=0;nodeI<m_mNode;nodeI++){
            for(int m=0;m<m_mDof_node;m++){
                for(int n=0;n<m_mDof_node;n++){
                    duIncdx(m,n)+=dNdx[nodeI](n)*elmtDofInc[nodeI](m);
                }
            }
        }
        m_matPtr->getTangentModulus(&duIncdx,&D);
        *t_stfMatrix=B.transpose()*D*B*volume+m_KHG2;
    }
    else{
        const int mGMatrix=4;   /** G-Matrix's row count*/
        MatrixXd G(mGMatrix,m_mNode*m_mDof_node,0.0);
        getGMatrix(dNdx,m_mNode,m_mDof_node,&G);
        MatrixXd a (mGMatrix,mGMatrix,0.0);
        /** cal Finc*/
        Vector2d dNdX2[m_mNode];
        Rank2Tensor2d Finc;
        m_shpfun.setRefCoords(elmtCoord2);
        m_shpfun.getDer2Ref(dNdX2);
        for(int nodeI=0;nodeI<m_mNode;nodeI++){
            for(int m=0;m<m_mDof_node;m++){
                for(int n=0;n<m_mDof_node;n++){
                    Finc(m,n)+=dNdX2[nodeI](n)*elmtCoord1[nodeI](m);
                }
            }
        }
        /** cal a*/
        m_matPtr->getSpatialTangentModulus(&Finc,&a);
        *t_stfMatrix=G.transpose()*a*G*volume+m_KHG2;
    }
    return 0;
}
PetscErrorCode CPE4R::getElmtWeightedVolumeInt(PetscScalar **t_valQPPtr,PetscScalar **t_valNodePtr, int t_mCpnt){
    double J=0;
    m_matPtr->getMatVariable(ElementVariableType::JACOBIAN,&J);
    double volume=m_QPW*J*m_det_dx0dr;
    double NPtr[m_mNode];     
    this->m_shpfun.getShpfunVal(NPtr);
    for(int nodeI=0;nodeI<m_mNode;++nodeI){// loop over every node in a elmt
        for(int cpntI=0;cpntI<t_mCpnt;++cpntI){// loop over every component
            t_valNodePtr[nodeI][cpntI]=t_valQPPtr[0][cpntI]*NPtr[nodeI]*volume;
        }
    }
    return 0;
}
void CPE4R::updateConvergence(){
    m_matPtr->updateConvergence();
    m_ifHGUpdated=false;
}
void CPE4R::getElmtVariableArray(ElementVariableType elmtVarType,PetscScalar **elmtVarPtr){
    m_matPtr->getMatVariableArray(elmtVarType,*elmtVarPtr);
}
double CPE4R::getKMax(double lame, double G, Vector2d t_dNdx2[4], double V2){
    if(lame){}
    double Kmax=0;
    for(int nI=0;nI<m_mNode;++nI){
        for(int di=0;di<m_mDof_node;++di){
            Kmax+=t_dNdx2[nI](di)*t_dNdx2[nI](di);
        }
    }
    // Kmax*=(lame+2*G)*V2;
    Kmax*=4*G*V2;
    return Kmax;
}
void CPE4R::updateHourglass(Vector2d *t_elmtCoord2, Vector2d *t_dNdx2){
    m_shpfun.getHGShpVec(t_dNdx2,t_elmtCoord2,m_gamma2);
    double lame=m_matPtr->getLame();
    double G=m_matPtr->getG();
    double J=0;
    m_matPtr->getMatVariable(ElementVariableType::JACOBIAN,&J);
    double volume=m_QPW*J*m_det_dx0dr;
    m_Kmax2=getKMax(lame,G,t_dNdx2,volume);
    for(int di=0;di<m_mDof_node;++di){
        m_Q2[di]=m_Q1[di];
    }
    for(int nI=0;nI<m_mNode;++nI){
        for(int nJ=0;nJ<m_mNode;++nJ){
            for(int ni=0;ni<m_mDof_node;++ni){
                int rowI=nI*m_mDof_node+ni;
                int colJ=nJ*m_mDof_node+ni;
                m_KHG2(rowI,colJ)=0.25*m_HG_coeff*m_Kmax2*m_gamma2[nI]*m_gamma2[nJ];
            }
        }
    }
    m_ifHGUpdated=true;
}