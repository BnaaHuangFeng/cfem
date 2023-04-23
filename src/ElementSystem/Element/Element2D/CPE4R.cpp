#include "ElementSystem/Element/Element2D/CPE4R.h"
#include "MathUtils/VectorXd.h"
const int CPE4R::m_dim=2;                     /**< element dimension*/
const int CPE4R::m_mDof_node=2;               /**< dof num per node*/
const int CPE4R::m_mNode=4;                   /**< a element's nodes number*/
const int CPE4R::m_mQPoint=1;
const int CPE4R::m_QPW=4.0;
ShpfunQuad4 CPE4R::m_shpfun=ShpfunQuad4(Vector2d(0.0,0.0)); /**< shape function relative computer*/

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
    this->getBMatrix(dNdx,m_mNode,m_mDof_node,&BMatrix);
    /**
     * update material
    */
    if(!m_nLarge){  // for small strain
        Rank2Tensor2d duIncdx(0.0);
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
        Rank2Tensor2d FInc=TensorConst2D::I;
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
    m_matPtr->getElementVariable(ElementVariableType::CAUCHYSTRESS,&stress);
    if(!t_converged) return 1;
    // evaluate elemental volume
    double J=0;
    m_matPtr->getElementVariable(ElementVariableType::JACOBIAN,&J);
    double volume=m_QPW*J;
    *t_elmtInnerForce=BMatrix.transpose()*(stress*volume);
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

    if(!m_nLarge){  // for small strain
        const int mBMatrix=3;   /** B-matrix's row num*/
        MatrixXd B(mBMatrix,m_mNode*m_mDof_node,0.0);
        getBMatrix(dNdx,m_mNode,m_mDof_node,&B);
        ViogtRank4Tensor2D D(0.0);
        /** cal duInc/dx*/
        Rank2Tensor2d duIncdx(0.0);
        for(int nodeI=0;nodeI<m_mNode;nodeI++){
            for(int m=0;m<m_mDof_node;m++){
                for(int n=0;n<m_mDof_node;n++){
                    duIncdx(m,n)+=dNdx[nodeI](n)*elmtDofInc[nodeI](m);
                }
            }
        }
        m_matPtr->getTangentModulus(&duIncdx,&D);
        *t_stfMatrix=B.transpose()*D*B;
    }
    else{
        const int mGMatrix=4;   /** G-Matrix's row count*/
        MatrixXd G(mGMatrix,m_mNode*m_mDof_node,0.0);
        getGMatrix(dNdx,m_mNode,m_mDof_node,&G);
        MatrixXd a (mGMatrix,mGMatrix,0.0);
        /** cal Finc*/
        Vector2d dNdX2[m_mNode];
        Rank2Tensor2d Finc(0.0);
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
        *t_stfMatrix=G.transpose()*a*G;
    }
    return 0;
}