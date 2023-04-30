#include "ElementSystem/ElementSystem.h"
#include "ElementSystem/Element/Element2D/ElementPack2d.h"
#include "MaterialSystem/MaterialPack2d.h"
#include "MathUtils/VectorXd.h"
ElementSystem::ElementSystem():
    m_timerPtr(nullptr),m_ifElmtDesRead(false),m_ifMatDesRead(false),
    m_ifSetMeshSysPtr(nullptr),m_ifAssignElmtType(nullptr),m_ifAssignMatype(nullptr){
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);    
}
ElementSystem::ElementSystem(Timer* timerPtr,ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr):
    m_ifSetMeshSysPtr(nullptr),m_ifAssignElmtType(nullptr),m_ifAssignMatype(nullptr){
    m_timerPtr=timerPtr;
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum); 
    /** int elemt description  **/
    /****************************/
    readElmtDes(elmtDesPtr);
    /** int material description  **/
    /*******************************/
    readMatDes(matDesPtr);
}
ElementSystem::~ElementSystem(){
    /** delete every elmt item **/
    /****************************/
    const int mElmts=m_elmtPtrs.size();
    for(int eI=0;eI<mElmts;eI++){//loop over elmt ptrs
        if(m_elmtPtrs[eI]) delete m_elmtPtrs[eI];
    }
}
PetscErrorCode ElementSystem::init(ElementDescription *elmtDesPtr,MaterialDescription *matDesPtr,MeshSystem *meshPtr){
    /** int elemt description  **/
    /****************************/
    readElmtDes(elmtDesPtr);
    /** int material description  **/
    /*******************************/
    readMatDes(matDesPtr);  
    init(meshPtr);
}
PetscErrorCode ElementSystem::init(MeshSystem *meshPtr){
    setMeshSysPtr(meshPtr);
    if(!m_ifElmtDesRead){
        MessagePrinter::printErrorTxt("need to read element description before init element system!");
        MessagePrinter::exitcfem();
    }
    else if(!m_ifMatDesRead){
        MessagePrinter::printErrorTxt("need to read material description before init element system!");
        MessagePrinter::exitcfem();
    }
    else{
        m_elmtPtrs.resize(meshPtr->m_mElmts_p,nullptr);
        assignElmtType();
        assignMatType();
    }
    return 0;
}
bool ElementSystem::checkElmtsAssigment(){
    PetscInt mElmt_p=m_elmtPtrs.size(); /**< elmt num in this rank*/
    for(PetscInt i=0;i<mElmt_p;++i){
        if(!m_elmtPtrs[i]){
            MessagePrinter::printErrorTxt("a elemnt (global id = "+
                                to_string(m_meshSysPtr->elmtRId2GId(i))+
                                ") has lost the element type assigment.");
            MessagePrinter::exitcfem();
        }
        if(!(m_elmtPtrs[i]->m_matPtr)){
            MessagePrinter::printErrorTxt("a elemnt (global id = "+
                                to_string(m_meshSysPtr->elmtRId2GId(i))+
                                ") has lost the material assigment.");
            MessagePrinter::exitcfem();            
        }
    }
    return true;
}
PetscErrorCode ElementSystem::assignElmtType(){
    const int mElmtType=m_elmtTypeNames.size();
    for(int elmtTypeI=0;elmtTypeI<mElmtType;++elmtTypeI){// loop over every elmt type
        vector<PetscInt> &elmtSet=m_meshSysPtr->m_setManager.getSet(
            m_elmtAssignSetNames[elmtTypeI],SetType::ELEMENT);
        PetscInt mElmtInSet=elmtSet.size(); /**< num of elmts in this set*/
        for(PetscInt i=0;i<mElmtInSet;i++){// loop over every elmt of this elmt type
            switch (m_elmtTypes[elmtTypeI])
            {
            case ElementType::CPE4R:
                m_elmtPtrs[elmtSet[i]]=new CPE4R(m_nLarge);
                break;
            default:
                MessagePrinter::printErrorTxt("can not create a element of unsupported element type.");
                MessagePrinter::exitcfem();
                break;
            }
            m_elmtPtrs[elmtSet[i]]->initElement(elmtSet[i],m_nLarge,nullptr);
        }
    }
    m_ifAssignElmtType=true;
    return 0;
}
PetscErrorCode ElementSystem::assignMatType(){
    const int mMatType=m_matTypeNames.size();
    for(int matTypeId=0;matTypeId<mMatType;++matTypeId){// loop over every material type
        vector<PetscInt> &elmtSet=m_meshSysPtr->m_setManager.getSet(
            m_materialAssignSetNames[matTypeId],SetType::ELEMENT);
        PetscInt mElmtInSet=elmtSet.size(); /**< num of elmts in this materia; type*/
        for(PetscInt i=0;i<mElmtInSet;i++){
            switch(m_matTypes[matTypeId]){
                case MaterialType::LINEARELASTIC:
                    m_elmtPtrs[elmtSet[i]]->m_matPtr=new LinearElasticMat2D(m_nLarge);
            }
            m_elmtPtrs[elmtSet[i]]->m_matPtr->initProperty(&(m_properties[matTypeId]));
        }
    }
    m_ifAssignMatype=true;
    return 0;
}
PetscErrorCode ElementSystem::checkInit(){
    if(!m_ifElmtDesRead){
        MessagePrinter::printErrorTxt("ElemntSystem: it has not read element description.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifMatDesRead){
        MessagePrinter::printErrorTxt("ElemntSystem: it has not read material description.");
        MessagePrinter::exitcfem();
    }
    if(!m_ifSetMeshSysPtr){
        MessagePrinter::printErrorTxt("ElemntSystem: its relied mesh system has not yet been set.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifAssignElmtType){
        MessagePrinter::printErrorTxt("ElemntSystem: it has not assign element type.");
        MessagePrinter::exitcfem();          
    }
    if(!m_ifAssignMatype){
        MessagePrinter::printErrorTxt("ElemntSystem: it has not assign material type.");
        MessagePrinter::exitcfem();          
    }
    checkElmtsAssigment();
    return 0;    
}
PetscErrorCode ElementSystem::assembleAMatrix(Vec *t_uInc1Ptr, Mat *t_AMatrixPtr){
    const int MNodeElmt2d=9, MNodeElmt3d=27;
    Vector2d coord2Ptr2d[MNodeElmt2d], uIncPtr2d[MNodeElmt2d];
    Vector3d coord2Ptr3d[MNodeElmt3d], uIncPtr3d[MNodeElmt3d];
    /**open access to node variable Vec***********************************************/
    /*********************************************************************************/
    m_meshSysPtr->openNodeVariableVec(NodeVariableType::COORD,&(m_meshSysPtr->m_nodes_coord2),2,VecAccessMode::READ);
    m_meshSysPtr->openNodeVariableVec(NodeVariableType::UINC,t_uInc1Ptr,1,VecAccessMode::READ);
    for(PetscInt eI=0;eI<m_meshSysPtr->m_mElmts_p;eI++){// loop over every element in this rank
        element *elmtPtr=m_elmtPtrs[eI];
        int mDofInElmt=elmtPtr->getDofNum();
        int dim=elmtPtr->getDim();
        int mNode=elmtPtr->getNodeNum();
        int mDofPerNode=elmtPtr->getDofPerNode();
        Vector *coord2Ptr, *uIncPtr;
        if(dim==2){
            coord2Ptr=coord2Ptr2d;
            uIncPtr=uIncPtr2d;
        }
        else if(dim==3){
            MessagePrinter::printErrorTxt("dim = 3 is not supported now");
            MessagePrinter::exitcfem();
        }
        else{
            MessagePrinter::printErrorTxt("dim = "+to_string(dim)+" is not supported now");
            MessagePrinter::exitcfem();
        }
        m_meshSysPtr->getElmtNodeCoord(elmtPtr->m_elmt_rId,2,coord2Ptr);
        m_meshSysPtr->getElmtNodeUInc(elmtPtr->m_elmt_rId,1,uIncPtr);
        MatrixXd AMatrixElmt=MatrixXd(mDofInElmt,mDofInElmt,0.0);   /**< elmt's jacobian matrix*/
        bool ifMatUpdateConvergerd=false;
        elmtPtr->getElmtStfMatrix(coord2Ptr,uIncPtr,&AMatrixElmt);
        m_meshSysPtr->addElmtAMatrix(elmtPtr->m_elmt_rId,&AMatrixElmt,t_AMatrixPtr);
    }
    /** assemble or restore global Mat**************************************/
    /***********************************************************************/
    PetscCall(MatAssemblyBegin(*t_AMatrixPtr,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*t_AMatrixPtr,MAT_FINAL_ASSEMBLY));
    m_meshSysPtr->closeNodeVariableVec(NodeVariableType::COORD,&(m_meshSysPtr->m_nodes_coord2),2,VecAccessMode::READ);
    m_meshSysPtr->closeNodeVariableVec(NodeVariableType::UINC,t_uInc1Ptr,1,VecAccessMode::READ);
}
PetscErrorCode ElementSystem::assemblRVec(Vec *t_uInc1Ptr, Vec *t_RVecPtr){
    const int MNodeElmt2d=9, MNodeElmt3d=27;
    Vector2d coord2Ptr2d[MNodeElmt2d], uIncPtr2d[MNodeElmt2d], fIVector2d[MNodeElmt2d];
    Vector3d coord2Ptr3d[MNodeElmt3d], uIncPtr3d[MNodeElmt3d], fIVector3d[MNodeElmt3d];
    /**open access to node variable Vec***********************************************/
    /*********************************************************************************/
    m_meshSysPtr->openNodeVariableVec(NodeVariableType::COORD,&(m_meshSysPtr->m_nodes_coord2),2,VecAccessMode::READ);
    m_meshSysPtr->openNodeVariableVec(NodeVariableType::UINC,t_uInc1Ptr,1,VecAccessMode::READ);
    m_meshSysPtr->openNodeVariableVec(NodeVariableType::RESIDUAL,t_RVecPtr,1,VecAccessMode::WRITE);
    for(PetscInt eI=0;eI<m_meshSysPtr->m_mElmts_p;eI++){// loop over every element in this rank
        element *elmtPtr=m_elmtPtrs[eI];
        int mDofInElmt=elmtPtr->getDofNum();
        int dim=elmtPtr->getDim();
        int mNode=elmtPtr->getNodeNum();
        int mDofPerNode=elmtPtr->getDofPerNode();
        Vector *coord2Ptr, *uIncPtr, *fIVector;
        if(dim==2){
            coord2Ptr=coord2Ptr2d;
            uIncPtr=uIncPtr2d;
            fIVector=fIVector2d;
        }
        else if(dim==3){
            MessagePrinter::printErrorTxt("dim = 3 is not supported now");
            MessagePrinter::exitcfem();
        }
        else{
            MessagePrinter::printErrorTxt("dim = "+to_string(dim)+" is not supported now");
            MessagePrinter::exitcfem();
        }
        m_meshSysPtr->getElmtNodeCoord(elmtPtr->m_elmt_rId,2,coord2Ptr);
        m_meshSysPtr->getElmtNodeUInc(elmtPtr->m_elmt_rId,1,uIncPtr);
        VectorXd fI=VectorXd(mDofInElmt,0.0);
        bool ifMatUpdateConvergerd=false;
        elmtPtr->getElmtInnerForce(coord2Ptr,uIncPtr,&fI,&ifMatUpdateConvergerd);
        if(!ifMatUpdateConvergerd) return 7890; // material updation failed
        for(int nodeI=0;nodeI<mNode;nodeI++){
            for(int dofI=0;dofI<mDofPerNode;dofI++){
                fIVector[nodeI](dofI)=fI(nodeI*mDofPerNode+dofI);
            }
        }
        m_meshSysPtr->addElmtResidual(elmtPtr->m_elmt_rId,fIVector,t_RVecPtr);
    }
    /** assemble or restore global Vec**************************************/
    /***********************************************************************/
    m_meshSysPtr->closeNodeVariableVec(NodeVariableType::COORD,&(m_meshSysPtr->m_nodes_coord2),2,VecAccessMode::READ);
    m_meshSysPtr->closeNodeVariableVec(NodeVariableType::UINC,t_uInc1Ptr,1,VecAccessMode::READ);
    m_meshSysPtr->closeNodeVariableVec(NodeVariableType::RESIDUAL,t_RVecPtr,1,VecAccessMode::WRITE);
}
void ElementSystem::readElmtDes(ElementDescription *elmtDesPtr){
    if(m_ifElmtDesRead)return;
    m_elmtTypeNames=elmtDesPtr->s_names;
    m_elmtTypes=elmtDesPtr->s_elmtTypes;
    m_elmtAssignSetNames=elmtDesPtr->s_setNames;
    m_nLarge=elmtDesPtr->s_nLarge;
    m_ifElmtDesRead=true;
}
void ElementSystem::readMatDes(MaterialDescription *matDesPtr){
    if(m_ifMatDesRead)return;
    m_matTypeNames=matDesPtr->s_names;
    m_matTypes=matDesPtr->s_matType;
    m_properties=matDesPtr->s_properties;
    m_materialAssignSetNames=matDesPtr->s_setName;    
    m_ifMatDesRead=true;
}