#include "PostProcessSystem/PostStructured2d.h"
#include "MaterialSystem/ElmtVarInfo.h"
#include "MeshSystem/StructuredMesh2D.h"
PostStructured2d::PostStructured2d(OutputDescription *t_outputDesPtr):PostProcessSystem(t_outputDesPtr){};

PostStructured2d::PostStructured2d(OutputDescription *t_outputDesPtr,MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr, LoadController *t_loadCtrlPtr):
    PostProcessSystem(t_outputDesPtr,t_meshSysPtr,t_elmtSysPtr,t_loadCtrlPtr){
}
PostStructured2d::~PostStructured2d(){
    DMDestroy(&m_dmScalar);
    DMDestroy(&m_dmRank2Tensor2d);
    DMDestroy(&m_dmRank2Tensor3d);
}
PetscErrorCode PostStructured2d::clear(){
    PetscCall(DMDestroy(&m_dmScalar));
    PetscCall(DMDestroy(&m_dmRank2Tensor2d));
    PetscCall(DMDestroy(&m_dmRank2Tensor3d));
    return 0;
}
PetscErrorCode PostStructured2d::init(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr){
    if(!m_ifMeshSysSet){
        m_meshSysPtr=t_meshSysPtr;
        m_ifMeshSysSet=true;
    }
    if(!m_ifElmtSysSet){
        m_elmtSysPtr=t_elmtSysPtr;
        m_ifElmtSysSet=true;
    }
    init();
    return 0;
}
PetscErrorCode PostStructured2d::init(){
    if(!m_ifMeshSysSet){
        MessagePrinter::printErrorTxt("PostStructured2d: need to set mesh system it rely before init PostStructured2d");
        MessagePrinter::exitcfem();
    }
    if(!m_ifElmtSysSet){
        MessagePrinter::printErrorTxt("PostStructured2d: need to set element system it rely before init PostStructured2d");
        MessagePrinter::exitcfem();
    }    
    initDm();
    m_ifDmInit=true;
    /** Create global vector*/
    PetscCall(DMCreateGlobalVector(m_dmScalar,&m_proj_weight));
    PetscCall(VecZeroEntries(m_proj_weight));
    PetscCall(DMCreateGlobalVector(m_meshSysPtr->m_dm,&m_his_node_vec));
    PetscCall(VecZeroEntries(m_his_node_vec));
    initBuffer();
    return 0;
}
PetscErrorCode PostStructured2d::initBuffer(){
    int mHisNodeVar=m_infoForHisOut.varInNodeVec.size();
    for(int varI=0;varI<mHisNodeVar;++varI){// loop over every historic nodal variavle
        PetscInt mCpnt=0;           /**variable component num*/
        PetscInt varNum=0;          /**variable num of a frame*/
        switch (m_infoForHisOut.varCpntIndInNodeVec[varI])
        {
        case -1:
            mCpnt=m_infoForHisOut.varCpntInNodeVec[varI];
            break;
        default:
            mCpnt=1;
            break;
        }
        switch (m_infoForHisOut.outputFormInNodeVec[varI])
        {
        case VarOutputForm::ALL:
            varNum=m_meshSysPtr->m_setManager.getSet(m_infoForHisOut.setNameInNodeVec[varI],SetType::NODE).size();
            break;
        case VarOutputForm::ANY:
        case VarOutputForm::SUM:
            varNum=1;
            break;
        default:
            break;
        }
        PetscScalar *buffer=new PetscScalar[m_hisBuffLen*varNum*mCpnt];
        PetscScalar *timeBuffer=new PetscScalar[m_hisBuffLen];
        m_infoForHisOut.bufferInNodeVec.push_back(buffer);
        m_infoForHisOut.timeBuffInNodeVec.push_back(timeBuffer);
        m_infoForHisOut.bufferFrameNumInNodeVec.push_back(0);
        m_infoForHisOut.dataNumPerFrameInNodeVec.push_back(varNum);
        m_infoForHisOut.mCpntPerDataInNodeVec.push_back(mCpnt);
    }
    int mHisElmtVar=m_infoForHisOut.varInElmtVec.size();
    for(int varI=0;varI<mHisElmtVar;++varI){// loop over every historic elemental variavle
        PetscInt mCpnt=0;           /**variable component num*/
        PetscInt varNum=0;          /**variable num of a frame*/
        switch (m_infoForHisOut.varCpntIndInElmtVec[varI])
        {
        case -1:
            mCpnt=m_infoForHisOut.varCpntInElmtVec[varI];
            break;
        default:
            mCpnt=1;
            break;
        }
        switch (m_infoForHisOut.outputFormInElmtVec[varI])
        {
        case VarOutputForm::ALL:
            varNum=m_meshSysPtr->m_setManager.getSet(m_infoForHisOut.setNameInElmtVec[varI],SetType::ELEMENT).size();
            break;
        case VarOutputForm::ANY:
        case VarOutputForm::SUM:
            varNum=1;
            break;
        default:
            break;
        }
        PetscScalar *buffer=new PetscScalar[m_hisBuffLen*varNum*mCpnt];
        PetscScalar *timeBuffer=new PetscScalar[m_hisBuffLen];
        m_infoForHisOut.bufferInElmtVec.push_back(buffer);
        m_infoForHisOut.timeBuffInElmtVec.push_back(timeBuffer);
        m_infoForHisOut.bufferFrameNumInElmtVec.push_back(0);
        m_infoForHisOut.dataNumPerFrameInElmtVec.push_back(varNum);
        m_infoForHisOut.mCpntPerDataInElmtVec.push_back(mCpnt);
    }
    return 0;
}
PetscErrorCode PostStructured2d::checkInit(){
    if(!m_ifMeshSysSet){
        MessagePrinter::printErrorTxt("PostStructured2d: need to set mesh system it rely before init PostStructured2d");
        MessagePrinter::exitcfem();
    }
    if(!m_ifElmtSysSet){
        MessagePrinter::printErrorTxt("PostStructured2d: need to set element system it rely before init PostStructured2d");
        MessagePrinter::exitcfem();
    }    
    if(!m_ifDmInit){
        MessagePrinter::printErrorTxt("PostStructured2d: DMs were not inited.");
        MessagePrinter::exitcfem();        
    }
    return 0;
}
PetscErrorCode PostStructured2d::output(int t_increI, PetscScalar t_t){
    outputFieldVariable(t_increI,t_t);
    outputHisVariable(t_increI,t_t);
    return 0;
}
PetscErrorCode PostStructured2d::projElmtVariable(ElementVariableType varType){
    if(m_ifProjVec){
        MessagePrinter::printErrorTxt("PostStructured2d: need to destroy Vec before a new projection.");
        MessagePrinter::exitcfem();             
    }
    auto it=ElmtVarInfo::elmtVarCpntNum.find(varType);
    int mCpnt=it->second;
    DM *dmPtr=nullptr;
    const int mNode=4;
    const int maxQPoint=4;  /**< max qpoint num*/
    PetscScalar **elmtVarPtr=new PetscScalar *[maxQPoint];
    PetscScalar **nodeVarPtr=new PetscScalar *[mNode];
    PetscScalar **elmtOnePtr=new PetscScalar *[maxQPoint];
    PetscScalar **nodeWeightPtr=new PetscScalar *[mNode];
    for(int nodeI=0;nodeI<mNode;++nodeI){
        nodeVarPtr[nodeI]=new PetscScalar[mCpnt];
        nodeWeightPtr[nodeI]=new PetscScalar[1];
        nodeWeightPtr[nodeI][0]=0.0;
    }
    for(int qpI=0;qpI<maxQPoint;++qpI){
        elmtVarPtr[qpI]=new PetscScalar[mCpnt];
        elmtOnePtr[qpI]=new PetscScalar[1];
        elmtOnePtr[qpI][0]=1.0;
    }
    getDmPtrByCpntNum(&dmPtr,mCpnt);
    PetscCall(DMCreateGlobalVector(*dmPtr,&m_proj_vec));
    m_ifProjVec=true;
    openNodeVariableVec(&m_proj_vec,&m_proj_val_local,&m_array_proj_val,mCpnt,VecAccessMode::WRITE);
    openNodeVariableVec(&m_proj_weight,&m_proj_weight_local,&m_array_proj_weight,1,VecAccessMode::WRITE);
    for(vector<element *>::iterator elmtIt=m_elmtSysPtr->m_elmtPtrs.begin();
    elmtIt!=m_elmtSysPtr->m_elmtPtrs.end();++elmtIt){// loop over every elmt in this rank
        PetscInt elmtRId=(*elmtIt)->m_elmt_rId;
        (*elmtIt)->getElmtWeightedVolumeInt(elmtOnePtr,nodeWeightPtr,1);
        addElmtVec(elmtRId,m_array_proj_weight,nodeWeightPtr,1);
    }
    closeNodeVariableVec(&m_proj_weight,&m_proj_weight_local,&m_array_proj_weight,1,VecAccessMode::WRITE);
    openNodeVariableVec(&m_proj_weight,&m_proj_weight_local,&m_array_proj_weight,1,VecAccessMode::READ);
    const int offset[mNode][2]={{0,0},{1,0},{1,1},{0,1}};
    for(vector<element *>::iterator elmtIt=m_elmtSysPtr->m_elmtPtrs.begin();
    elmtIt!=m_elmtSysPtr->m_elmtPtrs.end();++elmtIt){// loop over every elmt in this rank
        PetscInt elmtRId=(*elmtIt)->m_elmt_rId;
        (*elmtIt)->getElmtVariableArray(varType,elmtVarPtr);
        (*elmtIt)->getElmtWeightedVolumeInt(elmtVarPtr,nodeVarPtr,mCpnt);
        PetscInt elmtXI, elmtYI;
        getElmtDmdaIndByRId(elmtRId,&elmtXI,&elmtYI);
        for(int nodeI=0;nodeI<mNode;nodeI++){// loop over node in a elmt
            PetscInt nodeXI, nodeYI;
            nodeXI=elmtXI+offset[nodeI][0]; nodeYI=elmtYI+offset[nodeI][1];
            for(int cpntI=0;cpntI<mCpnt;cpntI++){// loop over component in a node
                nodeVarPtr[nodeI][cpntI]=nodeVarPtr[nodeI][cpntI]/m_array_proj_weight[nodeYI][nodeXI][0];
            }
        }
        addElmtVec(elmtRId,m_array_proj_val,nodeVarPtr,mCpnt);
    }
    closeNodeVariableVec(&m_proj_vec,&m_proj_val_local,&m_array_proj_val,mCpnt,VecAccessMode::WRITE);
    closeNodeVariableVec(&m_proj_weight,&m_proj_weight_local,&m_array_proj_weight,1,VecAccessMode::READ);
    return 0;
}
PetscErrorCode PostStructured2d::projVecClean(){
    if(m_ifProjVec){
        PetscCall(VecDestroy(&m_proj_vec));
        m_ifProjVec=false;
    }
    return 0;
}

PetscErrorCode PostStructured2d::genNodeVariable(NodeVariableType varType){
    switch (varType)
    {
    case NodeVariableType::U:
        openNodeVariableVec(&m_meshSysPtr->m_nodes_u2,&m_his_node_vec_local,&m_array_his_node,2,VecAccessMode::READ);
        break;
    case NodeVariableType::RF:
        m_elmtSysPtr->assemblRVec(&m_meshSysPtr->m_node_residual2,&m_his_node_vec);
        m_loadCtrlPtr->applyLoad(m_loadCtrlPtr->m_factor2,&m_his_node_vec);
        openNodeVariableVec(&m_his_node_vec,&m_his_node_vec_local,&m_array_his_node,2,VecAccessMode::READ);
        // for debug
        // PetscCall(VecView(m_his_node_vec,PETSC_VIEWER_STDOUT_WORLD));
        break;
    default:
        break;
    }
    return 0;
}
PetscErrorCode PostStructured2d::restoreNodeVariable(NodeVariableType varType){
    switch (varType)
    {
    case NodeVariableType::U:
        closeNodeVariableVec(&m_meshSysPtr->m_nodes_u2,&m_his_node_vec_local,&m_array_his_node,2,VecAccessMode::READ);
        break;
    case NodeVariableType::RF:
        closeNodeVariableVec(&m_his_node_vec,&m_his_node_vec_local,&m_array_his_node,2,VecAccessMode::READ);
    default:
        break;
    }
    m_array_his_node=nullptr;
    return 0;    
}
PetscErrorCode PostStructured2d::initDm(){
    PetscCall(DMDAGetLocalInfo(m_meshSysPtr->m_dm,&m_dmInfo));
    // get sequential Vec storing ym of every rank
    VecScatter scatter;
    Vec ymVecPar;
    Vec ymVecSeq;
    PetscScalar *ymArray=new PetscScalar[m_rankNum];
    PetscCall(VecCreate(PETSC_COMM_WORLD,&ymVecPar));
    PetscCall(VecSetType(ymVecPar,VECMPI));
    PetscCall(VecSetSizes(ymVecPar,1,m_rankNum)); 
    PetscCall(VecSetValue(ymVecPar,m_rank,m_dmInfo.ym,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(ymVecPar));
    PetscCall(VecAssemblyEnd(ymVecPar));
    PetscCall(VecScatterCreateToAll(ymVecPar,&scatter,&ymVecSeq));
    PetscCall(VecScatterBegin(scatter,ymVecPar,ymVecSeq,INSERT_VALUES,SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter,ymVecPar,ymVecSeq,INSERT_VALUES,SCATTER_FORWARD));
    PetscInt *rankIds=new PetscInt[m_rankNum];
    for(int i=0;i<m_rankNum;++i){
        rankIds[i]=i;
    }
    PetscCall(VecGetValues(ymVecSeq,m_rankNum,rankIds,ymArray));
    for(int i=0;i<m_rankNum;++i){
        rankIds[i]=(PetscInt)ymArray[i];
    }
    // crate m_dmScalar*************************************************************/
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
            DMDA_STENCIL_BOX,           /** stencil type. Use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.*/
            m_dmInfo.mx,m_dmInfo.my,    /** global dimension in each direction of the array*/
            1,m_rankNum,                /** corresponding number of processors in each dimension (or PETSC_DECIDE to have calculated)*/
            1,                          /** number of degrees of freedom per node*/
            1,                          /** stencil width*/
            &(m_dmInfo.mx),rankIds,     /** arrays containing the number of nodes in each cell along the x and y coordinates, or NULL.
                                         *  If non-null, these must be of length as m and n, and the corresponding m and n cannot be 
                                         *  PETSC_DECIDE. The sum of the lx[] entries must be M, and the sum of the ly[] entries must be N.*/
            &m_dmScalar                 /** DM pointer*/
            ));   
    PetscCall(DMSetFromOptions(m_dmScalar));
    PetscCall(DMSetUp(m_dmScalar));

    // crate m_dmRank2Tensor2d*************************************************************/
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
            DMDA_STENCIL_BOX,           /** stencil type. Use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.*/
            m_dmInfo.mx,m_dmInfo.my,    /** global dimension in each direction of the array*/
            1,m_rankNum,                /** corresponding number of processors in each dimension (or PETSC_DECIDE to have calculated)*/
            3,                          /** number of degrees of freedom per node*/
            1,                          /** stencil width*/
            &(m_dmInfo.mx),rankIds,     /** arrays containing the number of nodes in each cell along the x and y coordinates, or NULL.
                                         *  If non-null, these must be of length as m and n, and the corresponding m and n cannot be 
                                         *  PETSC_DECIDE. The sum of the lx[] entries must be M, and the sum of the ly[] entries must be N.*/
            &m_dmRank2Tensor2d          /** DM pointer*/
            ));   
    PetscCall(DMSetFromOptions(m_dmRank2Tensor2d));
    PetscCall(DMSetUp(m_dmRank2Tensor2d));

    // crate DMDA*************************************************************/
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
            DMDA_STENCIL_BOX,           /** stencil type. Use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.*/
            m_dmInfo.mx,m_dmInfo.my,    /** global dimension in each direction of the array*/
            1,m_rankNum,                /** corresponding number of processors in each dimension (or PETSC_DECIDE to have calculated)*/
            6,                          /** number of degrees of freedom per node*/
            1,                          /** stencil width*/
            &(m_dmInfo.mx),rankIds,     /** arrays containing the number of nodes in each cell along the x and y coordinates, or NULL.
                                         *  If non-null, these must be of length as m and n, and the corresponding m and n cannot be 
                                         *  PETSC_DECIDE. The sum of the lx[] entries must be M, and the sum of the ly[] entries must be N.*/
            &m_dmRank2Tensor3d          /** DM pointer*/
            ));   
    PetscCall(DMSetFromOptions(m_dmRank2Tensor3d));
    PetscCall(DMSetUp(m_dmRank2Tensor3d));
    delete []ymArray;
    delete []rankIds;
    PetscCall(VecDestroy(&ymVecPar));
    PetscCall(VecDestroy(&ymVecSeq));
    PetscCall(VecScatterDestroy(&scatter));
    return 0;
}
PetscErrorCode PostStructured2d::openNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrPtr, PetscInt mCpnt, VecAccessMode mode){
    DM *dmPtr;
    getDmPtrByCpntNum(&dmPtr,mCpnt);
    if(*arrayPtrPtr){
        MessagePrinter::printErrorTxt("PostStructured2d: arrayPtr need to pointer to null first");
        MessagePrinter::exitcfem();          
    }
    else{
        PetscCall(DMGetLocalVector(*dmPtr,localVecPtr));
        if(mode==VecAccessMode::READ){
            PetscCall(DMDAVecGetArrayDOFRead(*dmPtr,*localVecPtr,arrayPtrPtr));
            PetscCall(DMGlobalToLocal(*dmPtr,*globalVecPtr,INSERT_VALUES,*localVecPtr));
        }
        else if(mode==VecAccessMode::WRITE){
            PetscCall(DMDAVecGetArrayDOFWrite(*dmPtr,*localVecPtr,arrayPtrPtr));
            PetscCall(VecZeroEntries(*localVecPtr));  
        }        
    }
    return 0;
}
PetscErrorCode PostStructured2d::closeNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrPtr, PetscInt mCpnt, VecAccessMode mode){
    DM *dmPtr;
    getDmPtrByCpntNum(&dmPtr,mCpnt);
    if(*arrayPtrPtr){
        if(mode==VecAccessMode::READ){
            PetscCall(DMDAVecRestoreArrayDOFRead(*dmPtr,*localVecPtr,arrayPtrPtr));
            PetscCall(DMRestoreLocalVector(*dmPtr,localVecPtr));
        }        
        else if(mode==VecAccessMode::WRITE){
            PetscCall(DMDAVecRestoreArrayDOFWrite(*dmPtr,*localVecPtr,arrayPtrPtr));
            PetscCall(VecZeroEntries(*globalVecPtr));  
            PetscCall(DMLocalToGlobal(*dmPtr,*localVecPtr,ADD_VALUES,*globalVecPtr));
            PetscCall(DMRestoreLocalVector(*dmPtr,localVecPtr));            
        }
    }
    else{
        MessagePrinter::printErrorTxt("PostStructured2d: array is point to null, can not be restored");
        MessagePrinter::exitcfem();
    }
    *arrayPtrPtr=nullptr;
    return 0;
}



PetscErrorCode PostStructured2d::getDmPtrByCpntNum(DM **dmPtrAdr,int mCpnt){
    switch (mCpnt)
    {
    case 1:
        *dmPtrAdr=&m_dmScalar;
        break;
    case 4:
        *dmPtrAdr=&m_dmRank2Tensor2d;
        break;
    case 2:
        *dmPtrAdr=&m_meshSysPtr->m_dm;
        break;
    case 6:
        *dmPtrAdr=&m_dmRank2Tensor3d;
        break;
    default:
        MessagePrinter::printErrorTxt("PostStructured2d: unsupported projected vector components number.");
        MessagePrinter::exitcfem();      
        break;
    }
    return 0;
}
void PostStructured2d::getElmtDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr){
    *yIPtr=rId/(m_dmInfo.mx-1)+m_dmInfo.ys;
    *xIPtr=rId%(m_dmInfo.mx-1)+m_dmInfo.xs;    
}
void PostStructured2d::getNodeDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr){
    *yIPtr=rId/m_dmInfo.mx+m_dmInfo.ys;
    *xIPtr=rId%m_dmInfo.mx+m_dmInfo.xs;
}
PetscErrorCode PostStructured2d::addElmtVec(PetscInt rId,PetscScalar ***globalArray, PetscScalar **localArray,int mCpnt){
    PetscInt xI=0, yI=0;
    getElmtDmdaIndByRId(rId,&xI,&yI);
    addElmtVecByDmdaInd(xI,yI,globalArray,localArray,mCpnt);
    return 0;
}
PetscErrorCode PostStructured2d::addElmtVecByDmdaInd(PetscInt xI,PetscInt yI,PetscScalar ***globalArray, PetscScalar **localArray,int mCpnt){
    if(!globalArray){
        MessagePrinter::printErrorTxt("PostStructured2d: variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    if(yI+1>m_dmInfo.gys+m_dmInfo.gym-1||xI+1>m_dmInfo.gxs+m_dmInfo.gxm-1){
        snprintf(MessagePrinter::charBuff,MessagePrinter::buffLen,
                "DMDA (row,col)=(%d,%d) is out of range, max index = (%d,%d)",
                            xI+1,yI+1,m_dmInfo.gxs+m_dmInfo.gxm-1,m_dmInfo.gys+m_dmInfo.gym-1);
        MessagePrinter::printRankError(MessagePrinter::charBuff);      
        MessagePrinter::exitcfem();            
    }
    for(int cpntI=0;cpntI<mCpnt;cpntI++){
        globalArray[yI][xI][cpntI]+=localArray[0][cpntI];
        globalArray[yI][xI+1][cpntI]+=localArray[1][cpntI];
        globalArray[yI+1][xI+1][cpntI]+=localArray[2][cpntI];
        globalArray[yI+1][xI][cpntI]+=localArray[3][cpntI];
    }
    return 0;     
}