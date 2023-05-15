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
 * REFERENCE: Expression (4.42)*/
#include"MeshSystem/StructuredMesh2D.h"
#include"MathUtils/Vector2d.h"
#include"Utils/MessagePrinter.h"
#include"nlohmann/json.hpp"
#include"petsc.h"
#include<cmath>
#include<fstream>
#include "SolutionSystem/ArcLengthSolver.h"
StructuredMesh2D::StructuredMesh2D():
        m_array_nodes_coord0(nullptr),m_array_nodes_coord2(nullptr),
        m_array_nodes_uInc1(nullptr),m_array_nodes_uInc2(nullptr),
        m_array_nodes_u2(nullptr),
        m_array_nodes_residual1(nullptr),m_array_nodes_residual2(nullptr){
    m_dim=2;
    m_mDof_node=2;
    m_meshMode=MeshMode::STRUCTURED;
};
StructuredMesh2D::StructuredMesh2D(Timer *timerPtr):MeshSystem(timerPtr),
        m_array_nodes_coord0(nullptr),m_array_nodes_coord2(nullptr),
        m_array_nodes_uInc1(nullptr),m_array_nodes_uInc2(nullptr),
        m_array_nodes_u2(nullptr),
        m_array_nodes_residual1(nullptr),m_array_nodes_residual2(nullptr){
    m_dim=2;
    m_mDof_node=2;   
    m_meshMode=MeshMode::STRUCTURED; 
}
PetscErrorCode StructuredMesh2D::MeshSystemInit(MeshDescription *t_meshDesPtr){
    PetscErrorCode err=0;
    MeshMode meshMode=t_meshDesPtr->s_mode;
    MeshShape meshShape=t_meshDesPtr->s_shape;
    m_ifSaveMesh=t_meshDesPtr->s_ifSaveMesh;
    m_inputMeshFile_Name=t_meshDesPtr->s_inputMeshFile_Name;
    m_outputMeshFile_Name=t_meshDesPtr->s_outputMeshFile_Name;
    switch(meshMode){
        case MeshMode::STRUCTURED:
            err=initStructuredMesh(meshShape,t_meshDesPtr);
            break;
        case MeshMode::UNSTRUCTURED:
            break;
        default:
            break;
    }
    return err;
}

PetscErrorCode StructuredMesh2D::outputMeshFile(){
    // here the m_inputFileName must be "xxxx.json", which should has the .json extension
    std::ofstream meshout;
    if(m_rank==0){
        openMeshOutputFile(&meshout,std::ios::out);
    //****************************************
    //*** print out header information
    //****************************************
        meshout<<"<?xml version=\"1.0\"?>\n";
        meshout<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        meshout<<"<UnstructuredGrid>\n";
        meshout<<"<Piece NumberOfPoints=\""<<m_mNodes
            <<"\" NumberOfCells=\""<<m_mElmts<<"\">\n";

        meshout<<"<Points>\n";
        meshout<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";
        meshout.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    //*****************************
    // print out node coordinates
    //*****************************
    PetscScalar ***aCoord;
    PetscCall(DMDAVecGetArrayDOFRead(m_dm,m_nodes_coord0,&aCoord));
    for(int rankI=0;rankI<m_rankNum;rankI++){
        if(m_rank==rankI){// loop over all rank. print nodes coords if it's this rank's turn
            openMeshOutputFile(&meshout,std::ios::app);
            meshout<<std::scientific<<std::setprecision(6);
            PetscInt xs=m_daInfo.xs,    xe=m_daInfo.xs+m_daInfo.xm,
                     ys=m_daInfo.ys,    ye=m_daInfo.ys+m_daInfo.ym;
            for(PetscInt yI=ys;yI<ye;yI++){//loop over row in this rank
                for(PetscInt xI=xs;xI<xe;xI++){//loop over col in this rank
                    meshout<<aCoord[yI][xI][0]<<" ";
                    meshout<<aCoord[yI][xI][1]<<" ";
                    meshout<<0.0<<"\n";
                }
            }
            meshout.close();
        }
        PetscCall(PetscBarrier(NULL));  // ensure output coords in order
    }
    if(m_rank==0){
        openMeshOutputFile(&meshout,std::ios::app);
        meshout<<"</DataArray>\n";
        meshout<<"</Points>\n";
    //***************************************
    //*** For cell information
    //***************************************
        meshout<<"<Cells>\n";
        meshout<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        meshout.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    for(int rankI=0;rankI<m_rankNum;rankI++){
        if(m_rank==rankI){// loop over all rank. print element connectivity if it's this rank's turn
            openMeshOutputFile(&meshout,std::ios::app);
            meshout<<std::scientific<<std::setprecision(6);
            for(PetscInt eI=0;eI<m_mElmts_p;eI++){ // loop over elmts in this rank
                PetscInt mNodeInElmt=m_elmt_cnn[eI].size();
                for(PetscInt nI=0;nI<mNodeInElmt;nI++){// loop over nodes in this elmt
                    meshout<<m_elmt_cnn[eI][nI]<<" ";
                }
                meshout<<"\n";
            }
            meshout.close();
        }
        PetscCall(PetscBarrier(NULL));  // ensure output connectivity in order
    }
    if(m_rank==0){
        openMeshOutputFile(&meshout,std::ios::app);
        meshout<<"</DataArray>\n";
    //***************************************
    //*** for offset
    //***************************************
        meshout<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        meshout.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    const PetscInt mNodeInElmt=4;
    PetscInt offsetBase=((m_daInfo.mx-1)*m_daInfo.ys+m_daInfo.xs)*mNodeInElmt;
    for(int rankI=0;rankI<m_rankNum;rankI++){
        if(m_rank==rankI){// loop over all rank. print element connectivity if it's this rank's turn
            openMeshOutputFile(&meshout,std::ios::app);
            for(PetscInt eI=0;eI<m_mElmts_p;eI++){ // loop over elmts in this rank
                offsetBase+=mNodeInElmt;
                meshout<<offsetBase<<"\n";
            }
            
            meshout.close();
        }
        PetscCall(PetscBarrier(NULL));  // ensure output offset in order
    }
    if(m_rank==0){
        openMeshOutputFile(&meshout,std::ios::app);
        meshout<<"</DataArray>\n";
    //***************************************
    //*** for VTKCellType
    //***************************************
        meshout<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(PetscInt eI=0;eI<m_mElmts;eI++){ // loop over elmts in all rank
            meshout<<this->vtkType<<"\n";
        }
        meshout<<"</DataArray>\n";
        meshout<<"</Cells>\n";
        meshout<<"</Piece>\n";
        meshout<<"</UnstructuredGrid>\n";
        meshout<<"</VTKFile>"<<endl;
        meshout.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    PetscCall(DMDAVecRestoreArrayDOFRead(m_dm,m_nodes_coord0,&aCoord));
    return 0;
}

PetscErrorCode StructuredMesh2D::openNodeVariableVec(NodeVariableType vType, Vec *variableVecPtr, int state, VecAccessMode mode){
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(vType,state);
    Vec &localVec=getNodeLocalVariableVecRef(vType,state);
    if(arrayPtrRef){
        MessagePrinter::printErrorTxt("Node variable array need to be restored before setting up");
        MessagePrinter::exitcfem();             
    }
    else{
        PetscCall(DMGetLocalVector(m_dm,&localVec));
        if(mode==VecAccessMode::READ){
            PetscCall(DMDAVecGetArrayDOFRead(m_dm,localVec,&arrayPtrRef));
            PetscCall(DMGlobalToLocal(m_dm,*variableVecPtr,INSERT_VALUES,localVec));
        }
        else if(mode==VecAccessMode::WRITE){
            PetscCall(DMDAVecGetArrayDOFWrite(m_dm,localVec,&arrayPtrRef));
            PetscCall(VecZeroEntries(localVec));  
        }
    }
    return 0;
}

PetscErrorCode StructuredMesh2D::closeNodeVariableVec(NodeVariableType vType, Vec *variableVecPtr, int state, VecAccessMode mode){
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(vType,state);
    Vec &localVec=getNodeLocalVariableVecRef(vType,state);;
    if(arrayPtrRef){
        if(mode==VecAccessMode::READ){
            PetscCall(DMDAVecRestoreArrayDOFRead(m_dm,localVec,&arrayPtrRef));
            PetscCall(DMRestoreLocalVector(m_dm,&localVec));
        }
        else if(mode==VecAccessMode::WRITE){
            PetscCall(DMDAVecRestoreArrayDOFWrite(m_dm,localVec,&arrayPtrRef));
            PetscCall(VecZeroEntries(*variableVecPtr));  
            PetscCall(DMLocalToGlobal(m_dm,localVec,ADD_VALUES,*variableVecPtr));
            PetscCall(DMRestoreLocalVector(m_dm,&localVec));
        }  
    }
    else{
        MessagePrinter::printErrorTxt("array is point to null, can not be restored");
        MessagePrinter::exitcfem();
    }
    arrayPtrRef=nullptr;
    return 0;
}

PetscErrorCode StructuredMesh2D::getNodeCoord(PetscInt nodeRId,int state,Vector *coordsPtr){
    checkElmtRId(nodeRId);
    PetscInt xI=0,yI=0;
    getNodeDmdaIndByRId(nodeRId,&xI,&yI);
    getNodeVariableByDmdaInd(NodeVariableType::COORD,xI,yI,state,coordsPtr);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeCoord(PetscInt elmtRId,int state,Vector *coordsPtr,PetscInt *nodeNum){
    checkElmtRId(elmtRId);
    PetscInt xI=0,yI=0;
    getElmtDmdaIndByRId(elmtRId,&xI,&yI);
    getElmtNodeVariableByDmdaInd(NodeVariableType::COORD,xI,yI,state,coordsPtr,nodeNum);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeUInc(PetscInt elmtRId,int state,Vector *uIncPtr,PetscInt *nodeNum){
    checkElmtRId(elmtRId);
    PetscInt xI=0,yI=0;
    getElmtDmdaIndByRId(elmtRId,&xI,&yI);
    getElmtNodeVariableByDmdaInd(NodeVariableType::UINC,xI,yI,state,uIncPtr,nodeNum);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeResidual(PetscInt elmtRId,int state,Vector *residualPtr,PetscInt *nodeNum){
    checkElmtRId(elmtRId);
    PetscInt xI=0,yI=0;
    getElmtDmdaIndByRId(elmtRId,&xI,&yI);
    getElmtNodeVariableByDmdaInd(NodeVariableType::RESIDUAL,xI,yI,state,residualPtr,nodeNum);
    return 0;    
}
PetscErrorCode StructuredMesh2D::updateConfig(SNES *sensPtr){
    PetscCall(SNESGetSolution(*sensPtr,&m_nodes_uInc2));
    PetscCall(VecAYPX(m_nodes_coord2,1.0,m_nodes_uInc2));
    PetscCall(VecAYPX(m_nodes_u2,1.0,m_nodes_uInc2));
    PetscCall(VecZeroEntries(m_node_residual2));
    return 0;
}
PetscErrorCode StructuredMesh2D::updateConfig(void *solver,AlgorithmType algo){
    switch(algo){
        case AlgorithmType::STANDARD:{
            PetscCall(SNESGetSolution(*(SNES *)solver,&m_nodes_uInc2));
            break;
        }
        case AlgorithmType::ARCLENGTH_CYLENDER:{
            ArcLengthSolver *arcSolver=(ArcLengthSolver *)solver;
            arcSolver->getSolution(&m_nodes_uInc2);
        }
    }
    PetscCall(VecAYPX(m_nodes_coord2,1.0,m_nodes_uInc2));
    PetscCall(VecAYPX(m_nodes_u2,1.0,m_nodes_uInc2));
    PetscCall(VecZeroEntries(m_node_residual2));
    return 0;
}
PetscErrorCode StructuredMesh2D::addElmtAMatrix(PetscInt rid,MatrixXd *matrixPtr,Mat *APtr){
    checkElmtRId(rid);
    PetscInt xI=0,yI=0;
    getElmtDmdaIndByRId(rid,&xI,&yI);
    addElmtAMatrixByDmdaInd(xI,yI,matrixPtr,APtr);
    return 0;
}

PetscErrorCode StructuredMesh2D::addElmtResidual(PetscInt rid,Vector *residualPtr, Vec *fPtr){
    checkElmtRId(rid);
    PetscInt xI=0,yI=0;
    getElmtDmdaIndByRId(rid,&xI,&yI);
    addElmtResidualByDmdaInd(xI,yI,residualPtr,fPtr);
    return 0;
}
PetscErrorCode StructuredMesh2D::initStructuredMesh(MeshShape t_meshShape,MeshDescription *t_meshDesPtr){
    m_timerPtr->startTimer();
    MessagePrinter::printNormalTxt("Start to init the sturcted 2D mesh system");
    switch (t_meshShape)
    {
    case MeshShape::RECTANGULAR:
        m_geoParam[0]=t_meshDesPtr->s_size_json.at("xmax");
        m_geoParam[1]=t_meshDesPtr->s_size_json.at("ymax");
        break;
    case MeshShape::SIN:
    case MeshShape::HALFSIN:
    case MeshShape::COS:
    case MeshShape::HALFCOS:
        m_geoParam[0]=t_meshDesPtr->s_size_json.at("span");
        m_geoParam[1]=t_meshDesPtr->s_size_json.at("amplitude");
        m_geoParam[2]=t_meshDesPtr->s_size_json.at("width");
        break;
    case MeshShape::HALFCOSPLUSSTEP:
        m_geoParam[0]=t_meshDesPtr->s_size_json.at("span");
        m_geoParam[1]=t_meshDesPtr->s_size_json.at("amplitude");
        m_geoParam[2]=t_meshDesPtr->s_size_json.at("width");
        m_geoParam[3]=t_meshDesPtr->s_size_json.at("step");
        break;
    default:
        m_geoParam[0]=t_meshDesPtr->s_size_json.at("xmax");
        m_geoParam[1]=t_meshDesPtr->s_size_json.at("ymax");
        break;
    }
    PetscInt nx=t_meshDesPtr->s_nx,ny=t_meshDesPtr->s_ny;
    PetscInt *localNy=new PetscInt[m_rankNum];
    PetscInt rankave=ny/m_rankNum;          /**< nodes' num in y direction per processor*/
    PetscInt rankrem=ny-rankave*m_rankNum;  /**< num of the unassigned nodes, assign them to the last m_rank*/
    for(int i=0;i<m_rankNum;i++){
        localNy[i]=rankave;
    }
    localNy[m_rankNum-1]+=rankrem;
    /** crate DMDA*************************************************************/
    PetscCall(DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
            DMDA_STENCIL_BOX,           /** stencil type. Use either DMDA_STENCIL_BOX or DMDA_STENCIL_STAR.*/
            nx,ny,                      /** global dimension in each direction of the array*/
            1,m_rankNum,                /** corresponding number of processors in each dimension (or PETSC_DECIDE to have calculated)*/
            m_mDof_node,                /** number of degrees of freedom per node*/
            1,                          /** stencil width*/
            &nx,localNy,                /** arrays containing the number of nodes in each cell along the x and y coordinates, or NULL.
                                         *  If non-null, these must be of length as m and n, and the corresponding m and n cannot be 
                                         *  PETSC_DECIDE. The sum of the lx[] entries must be M, and the sum of the ly[] entries must be N.*/
            &m_dm                       /** DM pointer*/
            ));
    PetscCall(DMSetFromOptions(m_dm));
    PetscCall(DMSetUp(m_dm));
    /**************************************************************************/
    /** Create Vec of nodes' variable parallel data managed by DMDA)***********/
    /**************************************************************************/
    PetscCall(DMCreateGlobalVector(m_dm,&m_nodes_coord0));
    PetscCall(DMCreateGlobalVector(m_dm,&m_nodes_coord2));
    PetscCall(DMCreateGlobalVector(m_dm,&m_nodes_u2));
    PetscCall(DMCreateGlobalVector(m_dm,&m_nodes_uInc2));
    PetscCall(DMCreateGlobalVector(m_dm,&m_node_residual2));
    PetscCall(DMCreateGlobalVector(m_dm,&m_node_load));
    PetscCall(VecZeroEntries(m_nodes_coord0));
    PetscCall(VecZeroEntries(m_nodes_coord2));
    PetscCall(VecZeroEntries(m_nodes_u2));
    PetscCall(VecZeroEntries(m_nodes_uInc2));
    PetscCall(VecZeroEntries(m_node_residual2));
    PetscCall(VecZeroEntries(m_node_load));
    /**************************************************************************/
    /** Create AMatrix managed by DMDA)****************************************/
    /**************************************************************************/
    PetscCall(DMCreateMatrix(m_dm,&m_AMatrix2));
    PetscCall(MatSetFromOptions(m_AMatrix2));
    PetscCall(MatZeroEntries(m_AMatrix2));
    /**************************************************************************/
    /** cal node and element num***********************************************/
    /**************************************************************************/
    PetscCall(DMDAGetLocalInfo(m_dm,&m_daInfo));
    m_mNodes=nx*ny;
    m_mNodes_p=nx*localNy[m_rank];
    m_mElmts=(nx-1)*(ny-1);
    PetscInt nx_elmt,ny_elmt;   /**< element num in x,y direction of this m_rank*/
    nx_elmt=nx-1;
    if(m_rank==m_rankNum-1)ny_elmt=localNy[m_rank]-1;
    else ny_elmt=localNy[m_rank];
    m_mElmts_p=nx_elmt*ny_elmt;
    /**************************************************************************/
    /** set node and elemnt's global id, and element's connectivity************/
    /**************************************************************************/
    m_node_gId.resize(m_mNodes_p);
    m_dof_gId.resize(m_mNodes_p*m_mDof_node);
    m_elmt_gId.resize(m_mElmts_p);
    m_elmt_cnn.resize(m_mElmts_p);
    const int mNodePElmt=4; /**< node num per element*/
    PetscInt dofRId=0,nodeRId=0,elmtRId=0;                  /**< current ndoe/element's m_rank id*/
    PetscInt nodeGId=m_daInfo.mx*m_daInfo.ys;               /**< current node's global id*/
    PetscInt dofGId=m_daInfo.mx*m_daInfo.ys*m_mDof_node;    /**< current dof's global id*/
    PetscInt elmtGId=(m_daInfo.mx-1)*m_daInfo.ys;           /**< current node's global id*/
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\033[1;33m"));// set color to yellow
    PetscCall(PetscSynchronizedPrintf(              /**< print the node and element's global id range of every m_rank*/
            PETSC_COMM_WORLD,
            "[%2d]: owned node global id: %d -> %d \n      owned elment global id: %d -> %d\n",
            m_rank, nodeGId, nodeGId+m_mNodes_p, elmtGId, elmtGId+m_mElmts_p
    ));
    PetscCall(PetscSynchronizedPrintf(              /**< print the node and element's DMDA index range of every m_rank*/
            PETSC_COMM_WORLD,
            "[%2d]: owned node DMDA index: [%d, %d] X [%d, %d] \n      owned elment DMDA index: [%d, %d] X [%d, %d]\n",
            m_rank, m_daInfo.xs, m_daInfo.ys, m_daInfo.xs+m_daInfo.xm, m_daInfo.ys+m_daInfo.ym,                     // node DMDA index range
                    m_daInfo.xs, m_daInfo.ys, m_daInfo.xs+m_daInfo.xm-1, m_daInfo.ys+m_daInfo.ym-(m_rank==m_rankNum-1)
    ));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\033[0m"));// recover color
    PetscScalar ***aCoord0;         /** ptr to the node's coords in ref config*/
    Vector2d CoordLocal;          /** store coord0 of a node*/
    PetscCall(DMDAVecGetArrayDOF(m_dm,m_nodes_coord0,&aCoord0));
    /**************************************************************************/
    /** create node set "left", "right" vector ********************************/
    /**************************************************************************/
    vector<PetscInt> leftSet;
    vector<PetscInt> rightSet;
    vector<PetscInt> bottomSet;
    vector<PetscInt> topSet;
    vector<PetscInt> stepSet;
    for(int yI=m_daInfo.ys;yI<m_daInfo.ys+m_daInfo.ym;yI++){    // loop over every row
        for(int xI=m_daInfo.xs;xI<m_daInfo.xs+m_daInfo.xm;xI++){// loop over every col
            for(int dI=0;dI<m_mDof_node;++dI){// loop over dof in a node
                m_dof_gId[dofRId]=dofGId;
                ++dofRId;   ++dofGId;
            }
            m_node_gId[nodeRId]=nodeGId;                        // set node's global id
            if(yI==0) bottomSet.push_back(nodeRId);             // append node to set bottom
            if(yI==m_daInfo.my-1) topSet.push_back(nodeRId);    // append node to set top
            if(xI==0) leftSet.push_back(nodeRId);               // append node to set left
            if(xI==m_daInfo.mx-1) rightSet.push_back(nodeRId);  // append node to set right
            if(yI<m_daInfo.my-1&&xI<m_daInfo.mx-1){
                m_elmt_gId[elmtRId]=elmtGId;    // set element's global id
                /** Set element's connectivity*/
                m_elmt_cnn[elmtRId].resize(mNodePElmt);
                m_elmt_cnn[elmtRId][0]=nodeGId;
                m_elmt_cnn[elmtRId][1]=nodeGId+1;
                m_elmt_cnn[elmtRId][2]=nodeGId+m_daInfo.mx+1;
                m_elmt_cnn[elmtRId][3]=nodeGId+m_daInfo.mx;
                ++elmtRId;
                ++elmtGId;
            }
            /** set the nodes' coords in ref config*/
            switch (t_meshShape)
            {
            case MeshShape::RECTANGULAR:
                getRectCoord0ByDmdaInd(xI,yI,&CoordLocal);   
                break;
            case MeshShape::SIN:
                getSinCoord0ByDmdaInd(xI,yI,&CoordLocal);
                break;
            case MeshShape::HALFSIN:
                getHalfSinCoord0ByDmdaInd(xI,yI,&CoordLocal);
                break;
            case MeshShape::COS:
                getCosCoord0ByDmdaInd(xI,yI,&CoordLocal);
                break;
            case MeshShape::HALFCOS:
                getHalfCosCoord0ByDmdaInd(xI,yI,&CoordLocal);
                break;   
            case MeshShape::HALFCOSPLUSSTEP:
                getHalfCosPlusStepCoord0ByDmdaInd(xI,yI,&CoordLocal);
                if(CoordLocal(0)==m_geoParam[1]/2.0+m_geoParam[2]){// the node belong to the step edge
                    stepSet.push_back(nodeRId);
                }
                break;                       
            default:
                MessagePrinter::printErrorTxt("unsupported mesh shape.");
                MessagePrinter::exitcfem();
                break;
            }
            for(int dofI=0;dofI<m_mDof_node;dofI++){
                aCoord0[yI][xI][dofI]=CoordLocal(dofI);
            }
            ++nodeRId;
            ++nodeGId;
        }
    }
    // ************************************************************************/
    // create set "all"                                                     ***/
    // ************************************************************************/
    vector<PetscInt> elmtAll;
    vector<PetscInt> nodeAll;
    elmtAll.resize(m_mElmts_p);
    nodeAll.resize(m_mNodes_p);
    PetscInt rid=0;
    for(vector<PetscInt>::iterator it=nodeAll.begin();it!=nodeAll.end();++it){
        *it=rid++;
    }
    rid=0;
    for(vector<PetscInt>::iterator it=elmtAll.begin();it!=elmtAll.end();++it){
        *it=rid++;
    }
    m_setManager.createSet("all",SetType::ELEMENT,&elmtAll);
    m_setManager.createSet("all",SetType::NODE,&nodeAll);
    m_setManager.createSet("left",SetType::NODE,&leftSet);
    m_setManager.createSet("right",SetType::NODE,&rightSet);
    m_setManager.createSet("top",SetType::NODE,&topSet);
    m_setManager.createSet("bottom",SetType::NODE,&bottomSet);
    switch(t_meshShape){
        case MeshShape::HALFCOSPLUSSTEP:
            m_setManager.createSet("step",SetType::NODE,&stepSet);
            break;
        default:
            break;
    }
    PetscCall(DMDAVecRestoreArrayDOF(m_dm,m_nodes_coord0,&aCoord0));
    PetscCall(VecAssemblyBegin(m_nodes_coord0));
    PetscCall(VecAssemblyEnd(m_nodes_coord0));
    /**************************************************************************/
    /**copy coord0 to coord1, coord2                                        ***/
    /**************************************************************************/
    PetscCall(VecCopy(m_nodes_coord0,m_nodes_coord2));
    delete[] localNy; 
    m_timerPtr->endTimer();
    m_timerPtr->printElapseTime("Mesh system init is done",false);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeVariableByDmdaInd(NodeVariableType vType ,PetscInt xI,PetscInt yI,int state,Vector *variablePtr,PetscInt *nodeNum){
    if(nodeNum) *nodeNum=4;
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(vType,state);
    Vector2d *Vector2dPtr =(Vector2d *)variablePtr;
    if(!arrayPtrRef){
        MessagePrinter::printErrorTxt("variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    for(int dofI=0;dofI<m_mDof_node;dofI++){
        Vector2dPtr[0](dofI)=arrayPtrRef[yI][xI][dofI];
        Vector2dPtr[1](dofI)=arrayPtrRef[yI][xI+1][dofI];
        Vector2dPtr[2](dofI)=arrayPtrRef[yI+1][xI+1][dofI];
        Vector2dPtr[3](dofI)=arrayPtrRef[yI+1][xI][dofI];
    }
    return 0;
}

PetscErrorCode StructuredMesh2D::getNodeVariableByDmdaInd(NodeVariableType vType, PetscInt xI,PetscInt yI,int state,Vector *variablePtr){
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(vType,state);
    Vector2d *Vector2dPtr =(Vector2d *)variablePtr;
    if(!arrayPtrRef){
        MessagePrinter::printErrorTxt("variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    for(int dofI=0;dofI<m_mDof_node;dofI++){
        (*Vector2dPtr)(dofI)=arrayPtrRef[yI][xI][dofI];
    }
    return 0;   
}

PetscErrorCode StructuredMesh2D::addElmtAMatrixByDmdaInd(PetscInt xI,PetscInt yI,MatrixXd *matrixPtr,Mat *APtr){
    const int mNodePerElmt=m_dim*m_dim;             /**< node num per elmt*/
    const int mDofPerElmt=mNodePerElmt*m_mDof_node; /**< dof num per elmt*/
    MatStencil row[mDofPerElmt], col[mDofPerElmt];  /**< row & col stencil for elmt's K*/
    PetscScalar entry[mDofPerElmt*mDofPerElmt];     /**< array storing elmt's K's entries*/
    static const int relPositon[4][2]=              /**< (node id in elmt, direction of DMDA) -> relative positon to elmt in specific direction*/
    {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    PetscInt rowDofI, colDofI;                      /**< tmp index for loop over Ke's elmt dof in row/col*/
    PetscInt nodeDofI, nodeDofI2;                   /**< tmp index for loop over dof of a node*/
    if(yI+1>m_daInfo.gys+m_daInfo.gym-1||xI+1>m_daInfo.gxs+m_daInfo.gxm-1){
        snprintf(MessagePrinter::charBuff,MessagePrinter::buffLen,
                "DMDA (row,col)=(%d,%d) is out of range, max index = (%d,%d)",
                            xI+1,yI+1,m_daInfo.gxs+m_daInfo.gxm-1,m_daInfo.gys+m_daInfo.gym-1);
        MessagePrinter::printRankError(MessagePrinter::charBuff);      
        MessagePrinter::exitcfem();            
    }
    int rowNodeI = 0; /**< node id of Ke's row in a elmt*/
    int colNodeI = 0; /**< node id of Ke's col in a elmt*/
    rowDofI=0;
    for(rowNodeI=0;rowNodeI<mNodePerElmt;++rowNodeI){
        PetscInt rowI, rowJ;    /**< tmp index for loop over Ke's node's row*/
        rowI=xI+relPositon[rowNodeI][0];  rowJ=yI+relPositon[rowNodeI][1];
        for(nodeDofI=0;nodeDofI<m_mDof_node;++nodeDofI){
            row[rowDofI].i=rowI;        col[rowDofI].i=rowI;       
            row[rowDofI].j=rowJ;        col[rowDofI].j=rowJ;    
            row[rowDofI].c=nodeDofI;    col[rowDofI].c=nodeDofI;
            colDofI=0;
            for(colNodeI=0;colNodeI<mNodePerElmt;++colNodeI){
                for(nodeDofI2=0;nodeDofI2<m_mDof_node;++nodeDofI2){
                    entry[mDofPerElmt*rowDofI+colDofI]=(*matrixPtr)(rowDofI,colDofI);
                    ++colDofI;
                }
            }
            ++rowDofI;      
        }
    }
    PetscCall(MatSetValuesStencil(*APtr,mDofPerElmt,row,mDofPerElmt,col,entry,ADD_VALUES));
    return 0;
}

PetscErrorCode StructuredMesh2D::addElmtResidualByDmdaInd(PetscInt xI,PetscInt yI,Vector *residualPtr,Vec *fPtr){
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(NodeVariableType::RESIDUAL,1);
    Vector2d *residualPtr2d =(Vector2d *)residualPtr;
    if(fPtr){}
    if(!arrayPtrRef){
        MessagePrinter::printErrorTxt("variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    if(yI+1>m_daInfo.gys+m_daInfo.gym-1||xI+1>m_daInfo.gxs+m_daInfo.gxm-1){
        snprintf(MessagePrinter::charBuff,MessagePrinter::buffLen,
                "DMDA (row,col)=(%d,%d) is out of range, max index = (%d,%d)",
                            xI+1,yI+1,m_daInfo.gxs+m_daInfo.gxm-1,m_daInfo.gys+m_daInfo.gym-1);
        MessagePrinter::printRankError(MessagePrinter::charBuff);      
        MessagePrinter::exitcfem();            
    }
    for(int dofI=0;dofI<m_mDof_node;dofI++){
        arrayPtrRef[yI][xI][dofI]+=residualPtr2d[0](dofI);
        arrayPtrRef[yI][xI+1][dofI]+=residualPtr2d[1](dofI);
        arrayPtrRef[yI+1][xI+1][dofI]+=residualPtr2d[2](dofI);
        arrayPtrRef[yI+1][xI][dofI]+=residualPtr2d[3](dofI);
    }
    return 0;        
}

void StructuredMesh2D::openMeshOutputFile(ofstream *ofPtr,ios_base::openmode mode){
    char buff[110];
    string str;
    string _MeshFileName;
    if(m_outputMeshFile_Name.size()<2){
        _MeshFileName="mesh.vtu";
    }
    else{
        _MeshFileName=m_outputMeshFile_Name;
    }
    ofPtr->open(_MeshFileName,mode);
    if(!ofPtr->is_open()){
        snprintf(buff,110,"can\'t write mesh to vtu file(=%28s), please make sure you have the write permission",_MeshFileName.c_str());
        str=buff;
        MessagePrinter::printErrorTxt(str);
        MessagePrinter::exitcfem();
    }    
}

void StructuredMesh2D::getRectCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span, amplitude, width in order)*/
    (*aCoord)(0)=xI*m_geoParam[0]/(m_daInfo.mx-1);
    (*aCoord)(1)=yI*m_geoParam[1]/(m_daInfo.my-1);
}

void StructuredMesh2D::getSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span l, amplitude a, width w in order)*/
    (*aCoord)(1)=yI*m_geoParam[0]/(m_daInfo.my-1);
    /** x=a*sin(2pi*y/l)*/
    (*aCoord)(0)=m_geoParam[1]/2*sin(2.0*M_PI*(*aCoord)(1)/m_geoParam[0])+xI*m_geoParam[2]/(m_daInfo.mx-1);
}

void StructuredMesh2D::getHalfSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span l, amplitude a, width w in order)*/
    (*aCoord)(1)=yI*m_geoParam[0]/(m_daInfo.my-1);
    /** x=a*sin(pi*y/l)*/
    (*aCoord)(0)=m_geoParam[1]/2*sin(M_PI*(*aCoord)(1)/m_geoParam[0])+xI*m_geoParam[2]/(m_daInfo.mx-1);  
}
void StructuredMesh2D::getCosCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span l, amplitude a, width w in order)*/
    (*aCoord)(1)=-0.5*m_geoParam[0]+yI*m_geoParam[0]/(m_daInfo.my-1);
    /** x=a*cos(2pi*y/l)*/
    (*aCoord)(0)=m_geoParam[1]/2*cos(2.0*M_PI*(*aCoord)(1)/m_geoParam[0])+xI*m_geoParam[2]/(m_daInfo.mx-1);
}

void StructuredMesh2D::getHalfCosCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span l, amplitude a, width w in order)*/
    (*aCoord)(1)=-m_geoParam[0]+yI*m_geoParam[0]/(m_daInfo.my-1);
    /** x=a*cos(pi*y/l)*/
    (*aCoord)(0)=m_geoParam[1]/2*cos(M_PI*(*aCoord)(1)/m_geoParam[0])+xI*m_geoParam[2]/(m_daInfo.mx-1); 
}
void StructuredMesh2D::getHalfCosPlusStepCoord0ByDmdaInd(PetscInt xI,PetscInt yI, Vector *aCoord){
    PetscInt my_cos=(m_daInfo.my-1)*m_geoParam[0]/(m_geoParam[0]+m_geoParam[3])+1;
    PetscInt my_step=m_daInfo.my-my_cos;
    if(yI<my_cos){
        (*aCoord)(1)=-m_geoParam[0]+yI*m_geoParam[0]/(my_cos-1);
        /** x=a*cos(pi*y/l)*/
        (*aCoord)(0)=m_geoParam[1]/2*cos(M_PI*(*aCoord)(1)/m_geoParam[0])+xI*m_geoParam[2]/(m_daInfo.mx-1);     
        (*aCoord)(1)-=m_geoParam[3]; 
    }
    else{
        (*aCoord)(1)=-m_geoParam[3]+ (yI+1-my_cos)*m_geoParam[3]/my_step;
        (*aCoord)(0)=m_geoParam[1]/2+xI*m_geoParam[2]/(m_daInfo.mx-1);
    }
}
void StructuredMesh2D::getElmtDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr){
    *yIPtr=rId/(m_daInfo.mx-1)+m_daInfo.ys;
    *xIPtr=rId%(m_daInfo.mx-1)+m_daInfo.xs;
}

void StructuredMesh2D::getNodeDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr){
    *yIPtr=rId/m_daInfo.mx+m_daInfo.ys;
    *xIPtr=rId%m_daInfo.mx+m_daInfo.xs;
}

Vec & StructuredMesh2D::getNodeLocalVariableVecRef(NodeVariableType vType, int state){
    switch (vType)
    {
    case NodeVariableType::COORD:
        switch (state)
        {
        case 0:
            return m_nodes_coord0_local;
            break;
        case 2:
            return m_nodes_coord2_local;
            break;       
        default:
            MessagePrinter::printErrorTxt("required Coords vec state must be 0 or 2.");
            MessagePrinter::exitcfem();
            break;
        }
        break;
    case NodeVariableType::U:
        switch (state)
        {
        case 2:
            return m_nodes_u2_local;
            break;
        default:
            MessagePrinter::printErrorTxt("required u vec state must be 2.");
            MessagePrinter::exitcfem();        
            break;
        }
        break;
    case NodeVariableType::UINC:
        switch (state)
        {
        case 1:
            return m_nodes_uInc1_local;
            break;    
        case 2:
            return m_nodes_uInc2_local;
            break;    
        default:
            MessagePrinter::printErrorTxt("required incremental U vec state must be 1 or 2.");
            MessagePrinter::exitcfem();
            break;
        }
        break;      
    case NodeVariableType::RESIDUAL:
        if (state==1){
           return m_node_residual1_local;
        }
        if (state==2){
           return m_node_residual2_local;
        }
        else{
            MessagePrinter::printErrorTxt("required residual vec state must be 1 or 2.");
            MessagePrinter::exitcfem();           
        }
        break;
    case NodeVariableType::LOAD:
        if (state==1){
            return m_node_load_local;
        }
        else{
            MessagePrinter::printErrorTxt("required load vec state must be 1.");
            MessagePrinter::exitcfem();                
        }
        break;
    default:
        MessagePrinter::printErrorTxt("required node variable vec are not unsupported in current mesh system");
        MessagePrinter::exitcfem();       
        break;
    }
    return m_nodes_coord0_local;    
}

PetscScalar *** & StructuredMesh2D::getNodeVariablePtrRef(NodeVariableType vType, int state){
    
    switch (vType)
    {
    case NodeVariableType::COORD:
        switch (state)
        {
        case 0:
            return m_array_nodes_coord0;
            break;
        case 2:
            return m_array_nodes_coord2;
            break;       
        default:
            MessagePrinter::printErrorTxt("required Coords array address state must be 0 or 2.");
            MessagePrinter::exitcfem();
            break;
        }
        break;
    case NodeVariableType::U:
        if(state==2){
            return m_array_nodes_u2;
        }
        else{
            MessagePrinter::printErrorTxt("required u array address state must be 2.");
            MessagePrinter::exitcfem();            
        }
        break;
    case NodeVariableType::UINC:
        switch (state)
        {
        case 1:
            return m_array_nodes_uInc1;
            break;     
        case 2:
            return m_array_nodes_uInc2;
            break;
        default:
            MessagePrinter::printErrorTxt("required incremental U array address state must be 1 or 2.");
            MessagePrinter::exitcfem();
            break;
        }
        break;      
    case NodeVariableType::RESIDUAL:
        switch (state)
        {
        case 1:
            return m_array_nodes_residual1;
            break;     
        case 2:
            return m_array_nodes_residual2;
            break;
        default:
            MessagePrinter::printErrorTxt("required residual array address state must be 1 or 2.");
            MessagePrinter::exitcfem();    
            break;
        }
        break;
    case NodeVariableType::LOAD:
        if (state==1){
            return m_array_nodes_load;
        }
        else{
            MessagePrinter::printErrorTxt("required load array address state must be 1.");
            MessagePrinter::exitcfem();                
        }
        break;
    default:
        MessagePrinter::printErrorTxt("required node variable array address are not unsupported in current mesh system");
        MessagePrinter::exitcfem();       
        break;
    }
    return m_array_nodes_coord0;      
}
PetscErrorCode StructuredMesh2D::createGlobalVec(Vec *t_vecAdr){
    PetscCall(DMCreateGlobalVector(m_dm,t_vecAdr));
    PetscCall(VecZeroEntries(*t_vecAdr));
    return 0;
}
PetscErrorCode StructuredMesh2D::destroyGlobalVec(Vec *t_vecAdr){
    VecDestroy(t_vecAdr);
    return 0;
}
int StructuredMesh2D::nodeGId2RId(int gId){
    PetscInt startGId=m_daInfo.mx*m_daInfo.ys;       /**< node's start global id in this rank*/
    return gId-startGId;
}
int StructuredMesh2D::elmtGId2RId(int gId){
    PetscInt startGId=(m_daInfo.mx-1)*m_daInfo.ys;       /**< node's start global id in this rank*/
    return gId-startGId;    
}
PetscErrorCode StructuredMesh2D::printVaribale(NodeVariableType vType, Vec *variableVecPtr, int state, int comp){
    openNodeVariableVec(vType, variableVecPtr, state,VecAccessMode::READ);
    PetscScalar ***array=getNodeVariablePtrRef(vType,state);
    for(PetscMPIInt rankI=0;rankI<m_rankNum;++rankI){
        if(m_rank==rankI){
            if(m_rank==0){
                MessagePrinter::printStarsRank();
                printf("      ");
                for(PetscInt xI=m_daInfo.xs;xI<m_daInfo.xs+m_daInfo.xm;++xI){
                    printf("%12d ",xI);
                }
                printf("\n");
            }
            for(PetscInt yI=m_daInfo.ys;yI<m_daInfo.ys+m_daInfo.ym;++yI){
                printf("%4d: ",yI);
                for(PetscInt xI=m_daInfo.xs;xI<m_daInfo.xs+m_daInfo.xm;++xI){
                    printf("%12.5e ",array[yI][xI][comp]);
                }
                printf("\n");
            }
        }
        PetscCall(PetscBarrier(NULL));
    }
    closeNodeVariableVec(vType, variableVecPtr, state,VecAccessMode::READ);
    return 0;
}