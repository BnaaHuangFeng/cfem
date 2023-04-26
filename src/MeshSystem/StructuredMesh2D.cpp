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
    PetscCall(DMGetLocalVector(m_dm,&localVec));
    PetscCall(DMGlobalToLocal(m_dm,*variableVecPtr,INSERT_VALUES,localVec));
    if(arrayPtrRef){
        MessagePrinter::printErrorTxt("Node variable array need to be restored before setting up");
        MessagePrinter::exitcfem();             
    }
    else{
        if(mode==VecAccessMode::READ)
            PetscCall(DMDAVecGetArrayDOFRead(m_dm,localVec,&arrayPtrRef));
        else if(mode==VecAccessMode::WRITE)
            PetscCall(DMDAVecGetArrayDOFWrite(m_dm,localVec,&arrayPtrRef));
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
    PetscInt *xI=nullptr,*yI=nullptr;
    getNodeDmdaIndByRId(nodeRId,xI,yI);
    getNodeVariableByDmdaInd(NodeVariableType::COORD,*xI,*yI,state,coordsPtr);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeCoord(PetscInt elmtRId,int state,Vector *coordsPtr,PetscInt *nodeNum){
    checkElmtRId(elmtRId);
    PetscInt *xIPtr=nullptr, *yIPtr=nullptr;
    getElmtDmdaIndByRId(elmtRId,xIPtr,yIPtr);
    getElmtNodeVariableByDmdaInd(NodeVariableType::COORD,*xIPtr,*yIPtr,state,coordsPtr,nodeNum);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeUInc(PetscInt elmtRId,int state,Vector *uIncPtr,PetscInt *nodeNum){
    checkElmtRId(elmtRId);
    PetscInt *xIPtr=nullptr, *yIPtr=nullptr;
    getElmtDmdaIndByRId(elmtRId,xIPtr,yIPtr);
    getElmtNodeVariableByDmdaInd(NodeVariableType::UINC,*xIPtr,*yIPtr,state,uIncPtr,nodeNum);
    return 0;
}

PetscErrorCode StructuredMesh2D::getElmtNodeResidual(PetscInt elmtRId,int state,Vector *residualPtr,PetscInt *nodeNum){
    checkElmtRId(elmtRId);
    PetscInt *xIPtr=nullptr, *yIPtr=nullptr;
    getElmtDmdaIndByRId(elmtRId,xIPtr,yIPtr);
    getElmtNodeVariableByDmdaInd(NodeVariableType::RESIDUAL,*xIPtr,*yIPtr,state,residualPtr,nodeNum);
    return 0;    
}

PetscErrorCode StructuredMesh2D::addElmtAMatrix(PetscInt rid,MatrixXd *matrixPtr,Mat *APtr){
    checkElmtRId(rid);
    PetscInt *xIPtr=nullptr, *yIPtr=nullptr;
    getElmtDmdaIndByRId(rid,xIPtr,yIPtr);
    addElmtAMatrixByDmdaInd(*xIPtr,*yIPtr,matrixPtr,APtr);
    return 0;
}

PetscErrorCode StructuredMesh2D::addElmtResidual(PetscInt rid,Vector *residualPtr, Vec *fPtr){
    checkElmtRId(rid);
    PetscInt *xIPtr=nullptr, *yIPtr=nullptr;
    getElmtDmdaIndByRId(rid,xIPtr,yIPtr);
    addElmtResidualByDmdaInd(*xIPtr,*yIPtr,residualPtr,fPtr);
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
        m_geoParam[0]=t_meshDesPtr->s_size_json.at("span");
        m_geoParam[1]=t_meshDesPtr->s_size_json.at("amplitude");
        m_geoParam[2]=t_meshDesPtr->s_size_json.at("width");
        break;
    case MeshShape::HALFSIN:
        m_geoParam[0]=t_meshDesPtr->s_size_json.at("span");
        m_geoParam[1]=t_meshDesPtr->s_size_json.at("amplitude");
        m_geoParam[2]=t_meshDesPtr->s_size_json.at("width");
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
    PetscCall(DMCreateGlobalVector(m_dm,&m_nodes_uInc1));
    PetscCall(DMCreateGlobalVector(m_dm,&m_Residual));
    /**************************************************************************/
    /** Create AMatrix managed by DMDA)****************************************/
    /**************************************************************************/

    /**************************************************************************/
    /** cal node and element num***********************************************/
    /**************************************************************************/
    m_mNodes=nx*ny;
    m_mNodes_p=nx*localNy[m_rank];
    m_mElmts=(nx-1)*(ny-1);
    PetscInt nx_elmt,ny_elmt;   /**< element num in x,y direction of this m_rank*/
    nx_elmt=nx-1;
    if(m_rank==m_rankNum-1)ny_elmt=ny-1;
    else ny_elmt=ny;
    m_mElmts_p=nx_elmt*ny_elmt;
    /**************************************************************************/
    /** set node and elemnt's global id, and element's connectivity************/
    /**************************************************************************/
    PetscCall(DMDAGetLocalInfo(m_dm,&m_daInfo));
    m_node_gId.resize(m_mNodes_p);
    m_elmt_gId.resize(m_mElmts_p);
    m_elmt_cnn.resize(m_mElmts_p);

    const int mNodePElmt=4; /**< node num per element*/
    PetscInt nodeRId=0,elmtRId=0;                   /**< current ndoe/element's m_rank id*/
    PetscInt nodeGId=m_daInfo.mx*m_daInfo.ys;       /**< current node's global id*/
    PetscInt elmtGId=(m_daInfo.mx-1)*m_daInfo.ys;   /**< current node's global id*/
    PetscCall(PetscSynchronizedPrintf(              /**< print the node and element's global id range of every m_rank*/
            PETSC_COMM_WORLD,
            "[%d]: owned node global id: %d -> %d \t owned elment global id: %d -> %d\n",
            m_rank, nodeGId, nodeGId+m_mNodes_p, elmtGId, m_mElmts_p
    ));
    PetscCall(PetscSynchronizedPrintf(              /**< print the node and element's DMDA index range of every m_rank*/
            PETSC_COMM_WORLD,
            "[%d]: owned node DMDA index: [%d, %d] X [%d, %d] \t owned elment DMDA index: [%d, %d] X [%d, %d]\n",
            m_rank, m_daInfo.xs, m_daInfo.ys, m_daInfo.xs+m_daInfo.xm, m_daInfo.ys+m_daInfo.ym,                     // node DMDA index range
                    m_daInfo.xs, m_daInfo.ys, m_daInfo.xs+m_daInfo.xm-1, m_daInfo.ys+m_daInfo.ym-(m_rank==m_rankNum-1)
    ));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));
    PetscScalar ***aCoord0;         /** ptr to the node's coords in ref config*/
    Vector2d CoordLocal;          /** store coord0 of a node*/
    PetscCall(DMDAVecGetArrayDOF(m_dm,m_nodes_coord0,&aCoord0));
    for(int yI=m_daInfo.ys;yI<m_daInfo.ys+m_daInfo.ym;yI++){    // loop over every row
        for(int xI=m_daInfo.xs;xI<m_daInfo.xs+m_daInfo.xm;xI++){// loop over every col
            m_node_gId[nodeRId]=nodeGId;        // set node's global id
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
            default:
                break;
            }
            for(int dofI=0;dofI<m_mDof_node;dofI++){
                aCoord0[yI][xI][dofI]=CoordLocal(dofI);
            }
            ++nodeRId;
            ++nodeGId;
        }
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
    *nodeNum=4;
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(vType,state);
    if(!arrayPtrRef){
        MessagePrinter::printErrorTxt("variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    for(int dofI=0;dofI<m_mDof_node;dofI++){
        variablePtr[0](dofI)=arrayPtrRef[yI][xI][dofI];
        variablePtr[1](dofI)=arrayPtrRef[yI][xI+1][dofI];
        variablePtr[2](dofI)=arrayPtrRef[yI+1][xI+1][dofI];
        variablePtr[3](dofI)=arrayPtrRef[yI+1][xI][dofI];
    }
    return 0;
}

PetscErrorCode StructuredMesh2D::getNodeVariableByDmdaInd(NodeVariableType vType, PetscInt xI,PetscInt yI,int state,Vector *variablePtr){
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(vType,state);
    if(!arrayPtrRef){
        MessagePrinter::printErrorTxt("variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    for(int dofI=0;dofI<m_mDof_node;dofI++){
        (*variablePtr)(dofI)=arrayPtrRef[yI][xI][dofI];
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
    if(yI+m_dim-1>m_daInfo.gym-1||xI+m_dim-1>m_daInfo.gxm-1||
        yI<m_daInfo.gys||xI<m_daInfo.gxs){
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"DMDA (row,col)=(%d,%d) is out of range, max index = (%d,%d)",
                            xI+m_dim-1,yI+m_dim-1,m_daInfo.gxm-1,m_daInfo.gym-1));
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
            for(colNodeI=0;colNodeI<m_mDof_node;++colNodeI){
                for(nodeDofI2=0;nodeDofI2<m_mDof_node;++nodeDofI2){
                    entry[mDofPerElmt*rowDofI+colDofI]=(*matrixPtr)(rowDofI,colDofI);
                    ++colDofI;
                }
            }
            ++rowDofI;      
        }
    }
    // for(rowJ=yI;rowJ<yI+m_dim;++rowJ){
    //     for(rowI=xI;rowI<xI+m_dim;++rowI){
    //         for(nodeDofI=0;nodeDofI<m_mDof_node;++nodeDofI){
    //             row[rowDofI].i=rowI;        col[rowDofI].i=rowI;       
    //             row[rowDofI].j=rowJ;        col[rowDofI].j=rowJ;    
    //             row[rowDofI].c=nodeDofI;    col[rowDofI].c=nodeDofI;
    //             colDofI=0;
    //             for(colJ=yI;colJ<yI+m_dim;++colJ){
    //                 for(colI=xI;colI<xI+m_dim;++colI){
    //                     for(nodeDofI2=0;nodeDofI2<m_mDof_node;++nodeDofI2){
    //                         entry[mDofPerElmt*rowDofI+colDofI]=(*matrixPtr)(rowDofI,colDofI);
    //                         ++colDofI;
    //                     }
    //                 }
    //             }
    //             ++rowDofI;
    //         }
    //     }
    // }
    PetscCall(MatSetValuesStencil(*APtr,mDofPerElmt,row,mDofPerElmt,col,entry,ADD_VALUES));
    return 0;
}

PetscErrorCode StructuredMesh2D::addElmtResidualByDmdaInd(PetscInt xI,PetscInt yI,Vector *residualPtr,Vec *fPtr){
    PetscScalar ***& arrayPtrRef=getNodeVariablePtrRef(NodeVariableType::RESIDUAL,1);
    if(fPtr){}
    if(!arrayPtrRef){
        MessagePrinter::printErrorTxt("variable array need to point to Vec before use it.");
        MessagePrinter::exitcfem();
    }
    if(yI+1>m_daInfo.gym-1||xI+1>m_daInfo.gxm-1){
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,"DMDA (row,col)=(%d,%d) is out of range, max index = (%d,%d)",
                            xI+1,yI+1,m_daInfo.gxm-1,m_daInfo.gym-1));
        MessagePrinter::exitcfem();            
    }
    for(int dofI=0;dofI<m_mDof_node;dofI++){
        arrayPtrRef[yI][xI][dofI]+=residualPtr[0](dofI);
        arrayPtrRef[yI][xI+1][dofI]+=residualPtr[1](dofI);
        arrayPtrRef[yI+1][xI+1][dofI]+=residualPtr[2](dofI);
        arrayPtrRef[yI+1][xI][dofI]+=residualPtr[3](dofI);
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
    (*aCoord)(0)=xI*m_geoParam[0]/m_daInfo.mx;
    (*aCoord)(1)=yI*m_geoParam[1]/m_daInfo.my;
}

void StructuredMesh2D::getSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span l, amplitude a, width w in order)*/
    (*aCoord)(1)=yI*m_geoParam[0]/m_daInfo.my;
    /** x=a*sin(pi*y/l)*/
    (*aCoord)(0)=m_geoParam[1]*sin(M_PI*(*aCoord)(1)/m_geoParam[0])+xI*m_geoParam[2]/m_daInfo.mx;
}

void StructuredMesh2D::getHalfSinCoord0ByDmdaInd(PetscInt xI,PetscInt yI,Vector *aCoord){
    /**m_geoParam < geometry parameter to describe the regular mesh domain
     * (for rectangular domain, is x length,y length in order; 
     * for sin or half sin domain, is span l, amplitude a, width w in order)*/
    (*aCoord)(1)=yI*m_geoParam[0]/m_daInfo.my;
    /** x=a*sin(pi*y/2l)*/
    (*aCoord)(0)=m_geoParam[1]*sin(M_PI*(*aCoord)(1)/(2*m_geoParam[0]))+xI*m_geoParam[2]/m_daInfo.mx;  
}

void StructuredMesh2D::getElmtDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr){
    *yIPtr=rId/m_daInfo.mx+m_daInfo.ys;
    *xIPtr=rId%m_daInfo.mx+m_daInfo.xs;
}

void StructuredMesh2D::getNodeDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr){
    *yIPtr=rId/(m_daInfo.mx-1)+m_daInfo.ys;
    *xIPtr=rId%(m_daInfo.mx-1)+m_daInfo.xs;
}

Vec & StructuredMesh2D::getNodeLocalVariableVecRef(NodeVariableType vType, int state){
    switch (vType)
    {
    case NodeVariableType::COORD:
        switch (state)
        {
        case 0:
            return m_nodes_coord0;
            break;
        case 2:
            return m_nodes_coord2;
            break;       
        default:
            MessagePrinter::printErrorTxt("required Coords vec state must be 0 or 2.");
            MessagePrinter::exitcfem();
            break;
        }
        break;
    case NodeVariableType::UINC:
        switch (state)
        {
        case 1:
            return m_nodes_uInc1;
            break;    
        default:
            MessagePrinter::printErrorTxt("required incremental U vec state must be 1.");
            MessagePrinter::exitcfem();
            break;
        }
        break;      
    case NodeVariableType::RESIDUAL:
        if (state==1){
           return m_Residual;
        }
        else{
            MessagePrinter::printErrorTxt("required residual vec state must be 1.");
            MessagePrinter::exitcfem();           
        }
        break;
    default:
        MessagePrinter::printErrorTxt("required node variable vec are not unsupported in current mesh system");
        MessagePrinter::exitcfem();       
        break;
    }
    return m_nodes_coord0;    
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
    case NodeVariableType::UINC:
        switch (state)
        {
        case 1:
            return m_array_nodes_uInc1;
            break;     
        default:
            MessagePrinter::printErrorTxt("required incremental U array address state must be 1.");
            MessagePrinter::exitcfem();
            break;
        }
        break;      
    case NodeVariableType::RESIDUAL:
        if (state==1){
            return m_array_nodes_residual;
        }
        else{
            MessagePrinter::printErrorTxt("required residual array address state must be 1.");
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