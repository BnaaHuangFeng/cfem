#include "PostProcessSystem/PostStructured2d.h"
#include "MaterialSystem/ElmtVarInfo.h"
#include "PostProcessSystem/OutputVarInfo.h"
#include "MeshSystem/NodeVarInfo.h"
#include <fstream>
#include "MeshSystem/StructuredMesh2D.h"
using namespace std;
PetscErrorCode PostStructured2d::outputFieldVariable(int t_increI, PetscScalar t_t){
    int interval=m_outputDesPtr->s_FD.s_interval;
    if(t_increI%interval!=0)return 0;
    if(t_t){}
    string fileName;
    FieldOutputFormat fieldFormat=m_outputDesPtr->s_FD.s_format;
    fileName=fieldOutputFileName(t_increI,fieldFormat);
    std::ofstream out;
    StructuredMesh2D *mesh2dPtr=(StructuredMesh2D *)m_meshSysPtr;
    if(m_rank==0){
        openOutputFile(fileName,&out,std::ios::out);
    //****************************************
    //*** print out header information
    //****************************************
        out<<"<?xml version=\"1.0\"?>\n";
        out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        out<<"<UnstructuredGrid>\n";
        out<<"<Piece NumberOfPoints=\""<<m_meshSysPtr->m_mNodes
            <<"\" NumberOfCells=\""<<m_meshSysPtr->m_mElmts<<"\">\n";

        out<<"<Points>\n";
        out<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";
        out.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    //*****************************
    // print out node coordinates
    //*****************************
    PetscScalar ***aCoord;
    PetscCall(DMDAVecGetArrayDOFRead(m_meshSysPtr->m_dm,m_meshSysPtr->m_nodes_coord2,&aCoord));
    for(int rankI=0;rankI<m_rankNum;rankI++){
        if(m_rank==rankI){// loop over all rank. print nodes coords if it's this rank's turn
            openOutputFile(fileName,&out,std::ios::app);
            out<<std::scientific<<std::setprecision(6);
            PetscInt xs=m_dmInfo.xs,    xe=m_dmInfo.xs+m_dmInfo.xm,
                     ys=m_dmInfo.ys,    ye=m_dmInfo.ys+m_dmInfo.ym;
            for(PetscInt yI=ys;yI<ye;yI++){//loop over row in this rank
                for(PetscInt xI=xs;xI<xe;xI++){//loop over col in this rank
                    out<<aCoord[yI][xI][0]<<" ";
                    out<<aCoord[yI][xI][1]<<" ";
                    out<<0.0<<"\n";
                }
            }
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // ensure output coords in order
    }
    if(m_rank==0){
        openOutputFile(fileName,&out,std::ios::app);
        out<<"</DataArray>\n";
        out<<"</Points>\n";
    //***************************************
    //*** For cell information
    //***************************************
        out<<"<Cells>\n";
        out<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        out.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    for(int rankI=0;rankI<m_rankNum;rankI++){
        if(m_rank==rankI){// loop over all rank. print element connectivity if it's this rank's turn
            openOutputFile(fileName,&out,std::ios::app);
            out<<std::scientific<<std::setprecision(6);
            for(PetscInt eI=0;eI<m_meshSysPtr->m_mElmts_p;eI++){ // loop over elmts in this rank
                PetscInt mNodeInElmt=m_meshSysPtr->m_elmt_cnn[eI].size();
                for(PetscInt nI=0;nI<mNodeInElmt;nI++){// loop over nodes in this elmt
                    out<<m_meshSysPtr->m_elmt_cnn[eI][nI]<<" ";
                }
                out<<"\n";
            }
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // ensure output connectivity in order
    }
    if(m_rank==0){
        openOutputFile(fileName,&out,std::ios::app);
        out<<"</DataArray>\n";
    //***************************************
    //*** for offset
    //***************************************
        out<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        out.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    const PetscInt mNodeInElmt=4;
    PetscInt offsetBase=((m_dmInfo.mx-1)*m_dmInfo.ys+m_dmInfo.xs)*mNodeInElmt;
    for(int rankI=0;rankI<m_rankNum;rankI++){
        if(m_rank==rankI){// loop over all rank. print element connectivity if it's this rank's turn
            openOutputFile(fileName,&out,std::ios::app);
            for(PetscInt eI=0;eI<m_meshSysPtr->m_mElmts_p;eI++){ // loop over elmts in this rank
                offsetBase+=mNodeInElmt;
                out<<offsetBase<<"\n";
            }
            
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // ensure output offset in order
    }
    if(m_rank==0){
        openOutputFile(fileName,&out,std::ios::app);
        out<<"</DataArray>\n";
    //***************************************
    //*** for VTKCellType
    //***************************************
        static const int vtkType=9;             /**< vtk cell type*/
        out<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(PetscInt eI=0;eI<m_meshSysPtr->m_mElmts;eI++){ // loop over elmts in all rank
            out<<vtkType<<"\n";
        }
        out<<"</DataArray>\n";
        out<<"</Cells>\n";
        out.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    //***************************************
    //*** write field variable name
    //***************************************
    string scalarNameSeq="",vectorNameSeq="",tensorNameSeq="";
    for(string name : m_infoForFieldOut.scalarName) scalarNameSeq+=name+" ";
    for(string name : m_infoForFieldOut.vectorName) vectorNameSeq+=name+" ";
    for(string name : m_infoForFieldOut.tensorName) tensorNameSeq+=name+" ";
    if(m_rank==0){
        openOutputFile(fileName,&out,std::ios::app);
        out<<"<PointData ";
        if(scalarNameSeq!=""){
            out<<" Scalar=\""+scalarNameSeq+"\"";
        }
        if(vectorNameSeq!=""){
            out<<" Vector=\""+vectorNameSeq+"\"";
        }
        if(tensorNameSeq!=""){
            out<<" Tensor=\""+tensorNameSeq+"\"";
        }
        out<<">\n";
        out.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    //***************************************
    //*** write data of field variable in node
    //*************************************** 
    int mNodeVarType=m_infoForFieldOut.varInNodeVec.size();
    for(int i=0;i<mNodeVarType;++i){
        Vec *globalVecPtr=m_meshSysPtr->globalVecPtr(m_infoForFieldOut.varInNodeVec[i],2);
        mesh2dPtr->openNodeVariableVec(m_infoForFieldOut.varInNodeVec[i],globalVecPtr,2,VecAccessMode::READ);
        PetscScalar *** &varArray=mesh2dPtr->getNodeVariablePtrRef(m_infoForFieldOut.varInNodeVec[i],2);
        if(m_rank==0){
            openOutputFile(fileName,&out,std::ios::app);
            out<<"<DataArray type=\"Float64\" Name=\"" << m_infoForFieldOut.varNameInNodeVec[i] << "\"  ";
            out<<"NumberOfComponents=\""+to_string(m_infoForFieldOut.varCpntInNodeVec[i])+"\" format=\"ascii\">\n";
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.        
        for(int rankI=0;rankI<m_rankNum;rankI++){
            if(m_rank==rankI){// loop over all rank. print element connectivity if it's this rank's turn
                openOutputFile(fileName,&out,std::ios::app);
                out<<std::scientific<<std::setprecision(6);
                PetscInt xs=m_dmInfo.xs,    xe=m_dmInfo.xs+m_dmInfo.xm,
                         ys=m_dmInfo.ys,    ye=m_dmInfo.ys+m_dmInfo.ym;  
                for(PetscInt yI=ys;yI<ye;yI++){//loop over row in this rank
                    for(PetscInt xI=xs;xI<xe;xI++){//loop over col in this rank
                        out<<varArray[yI][xI][0]<<" ";
                        out<<varArray[yI][xI][1]<<" ";
                        out<<0.0<<"\n";
                    }
                }                              
                out.close();
            }
            PetscCall(PetscBarrier(NULL));
        }
        if(m_rank==0){
            openOutputFile(fileName,&out,std::ios::app);
            out << "</DataArray>\n\n";
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.          
        mesh2dPtr->closeNodeVariableVec(m_infoForFieldOut.varInNodeVec[i],globalVecPtr,2,VecAccessMode::READ);
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    //***************************************
    //*** write data of field variable in elmt
    //***************************************     
    int mElmtVarType=m_infoForFieldOut.varInElmtVec.size();
    for(int i=0;i<mElmtVarType;++i){
        projElmtVariable(m_infoForFieldOut.varInElmtVec[i]);
        openNodeVariableVec(&m_proj_vec,&m_proj_val_local,&m_array_proj_val,m_infoForFieldOut.varCpntInElmtVec[i],VecAccessMode::READ);
        if(m_rank==0){
            openOutputFile(fileName,&out,std::ios::app);
            out<<"<DataArray type=\"Float64\" Name=\"" << m_infoForFieldOut.varNameInElmtVec[i] << "\"  ";
            out<<"NumberOfComponents=\""+to_string(m_infoForFieldOut.varCpntInElmtVec[i])+"\" format=\"ascii\">\n";
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.        
        for(int rankI=0;rankI<m_rankNum;rankI++){
            if(m_rank==rankI){// loop over all rank. print element connectivity if it's this rank's turn
                openOutputFile(fileName,&out,std::ios::app);
                out<<std::scientific<<std::setprecision(6);
                PetscInt xs=m_dmInfo.xs,    xe=m_dmInfo.xs+m_dmInfo.xm,
                         ys=m_dmInfo.ys,    ye=m_dmInfo.ys+m_dmInfo.ym;  
                for(PetscInt yI=ys;yI<ye;yI++){//loop over row in this rank
                    for(PetscInt xI=xs;xI<xe;xI++){//loop over col in this rank
                        for(PetscInt cpntI=0;cpntI<m_infoForFieldOut.varCpntInElmtVec[i];++cpntI){
                            out<<m_array_proj_val[yI][xI][cpntI]<<" ";
                        }
                        out<<endl;
                    }
                }                              
                out.close();
            }
            PetscCall(PetscBarrier(NULL));
        }
        if(m_rank==0){
            openOutputFile(fileName,&out,std::ios::app);
            out << "</DataArray>\n\n";
            out.close();
        }
        PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.          
        closeNodeVariableVec(&m_proj_vec,&m_proj_val_local,&m_array_proj_val,m_infoForFieldOut.varCpntInElmtVec[i],VecAccessMode::READ);
        projVecClean();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    //***************************************
    //*** end writting
    //***************************************
    if(m_rank==0){
        openOutputFile(fileName,&out,std::ios::app);
        out<< "</PointData>\n";
        out<<"</Piece>\n";
        out<<"</UnstructuredGrid>\n";
        out<<"</VTKFile>"<<endl;
        out.close();
    }
    PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
    PetscCall(DMDAVecRestoreArrayDOFRead(m_meshSysPtr->m_dm,m_meshSysPtr->m_nodes_coord2,&aCoord));
    return 0;    
}
PetscErrorCode PostStructured2d::outputHisVariable(int t_increI, PetscScalar t_t){
    if(!t_increI)return 0;
    int mHisNodeVar=m_infoForHisOut.varInNodeVec.size();
    int mHisElmtVar=m_infoForHisOut.varInElmtVec.size();
    for(int nodeVarI=0;nodeVarI<mHisNodeVar;++nodeVarI){
        int interval=m_infoForHisOut.intervalInNodeVec[nodeVarI];
        if(t_increI%interval!=0&&t_t!=m_loadCtrlPtr->factorFinal())continue;
        HistoryVariableType hisVarType=m_infoForHisOut.hisVarInNodeVec[nodeVarI];
        NodeVariableType nodeVarType=m_infoForHisOut.varInNodeVec[nodeVarI];
        int processedDataize=m_infoForHisOut.dataNumPerFrameInNodeVec[nodeVarI];   // data num in this frame
        int mCpntPerData=m_infoForHisOut.mCpntPerDataInNodeVec[nodeVarI];
        int cpntInd=m_infoForHisOut.varCpntIndInNodeVec[nodeVarI];
        VarOutputForm outputForm=m_infoForHisOut.outputFormInNodeVec[nodeVarI];
        vector<int> set=m_meshSysPtr->m_setManager.getSet(
            m_infoForHisOut.setNameInNodeVec[nodeVarI],SetType::NODE);
        PetscInt setSize=set.size();
        PetscScalar *timeBuf=m_infoForHisOut.timeBuffInNodeVec[nodeVarI];
        PetscScalar *buf=m_infoForHisOut.bufferInNodeVec[nodeVarI];
        int buffFrameNum=m_infoForHisOut.bufferFrameNumInNodeVec[nodeVarI]; /**< current data num in buffer*/
        timeBuf[buffFrameNum]=t_t;
        genNodeVariable(nodeVarType);
        switch (outputForm){
            case VarOutputForm::ANY:
            case VarOutputForm::SUM:
                if(outputForm==VarOutputForm::ANY)setSize=min(1,setSize);
                switch(nodeVarType){
                    case NodeVariableType::U:
                    case NodeVariableType::RF:{
                        switch (cpntInd){
                            case -1:
                                buf[(processedDataize*mCpntPerData)*buffFrameNum+0]=0.0;
                                buf[(processedDataize*mCpntPerData)*buffFrameNum+1]=0.0;
                                buf[(processedDataize*mCpntPerData)*buffFrameNum+2]=0.0;
                                for(PetscInt varI=0;varI<setSize;++varI){
                                    PetscInt xI=0, yI=0;
                                    getNodeDmdaIndByRId(set[varI],&xI,&yI);
                                    buf[(processedDataize*mCpntPerData)*buffFrameNum+0]+=m_array_his_node[yI][xI][0];
                                    buf[(processedDataize*mCpntPerData)*buffFrameNum+1]+=m_array_his_node[yI][xI][1];                     
                                }
                                break;
                            default:
                                if(cpntInd<2){
                                    buf[(processedDataize*mCpntPerData)*buffFrameNum+0]=0.0;
                                    for(PetscInt varI=0;varI<setSize;++varI){
                                        PetscInt xI=0, yI=0;
                                        getNodeDmdaIndByRId(set[varI],&xI,&yI);
                                        buf[(processedDataize*mCpntPerData)*buffFrameNum+0]+=m_array_his_node[yI][xI][cpntInd];                    
                                    }                            
                                }
                                break;
                        }
                    }
                        break;
                    default:
                        MessagePrinter::printErrorTxt("PostStructured2d: Historic variable of this kind is not developed");
                        MessagePrinter::exitcfem();
                        break;
                }
                break;
            default:
                MessagePrinter::printErrorTxt("PostStructured2d: VarOutputForm of this kind is not developed");
                MessagePrinter::exitcfem();            
                break;
        }
        restoreNodeVariable(nodeVarType);
        ++buffFrameNum;
        ++m_infoForHisOut.bufferFrameNumInNodeVec[nodeVarI];
        if(buffFrameNum<m_hisBuffLen&&t_t!=m_loadCtrlPtr->factorFinal())continue;
        // process the buffer
        PetscScalar *processedData=nullptr;
        switch(outputForm){
            case VarOutputForm::ANY:
                processedData=buf;
                break;
            case VarOutputForm::SUM:{
                Mat framesMat;
                PetscInt mcol = mCpntPerData*buffFrameNum;
                PetscCall(MatCreateDense(PETSC_COMM_WORLD,1,mcol,m_rankNum,mcol,NULL,&framesMat));
                PetscCall(MatSetUp(framesMat));
                PetscInt *colLocalInd=new PetscInt[mcol];
                for(int i=0;i<mcol;++i)colLocalInd[i]=i;
                PetscInt rowLocalInd=0;
                PetscCall(MatSetValuesLocal(framesMat,1,&rowLocalInd,mcol,colLocalInd,buf,INSERT_VALUES));
                PetscCall(MatAssemblyBegin(framesMat,MAT_FINAL_ASSEMBLY));
                PetscCall(MatAssemblyEnd(framesMat,MAT_FINAL_ASSEMBLY));
                // PetscCall(MatView(framesMat,PETSC_VIEWER_STDOUT_WORLD));
                PetscCall(MatGetColumnSums(framesMat,buf));
                processedData=buf;
                delete[] colLocalInd;
                PetscCall(MatDestroy(&framesMat));
                }
                break;
            default:
                MessagePrinter::printErrorTxt("PostStructured2d: VarOutputForm of this kind is not developed");
                MessagePrinter::exitcfem();            
                break;            
        }
        // output the buff
        switch(outputForm){
            case VarOutputForm::ANY:{
                if(setSize>0){
                    std::ofstream out;
                    string fileName=hisOutputFileName(m_infoForHisOut.setNameInNodeVec[nodeVarI],hisVarType,m_infoForHisOut.outFormatInNodeVec[nodeVarI]);
                    fileName="rank-"+to_string(m_rank)+"-"+fileName;
                    openOutputFile(fileName,&out,std::ios::app);
                    out<<std::scientific<<std::setprecision(6);
                    for(int frameI=0;frameI<buffFrameNum;frameI++){
                            out<<timeBuf[frameI];
                        for(int cpntI=0;cpntI<mCpntPerData;++cpntI){
                            out<<", ";
                            out<<processedData[mCpntPerData*frameI+cpntI];
                        }
                        out<<endl;
                    }
                    out.close();
                }
                break;
            }    
            case VarOutputForm::SUM:{
                if(m_rank==0){
                    std::ofstream out;
                    string fileName=hisOutputFileName(m_infoForHisOut.setNameInNodeVec[nodeVarI],hisVarType,m_infoForHisOut.outFormatInNodeVec[nodeVarI]);
                    openOutputFile(fileName,&out,std::ios::app);
                    out<<std::scientific<<std::setprecision(6);
                    for(int frameI=0;frameI<buffFrameNum;frameI++){
                            out<<timeBuf[frameI];
                        for(int cpntI=0;cpntI<mCpntPerData;++cpntI){
                            out<<", ";
                            out<<processedData[mCpntPerData*frameI+cpntI];
                        }
                        out<<endl;
                    }
                    out.close();
                }
                PetscCall(PetscBarrier(NULL));  // wait for rank 0 writting.
                break;
            }
            default:
                break;
        }
        // clear buffer;
        m_infoForHisOut.bufferFrameNumInNodeVec[nodeVarI]=0;
    }
    if(mHisElmtVar){}
    return 0;
}
void PostStructured2d::openOutputFile(string fileName, ofstream *ofPtr,ios_base::openmode mode){
    char buff[110];
    string str;
    ofPtr->open(fileName,mode);
    if(!ofPtr->is_open()){
        snprintf(buff,110,"can\'t write mesh to vtu file(=%28s), please make sure you have the write permission",fileName.c_str());
        str=buff;
        MessagePrinter::printErrorTxt(str);
        MessagePrinter::exitcfem();
    }    
}