#include "PostProcessSystem/PostProcessSystem.h"
#include "MaterialSystem/ElmtVarInfo.h"
#include "PostProcessSystem/OutputVarInfo.h"
#include "MeshSystem/NodeVarInfo.h"
#include "cstdio"
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
PostProcessSystem::PostProcessSystem(OutputDescription *t_outputDesPtr):
    m_ifMeshSysSet(false),m_ifElmtSysSet(false),
    m_ifProjVec(false),m_ifHisNodeVec(false),m_ifHisElmtVec(false),
    m_meshSysPtr(nullptr),m_elmtSysPtr(nullptr),m_loadCtrlPtr(nullptr),
    m_array_proj_weight(nullptr),m_array_proj_val(nullptr),
    m_array_his_node(nullptr),m_array_his_elmt(nullptr){
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);
    m_prefix=t_outputDesPtr->s_outPrefix;
    m_ifOutputDesSet=false;
    readOutputDes(t_outputDesPtr);
};

PostProcessSystem::PostProcessSystem(OutputDescription *t_outputDesPtr,MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr, LoadController *t_loadCtrlPtr):
    m_ifMeshSysSet(false),m_ifElmtSysSet(false),
    m_ifProjVec(false),m_ifHisNodeVec(false),m_ifHisElmtVec(false),
    m_meshSysPtr(nullptr),m_elmtSysPtr(nullptr),m_loadCtrlPtr(nullptr),
    m_array_proj_weight(nullptr),m_array_proj_val(nullptr),
    m_array_his_node(nullptr),m_array_his_elmt(nullptr){
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);
    m_meshSysPtr=t_meshSysPtr;
    m_ifMeshSysSet=true;
    m_elmtSysPtr=t_elmtSysPtr;
    m_ifElmtSysSet=true;
    m_loadCtrlPtr=t_loadCtrlPtr;
    m_ifLoadCtrlSet=true;
    m_prefix=t_outputDesPtr->s_outPrefix;
    m_ifOutputDesSet=false;
    readOutputDes(t_outputDesPtr);
}
PostProcessSystem::~PostProcessSystem(){
    VecDestroy(&m_proj_weight);
    VecDestroy(&m_his_node_vec);
    if(m_ifProjVec)
        VecDestroy(&m_proj_vec);
    int mhisNodeVar=m_infoForHisOut.hisVarInNodeVec.size();
    for(int i=0;i<mhisNodeVar;++i){
        delete[]m_infoForHisOut.bufferInNodeVec[i];
        delete[]m_infoForHisOut.timeBuffInNodeVec[i];
    }
    int mhisElmtVar=m_infoForHisOut.hisVarInElmtVec.size();
    for(int i=0;i<mhisElmtVar;++i){
        delete[]m_infoForHisOut.bufferInElmtVec[i];
        delete[]m_infoForHisOut.timeBuffInElmtVec[i];        
    }

}
PetscErrorCode PostProcessSystem::readOutputDes(OutputDescription *t_outputDesPtr){
    if(m_ifOutputDesSet)return 0;
    m_outputDesPtr=t_outputDesPtr;
    // read field output variable info
    FieldOutputDescription *fieldDesPtr=&t_outputDesPtr->s_FD;
    for(FieldVariableType fieldVarType : fieldDesPtr->s_varTypes){
        VarMathType varMathType=FieldVarInfo::varMathType.find(fieldVarType)->second;
        VarPosition varPosition=FieldVarInfo::varPosition.find(fieldVarType)->second;
        int varTypeInd         =FieldVarInfo::varType.find(fieldVarType)->second;
        string varName;
        int mCpnt=0;
        switch (varPosition)
        {
        case VarPosition::NODE:
            varName=NodeVarInfo::nodeVarName.find((NodeVariableType)varTypeInd)->second;
            mCpnt=NodeVarInfo::nodeVarCpntNum.find((NodeVariableType)varTypeInd)->second;
            m_infoForFieldOut.varInNodeVec.push_back(NodeVariableType(varTypeInd));
            m_infoForFieldOut.varNameInNodeVec.push_back(varName);
            m_infoForFieldOut.varCpntInNodeVec.push_back(mCpnt);
            break;
        case VarPosition::ELEMENT:
            varName=ElmtVarInfo::elmtVarName.find((ElementVariableType)varTypeInd)->second;
            mCpnt=ElmtVarInfo::elmtVarCpntNum.find((ElementVariableType)varTypeInd)->second;
            m_infoForFieldOut.varInElmtVec.push_back(ElementVariableType(varTypeInd));
            m_infoForFieldOut.varNameInElmtVec.push_back(varName);
            m_infoForFieldOut.varCpntInElmtVec.push_back(mCpnt);
            break;
        default:
            break;
        }
        switch (varMathType)
        {
        case VarMathType::SCALAR:
            m_infoForFieldOut.scalarName.push_back(varName);
            break;
        case VarMathType::TENSORRANK2:
        case VarMathType::TENSORRANK4:
            m_infoForFieldOut.tensorName.push_back(varName);
            break;
        case VarMathType::VECTOR:
            m_infoForFieldOut.vectorName.push_back(varName);
            break;
        default:
            break;
        }
    }
    // read historic output variable info
    for(HistoryOutputDescription hisDes : t_outputDesPtr->s_HD){
        int intval=hisDes.s_interval;
        string setName=hisDes.s_setName;
        HistoryOutputFormat outFormat=hisDes.s_format;
        for(HistoryVariableType hisVarType : hisDes.s_varTypes){
            VarPosition varPosition=HisVarInfo::varPosition.find(hisVarType)->second;
            int varTypeInd         =HisVarInfo::varType.find(hisVarType)->second;
            string varName;
            int mCpnt=0;
            int cpntIndex=HisVarInfo::cpntInd.find(hisVarType)->second;
            VarOutputForm outputForm=HisVarInfo::varOutputForm.find(hisVarType)->second;
            string outFileName=m_prefix+"/"+hisOutputFileName(setName,hisVarType,outFormat);
            // delete old file if exists
            // if(m_rank==0)
                remove(outFileName.c_str());
            string outFileNameRank=m_prefix+"/rank-"+to_string(m_rank)+'-'+hisOutputFileName(setName,hisVarType,outFormat);
            remove(outFileNameRank.c_str());
            switch (varPosition)
            {
            case VarPosition::NODE:
                m_infoForHisOut.hisVarInNodeVec.push_back(hisVarType);
                varName=NodeVarInfo::nodeVarName.find((NodeVariableType)varTypeInd)->second;
                mCpnt=NodeVarInfo::nodeVarCpntNum.find((NodeVariableType)varTypeInd)->second;
                m_infoForHisOut.varInNodeVec.push_back(NodeVariableType(varTypeInd));
                m_infoForHisOut.varNameInNodeVec.push_back(varName);
                m_infoForHisOut.varCpntInNodeVec.push_back(mCpnt);
                m_infoForHisOut.varCpntIndInNodeVec.push_back(cpntIndex);
                m_infoForHisOut.intervalInNodeVec.push_back(intval);
                m_infoForHisOut.setNameInNodeVec.push_back(setName);
                m_infoForHisOut.outFormatInNodeVec.push_back(outFormat);
                m_infoForHisOut.outputFormInNodeVec.push_back(outputForm);
                m_infoForHisOut.fileNameInNodeVec.push_back(outFileName);
                break;
            case VarPosition::ELEMENT:
                m_infoForHisOut.hisVarInElmtVec.push_back(hisVarType);
                varName=ElmtVarInfo::elmtVarName.find((ElementVariableType)varTypeInd)->second;
                mCpnt=ElmtVarInfo::elmtVarCpntNum.find((ElementVariableType)varTypeInd)->second;
                m_infoForFieldOut.varInElmtVec.push_back(ElementVariableType(varTypeInd));
                m_infoForFieldOut.varNameInElmtVec.push_back(varName);
                m_infoForFieldOut.varCpntInElmtVec.push_back(mCpnt);
                m_infoForHisOut.varCpntIndInElmtVec.push_back(cpntIndex);
                m_infoForHisOut.intervalInElmtVec.push_back(intval);
                m_infoForHisOut.setNameInElmtVec.push_back(setName);
                m_infoForHisOut.outFormatInElmtVec.push_back(outFormat);  
                m_infoForHisOut.outputFormInElmtVec.push_back(outputForm);     
                m_infoForHisOut.fileNameInElmtVec.push_back(outFileName);      
                break;
            default:
                break;
            }            
        }
    }
    if(m_rank==0){
        if(access(m_prefix.c_str(),0)!=0){
            if(mkdir(m_prefix.c_str(),S_IRWXU|S_IRWXG|S_IRWXO)){
                MessagePrinter::printRankError("Directory "+m_prefix+" creating failed!");
                MessagePrinter::exitcfem();
            }
        }        
    }
    PetscCall(PetscBarrier(NULL)); // waiting for rank 0 to create folder
    m_ifOutputDesSet=true;
    return 0;
}
string PostProcessSystem::fieldOutputFileName(int t_increI,FieldOutputFormat format){
    string fileName=m_prefix;
    ostringstream ss;
    ss<<setfill('0')<<setw(8)<<t_increI;
    switch (format)
    {
    case FieldOutputFormat::VTU:
        ss<<".vtu";
        break;
    case FieldOutputFormat::VTK:
        ss<<".vtk";
        break;
    default:
        break;
    }
    fileName+="-"+ss.str();
    return fileName;
}
string PostProcessSystem::hisOutputFileName(string setName,HistoryVariableType vType, HistoryOutputFormat format){
    string fileName;
    string varName=HisVarInfo::varName.find(vType)->second;
    string append;
    switch (format)
    {
    case HistoryOutputFormat::CSV:
        append=".csv";
        break;
    case HistoryOutputFormat::MATLAB:
        append=".m";
        break;
    default:
        break;
    }
    fileName=setName+"-"+varName+append;
    return fileName;
}