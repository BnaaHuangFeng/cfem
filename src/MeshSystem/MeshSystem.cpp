#include"MeshSystem/MeshSystem.h"
#include"petsc.h"
#include"Utils/MessagePrinter.h"
#include<vector>
MeshSystem::MeshSystem(){
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);
    m_ifSaveMesh=true;
}
MeshSystem::MeshSystem(Timer *timerPtr){
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);
    m_ifSaveMesh=true;
    m_timerPtr=timerPtr;   
}
MeshSystem::~MeshSystem(){
    VecDestroy(&m_nodes_coord0);
    VecDestroy(&m_nodes_coord1);
    VecDestroy(&m_nodes_coord2);
    DMDestroy(&m_dm);
}
bool MeshSystem::getElmtCnnByRId(int rId,const vector<int> *elmtCnn){
    checkElmtRId(rId);
    elmtCnn=&(m_elmt_cnn[rId]);
    // get rid of unused warning
    if(elmtCnn){}
    return true;
}
void MeshSystem::checkNodeRId(int rId){
    if(rId<0||rId>m_mNodes_p){
        string errStr="[rank "+to_string(m_rank)+"]: "+
        "node id " +to_string(rId)+" in rank is out of range, max is "+to_string(m_mNodes_p);
        MessagePrinter::printErrorTxt(errStr);
        MessagePrinter::exitcfem();
    }
}
void MeshSystem::checkElmtRId(int rId){
    if(rId<0||rId>m_mElmts_p){
        string errStr="[rank "+to_string(m_rank)+"]: "+
        "element id " +to_string(rId)+" in rank is out of range, max is "+to_string(m_mElmts_p);
        MessagePrinter::printErrorTxt(errStr);
        MessagePrinter::exitcfem();
    }
}