#include "PostProcessSystem/PostProcessSystem.h"
PostProcessSystem::PostProcessSystem():
    m_ifMeshSysSet(false),m_ifElmtSysSet(false),m_ifProjValVec(false),
    m_meshSysPtr(nullptr),m_elmtSysPtr(nullptr),m_elmtVarType(ElementVariableType::NONE),
    m_array_proj_weight(nullptr),m_array_proj_val(nullptr){
        MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
        MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);
    };

PostProcessSystem::PostProcessSystem(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr):
    m_ifMeshSysSet(false),m_ifElmtSysSet(false),m_ifProjValVec(false),
    m_meshSysPtr(nullptr),m_elmtSysPtr(nullptr),m_elmtVarType(ElementVariableType::NONE),
    m_array_proj_weight(nullptr),m_array_proj_val(nullptr){
    m_meshSysPtr=t_meshSysPtr;
    m_ifMeshSysSet=true;
    m_elmtSysPtr=t_elmtSysPtr;
    m_ifElmtSysSet=true;
    MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&m_rankNum);
}
PostProcessSystem::~PostProcessSystem(){
    VecDestroy(&m_proj_weight);
    VecDestroy(&m_proj_val);
}