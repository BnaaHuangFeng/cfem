#pragma once
#include "MeshSystem/MeshSystem.h"
#include "ElementSystem/ElementSystem.h"
#include "LoadController/LoadController.h"
#include "InputSystem/DescriptionInfo.h"
#include "PostProcessSystem/OutputVarInfo.h"
struct InfoForFieldOutput{
    vector<string>              scalarName;
    vector<string>              vectorName;
    vector<string>              tensorName;

    vector<NodeVariableType>    varInNodeVec;
    vector<string>              varNameInNodeVec;
    vector<int>                 varCpntInNodeVec;

    vector<ElementVariableType> varInElmtVec;
    vector<string>              varNameInElmtVec;
    vector<int>                 varCpntInElmtVec;
};
struct InfoForHisOutput{
    vector<HistoryVariableType> hisVarInNodeVec;
    vector<NodeVariableType>    varInNodeVec;
    vector<string>              varNameInNodeVec;
    vector<int>                 varCpntInNodeVec;
    vector<int>                 varCpntIndInNodeVec;    /**< -1 mean all component*/
    vector<int>                 intervalInNodeVec;
    vector<string>              setNameInNodeVec;
    vector<HistoryOutputFormat> outFormatInNodeVec;
    vector<VarOutputForm>       outputFormInNodeVec;     
    vector<PetscScalar *>       bufferInNodeVec;
    vector<PetscScalar *>       timeBuffInNodeVec;
    vector<int>                 bufferFrameNumInNodeVec;
    vector<int>                 dataNumPerFrameInNodeVec;
    vector<int>                 mCpntPerDataInNodeVec;
    vector<string>              fileNameInNodeVec;
 
    vector<HistoryVariableType> hisVarInElmtVec;
    vector<ElementVariableType> varInElmtVec;
    vector<string>              varNameInElmtVec;
    vector<int>                 varCpntInElmtVec;   
    vector<int>                 varCpntIndInElmtVec;    /**< -1 mean all component*/
    vector<int>                 intervalInElmtVec;
    vector<string>              setNameInElmtVec;
    vector<HistoryOutputFormat> outFormatInElmtVec;
    vector<VarOutputForm>       outputFormInElmtVec;    
    vector<PetscScalar *>       bufferInElmtVec;
    vector<PetscScalar *>       timeBuffInElmtVec;
    vector<int>                 bufferFrameNumInElmtVec;
    vector<int>                 dataNumPerFrameInElmtVec;
    vector<int>                 mCpntPerDataInElmtVec;
    vector<string>              fileNameInElmtVec;
};
class PostProcessSystem
{
protected:
    bool            m_ifOutputDesSet;
    bool            m_ifMeshSysSet; 
    bool            m_ifElmtSysSet;
    bool            m_ifLoadCtrlSet;
    bool            m_ifProjVec;            /**< if the Vec m_proj_vec has content*/
    bool            m_ifHisNodeVec;         /**< if the Vec m_his_node_vec has content*/
    bool            m_ifHisElmtVec;         /**< if the Vec m_his_elmt_vec has content*/
    bool            m_ifDmInit;
    InfoForFieldOutput  m_infoForFieldOut;
    InfoForHisOutput    m_infoForHisOut;
    MeshSystem      *m_meshSysPtr;          /**< ptr to relied mesh system*/
    ElementSystem   *m_elmtSysPtr;          /**< ptr to relied elemt system*/
    LoadController  *m_loadCtrlPtr;         /**< ptr to relied load controller*/
    int             m_mCpnt;                /**< num of component of current projected vector*/
    OutputDescription *m_outputDesPtr;      /**< ptr to output description*/
    string          m_prefix;               /**< output file's prefix*/
    static const int m_hisBuffLen=20;       /**< historic variable buffer length*/
protected:
    string fieldOutputFileName(int t_increI,FieldOutputFormat format);
    string hisOutputFileName(string setName,HistoryVariableType vType, HistoryOutputFormat format);
    PetscErrorCode readOutputDes(OutputDescription *t_outputDesPtr);
public:
    PostProcessSystem(OutputDescription *t_outputDesPtr);
    PostProcessSystem(OutputDescription *t_outputDesPtr, MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr, LoadController *t_loadCtrlPtr);
    virtual ~PostProcessSystem();
    virtual PetscErrorCode clear()=0;
    virtual PetscErrorCode init(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr)=0;
    virtual PetscErrorCode init()=0;
    virtual PetscErrorCode checkInit()=0;
    virtual PetscErrorCode projElmtVariable(ElementVariableType varType)=0;
    virtual PetscErrorCode projVecClean()=0;  
    /**
     * calulate required nodal variable and let m_array_his_node ptr to the variable Vec
     * @param varType > NodeVariableType
    */
    virtual PetscErrorCode genNodeVariable(NodeVariableType varType)=0;
    /**
     * restore required nodal variable and let m_array_his_node ptr to null
    */
    virtual PetscErrorCode restoreNodeVariable(NodeVariableType varType)=0;
    /**
     * calulate required elemental variable and let m_array_his_elmt ptr to the variable Vec
    */
    // virtual PetscErrorCode genElmtVariable(ElementVariableType varType)=0;
    virtual PetscErrorCode openNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode)=0;
    virtual PetscErrorCode closeNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode)=0;
    /**
     * output the calculation result (including postprocess to get the output variable)
     * @param t_increI >  increment ID of the latest converged
     * @param t_t > accumulative time of thet latest coverged
    */
    virtual PetscErrorCode output(int t_increI, PetscScalar t_t)=0;
public:
    PetscMPIInt     m_rank;
    PetscMPIInt     m_rankNum;
    Vec             m_proj_weight;          /**< temporary global Vec storing sum of current projected weight*/
    Vec             m_proj_vec;             /**< temporary global Vec storing weighted Vec vals*/
    Vec             m_proj_weight_local;    /**< temporary local Vec storing sum of current projected weight*/
    Vec             m_proj_val_local;       /**< temporary local Vec storing weighted Vec vals*/
    Vec             m_his_node_vec;         /**< temporary global Vec storing historic variable in node position*/
    Vec             m_his_node_vec_local;   /**< not neccessary for m_his_node_vec*/
    Vec             m_his_elmt_vec;         /**< temporary global Vec storing historic variable in elmt position*/
    Vec             m_his_elmt_vec_local;   /**< not neccessary for m_his_elmt_vec*/
    PetscScalar     ***m_array_proj_weight;
    PetscScalar     ***m_array_proj_val;
    PetscScalar     ***m_array_his_node;
    PetscScalar     ***m_array_his_elmt;
};
