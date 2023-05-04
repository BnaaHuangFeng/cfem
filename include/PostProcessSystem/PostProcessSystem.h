#include "MeshSystem/MeshSystem.h"
#include "ElementSystem/ElementSystem.h"
class PostProcessSystem
{
protected:
    bool            m_ifMeshSysSet; 
    bool            m_ifElmtSysSet;
    bool            m_ifProjValVec;         /**< if the Vec m_proj_val has content*/
    bool            m_ifDmInit;
    MeshSystem      *m_meshSysPtr;          /**< ptr to relied mesh system*/
    ElementSystem   *m_elmtSysPtr;          /**< ptr to relied elemt system*/
    int             m_mCpnt;                /**< num of component of current projected vector*/
    ElementVariableType m_elmtVarType;      /**< elmt vatiable type of current projected vector*/

public:
    PostProcessSystem();
    PostProcessSystem(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr);
    virtual ~PostProcessSystem();
    virtual PetscErrorCode init(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr)=0;
    virtual PetscErrorCode init()=0;
    virtual PetscErrorCode checkInit()=0;
    virtual PetscErrorCode projElmtVariable(ElementVariableType varType)=0;
    virtual PetscErrorCode projVecClean(int mCpnt)=0;
    virtual PetscErrorCode openNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode)=0;
    virtual PetscErrorCode closeNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode)=0;
public:
    PetscMPIInt     m_rank;
    PetscMPIInt     m_rankNum;
    Vec             m_proj_weight;          /**< temporary global Vec storing sum of current projected weight*/
    Vec             m_proj_val;             /**< temporary global Vec storing weighted Vec vals*/
    Vec             m_proj_weight_local;    /**< temporary local Vec storing sum of current projected weight*/
    Vec             m_proj_val_local;       /**< temporary local Vec storing weighted Vec vals*/
    PetscScalar     ***m_array_proj_weight;
    PetscScalar     ***m_array_proj_val;
};
