#pragma once
#include "MeshSystem/MeshSystem.h"
#include "ElementSystem/ElementSystem.h"
#include "PostProcessSystem/PostProcessSystem.h"
class PostStructured2d:public PostProcessSystem
{
private:
    DM              m_dmScalar;
    DM              m_dmRank2Tensor2d;
    DM              m_dmRank2Tensor3d;
    DMDALocalInfo   m_dmInfo;               /**< DMDALocalInfo of mesh system's DM*/
private:
    PetscErrorCode initDm();
    /**
     * init historic variable buffer
    */
    PetscErrorCode initBuffer();
    /**
     * get the corresponding DM ptr by vector componenet num
     * @param dmPtrAdr < address to store the DM ptr
     * @param mCpnt > vector component num
    */
    PetscErrorCode getDmPtrByCpntNum(DM **dmPtrAdr,int mCpnt);
    /**
     * Get element's global DMDA id via its id in rank
     * @param rId > elmt's id in rank
     * @param xIPtr < ptr to dmda x index
     * @param yIPtr < ptr to dmda y index
    */
    void getElmtDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr);
    /**
     * Get node's global DMDA id via its id in rank
     * @param rId > elmt's id in rank
     * @param xIPtr > ptr to dmda x index
     * @param yIPtr > ptr to dmda y index
    */
    void getNodeDmdaIndByRId(PetscInt rId,PetscInt *xIPtr,PetscInt *yIPtr);
    /**
     * Add a element's Vector to global one (need to do 
     * closeNodeVariableVec after elmts in this rank have called this func in order to complete assembly)
     * @param rId > elmt's id in this rank
     * @param globalArray < ptr to global vector array (y ind in dm, x ind in dm, component id) -> val
     * @param localArray > ptr to elmt local array (node id in a elmt, componenet id) -> val
     * @param mCpnt > vector component num
    */
    PetscErrorCode addElmtVec(PetscInt rId,PetscScalar ***globalArray, PetscScalar **localArray,int mCpnt);
    /**
     * Add a element's residual (unbalanced forces (f^int-f^ext) ) Vector to global one by DMDA index (need to do 
     * closeNodeVariableVec after elmts in this rank have called this func in order to complete assembly)
     * @param xI > DMDA x index
     * @param yI > DMDA y index
     * @param globalArray < ptr to global vector array (y ind in dm, x ind in dm, component id) -> val
     * @param localArray > ptr to elmt local array (node id in a elmt, componenet id) -> val
     * @param mCpnt > vector component num
    */
    PetscErrorCode addElmtVecByDmdaInd(PetscInt xI,PetscInt yI,PetscScalar ***globalArray, PetscScalar **localArray,int mCpnt);  
    PetscErrorCode outputFieldVariable(int t_increI, PetscScalar t_t);
    PetscErrorCode outputHisVariable(int t_increI, PetscScalar t_t); 
    void openOutputFile(string fileName, ofstream *ofPtr,ios_base::openmode mode); 

public:
    PostStructured2d(OutputDescription *t_outputDesPtr);
    PostStructured2d(OutputDescription *t_outputDesPtr,MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr, LoadController *t_loadCtrlPtr);
    virtual ~PostStructured2d();
    PetscErrorCode clear();
    virtual PetscErrorCode init(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr);
    virtual PetscErrorCode init();
    virtual PetscErrorCode checkInit();
    /**
     * output the calculation result (including postprocess to get the output variable)
     * @param t_increI >  increment ID of the latest converged
     * @param t_t > accumulative time of thet latest coverged
    */
    virtual PetscErrorCode output(int t_increI, PetscScalar t_t);
    virtual PetscErrorCode projElmtVariable(ElementVariableType varType);
    virtual PetscErrorCode projVecClean();
    /**
     * calulate required nodal variable and let m_array_his_node ptr to the variable Vec
     * @param varType > NodeVariableType
    */
    virtual PetscErrorCode genNodeVariable(NodeVariableType varType);
    /**
     * restore required nodal variable and let m_array_his_node ptr to null
     * @param varType > NodeVariableType
    */
    virtual PetscErrorCode restoreNodeVariable(NodeVariableType varType);
    /**
     * calulate required elemental variable and let m_array_his_elmt ptr to the variable Vec
    */
    // virtual PetscErrorCode genElmtVariable(ElementVariableType varType);
    virtual PetscErrorCode openNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode);
    virtual PetscErrorCode closeNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode);

};