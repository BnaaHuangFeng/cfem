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
public:
    PostStructured2d();
    PostStructured2d(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr);
    virtual ~PostStructured2d();
    virtual PetscErrorCode init(MeshSystem *t_meshSysPtr, ElementSystem *t_elmtSysPtr);
    virtual PetscErrorCode init();
    virtual PetscErrorCode checkInit();
    virtual PetscErrorCode projElmtVariable(ElementVariableType varType);
    virtual PetscErrorCode projVecClean(int mCpnt);
    virtual PetscErrorCode openNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode);
    virtual PetscErrorCode closeNodeVariableVec(Vec *globalVecPtr, Vec *localVecPtr, PetscScalar ****arrayPtrRef, PetscInt mCpnt, VecAccessMode mode);
};
