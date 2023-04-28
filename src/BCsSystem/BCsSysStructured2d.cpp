#include "BCsSystem/BCsSysStructured2d.h"
PetscErrorCode BCsSysStructured2d::init(){
    if(m_ifHasReadBCDes){
        if(m_ifHasInit)return 0;
        for(vector<SingleBCDes>::iterator itBCBlock=m_bcDesPtr->begin();
        itBCBlock!=m_bcDesPtr->end();++itBCBlock){// loop over every single bcs description block
            vector<PetscInt> &set2Cstd=m_meshPtr->m_setManager.getSet(itBCBlock->s_setName,SetType::NODE);
            for(vector<PetscInt>::iterator itNodeI=set2Cstd.begin();
            itNodeI!=set2Cstd.end();++itNodeI){// loop over every constrained node in this bcs block
                PetscInt nodeGId=m_meshPtr->m_node_gId[*itNodeI];
                for(vector<int>::iterator itDofI=itBCBlock->s_presetDofIds.begin();
                itDofI!=itBCBlock->s_presetDofIds.end();++itDofI){// loop over every constrained dof in a constrained node
                    PetscInt dofGId=nodeGId*m_meshPtr->m_mDof_node+*itDofI;/**< constrained dof's global id*/\
                    PetscScalar presetval=itBCBlock->s_bcVals;
                    map<PetscInt,PetscScalar>::iterator itBcMap=m_bcValsMap.find(dofGId);
                    if(itBcMap==m_bcValsMap.end()){
                        MessagePrinter::printErrorTxt("dof "+to_string(dofGId)+" is repeatedly constrained!");
                        MessagePrinter::exitcfem();
                    }
                    else{
                        m_bcValsMap.insert(pair<PetscInt,PetscScalar>(dofGId,presetval));
                    }
                }
            }
        }
    }
    else{
        MessagePrinter::printErrorTxt("boundary condition description has not been read yet.");
        MessagePrinter::exitcfem();
    }
    PetscCall(DMCreateGlobalVector(m_meshPtr->m_dm,&m_u));
    PetscCall(DMCreateGlobalVector(m_meshPtr->m_dm,&m_b));
    m_ifHasInit=true;
    return 0;
}
PetscErrorCode BCsSysStructured2d::init(BCDescription *t_bcDesPtr){
    if(m_ifHasInit)return 0;
    if(!m_ifHasReadBCDes){
        m_bcDesPtr=t_bcDesPtr;
        m_ifHasReadBCDes=true;
    }
    init();
    return 0;
}
BCsSysStructured2d::~BCsSysStructured2d(){
    VecDestroy(&m_u);
    VecDestroy(&m_b);  
}
PetscErrorCode BCsSysStructured2d::applyBoundaryCondition(PetscScalar facInc,Vec *residualPtr, Mat *AMatrixPtr){
    PetscInt mConstrainedDof=m_bcValsMap.size();    /**< num of constrained dof in this rank*/
    PetscScalar pivot=1.0;
    PetscCall(VecZeroEntries(m_u));   PetscCall(VecZeroEntries(m_b));
    PetscInt *arrayConstrainedRows=new PetscInt[mConstrainedDof];
    PetscScalar *arrayPresetVals=new PetscScalar[mConstrainedDof];
    PetscScalar *arrayZero=nullptr;
    if(facInc!=0.0)
        arrayZero=new PetscScalar[mConstrainedDof];
    /**fill arrayConstrainedRows and arrayPresetVals*/
    /************************************************/
    PetscInt mDof=0;
    for(map<PetscInt,PetscScalar>::iterator itBCMap=m_bcValsMap.begin();
    itBCMap!=m_bcValsMap.end();++itBCMap){//loop over every constrained dof in this rank
        arrayConstrainedRows[mDof]=itBCMap->first;
        if(facInc!=0.0) arrayZero[mDof]=0.0;
        arrayPresetVals[mDof++]=-facInc*itBCMap->second;    /**< because residual=-RHS, so set to negative*/
    }
    /**set prescribed dof in Vec m_u*/
    PetscCall(VecSetValues(m_u,mDof,arrayConstrainedRows,arrayPresetVals,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(m_u));
    if(facInc!=0.0){
        PetscCall(VecSetValues(*residualPtr,mDof,arrayConstrainedRows,arrayZero,INSERT_VALUES));       
    }
    else{
        PetscCall(VecSetValues(*residualPtr,mDof,arrayConstrainedRows,arrayPresetVals,INSERT_VALUES));
    }
    PetscCall(VecAssemblyEnd(m_u));
    PetscCall(MatZeroRowsColumns(*AMatrixPtr,mDof,arrayConstrainedRows,pivot,m_u,m_b));
    if(facInc!=0.0)
        VecAYPX(*residualPtr,1.0,m_b); 
    delete[] arrayConstrainedRows;
    delete[] arrayPresetVals;
    if(facInc!=0.0) delete[] arrayZero;
    return 0;
}