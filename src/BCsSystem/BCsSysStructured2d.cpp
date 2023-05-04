#include "BCsSystem/BCsSysStructured2d.h"
BCsSysStructured2d::~BCsSysStructured2d(){
    if(m_arrayConstrainedRows)delete[] m_arrayConstrainedRows;
    if(m_arrayPresetVals)delete[] m_arrayPresetVals;
    if(m_arrayZero)delete[] m_arrayZero;
    m_arrayConstrainedRows=nullptr;
    m_arrayPresetVals=nullptr;
    m_arrayZero=nullptr;
    VecDestroy(&m_uIncInitial);
    if(m_drclt_method==DirichletMethod::SETLARGE)
        VecDestroy(&m_max_entry_vec);
}
PetscErrorCode BCsSysStructured2d::init(){
    if(m_ifHasReadBCDes){
        if(!m_ifSetMeshSysPtr){
            MessagePrinter::printErrorTxt("BCsSystem: its relied mesh system has not yet been set.");
            MessagePrinter::exitcfem();             
        }
        if(m_ifHasInit)return 0;
        // init m_bcValsMap *************//
        // ******************************//
        for(vector<SingleBCDes>::iterator itBCBlock=m_bcDesPtr->begin();
        itBCBlock!=m_bcDesPtr->end();++itBCBlock){// loop over every single bcs description block
            vector<PetscInt> &set2Cstd=m_meshSysPtr->m_setManager.getSet(itBCBlock->s_setName,SetType::NODE);
            for(vector<PetscInt>::iterator itNodeI=set2Cstd.begin();
            itNodeI!=set2Cstd.end();++itNodeI){// loop over every constrained node in this bcs block
                PetscInt nodeGId=m_meshSysPtr->m_node_gId[*itNodeI];
                for(vector<int>::iterator itDofI=itBCBlock->s_presetDofIds.begin();
                itDofI!=itBCBlock->s_presetDofIds.end();++itDofI){// loop over every constrained dof in a constrained node
                    PetscInt dofGId=nodeGId*m_meshSysPtr->m_mDof_node+*itDofI-1;/**< constrained dof's global id*/
                    PetscScalar presetval=itBCBlock->s_bcVals;
                    map<PetscInt,PetscScalar>::iterator itBcMap=m_bcValsMap.find(dofGId);
                    if(itBcMap!=m_bcValsMap.end()){
                        MessagePrinter::printErrorTxt("dof "+to_string(dofGId)+" is repeatedly constrained!");
                        MessagePrinter::exitcfem();
                    }
                    else{
                        m_bcValsMap.insert(pair<PetscInt,PetscScalar>(dofGId,presetval));
                    }
                }
            }
        }
        m_mConstrainedDof=m_bcValsMap.size();
        // init m_arrayConstrainedRows, m_arrayPresetVals, m_arrayZero
        m_arrayConstrainedRows=new PetscInt[m_mConstrainedDof];
        m_arrayPresetVals=new PetscScalar[m_mConstrainedDof];
        m_arrayZero=new PetscScalar[m_mConstrainedDof];
        m_penalty=0;
        PetscInt dofI=0;
        for(map<PetscInt,PetscScalar>::iterator itBCMap=m_bcValsMap.begin();
        itBCMap!=m_bcValsMap.end();++itBCMap){//loop over every constrained dof in this rank
            m_arrayConstrainedRows[dofI]=itBCMap->first;
            m_arrayZero[dofI]=0.0;
            m_arrayPresetVals[dofI++]=itBCMap->second;
        }
        if(dofI!=m_mConstrainedDof){
            MessagePrinter::printErrorTxt("constrained dof num comes from m_bcValsMap.size() is not equal to that comes from iterator.");
        }
    }
    else{
        MessagePrinter::printErrorTxt("boundary condition description has not been read yet.");
        MessagePrinter::exitcfem();
    }
    PetscCall(DMCreateGlobalVector(m_meshSysPtr->m_dm,&m_uIncInitial));
    if(m_drclt_method==DirichletMethod::SETLARGE){
        PetscCall(DMCreateGlobalVector(m_meshSysPtr->m_dm,&m_max_entry_vec));
    }
    m_ifHasInit=true;
    return 0;
}
PetscErrorCode BCsSysStructured2d::init(BCDescription *t_bcDesPtr){
    if(!m_ifSetMeshSysPtr){
        MessagePrinter::printErrorTxt("BCsSystem: its relied mesh system has not yet been set.");
        MessagePrinter::exitcfem();             
    }
    if(m_ifHasInit)return 0;
    if(!m_ifHasReadBCDes){
        m_bcDesPtr=t_bcDesPtr;
        m_ifHasReadBCDes=true;
    }
    init();
    return 0;
}
PetscErrorCode BCsSysStructured2d::init(BCDescription *t_bcDesPtr, MeshSystem *t_meshSysPtr){
    setMeshSysPtr(t_meshSysPtr);
    init(t_bcDesPtr);
    return 0;
}
PetscErrorCode BCsSysStructured2d::checkInit(){
    if(!m_ifSetMeshSysPtr){
        MessagePrinter::printErrorTxt("BCsSystem: its relied mesh system has not yet been set.");
        MessagePrinter::exitcfem();    
    }
    if(!m_ifHasReadBCDes){
        MessagePrinter::printErrorTxt("BCsSystem: boundary condition description has not been read yet.");
        MessagePrinter::exitcfem();        
    }
    if(!m_ifHasInit){
        MessagePrinter::printErrorTxt("BCsSystem: boundary condition has not yet been analysed and init.");
        MessagePrinter::exitcfem();              
    }
    return 0;
}
PetscErrorCode BCsSysStructured2d::setArrayPresetVals(PetscScalar facInc){
    PetscInt dofI=0;
    for(map<PetscInt,PetscScalar>::iterator itBCMap=m_bcValsMap.begin();
    itBCMap!=m_bcValsMap.end();++itBCMap){//loop over every constrained dof in this rank
        m_arrayPresetVals[dofI++]=facInc*itBCMap->second;    /**< because residual=-RHS, so set to negative*/
    }
    if(dofI!=m_mConstrainedDof){
        MessagePrinter::printErrorTxt("constrained dof num comes from m_bcValsMap.size() is not equal to that comes from iterator.");
    }    
    return 0;
}
PetscErrorCode BCsSysStructured2d::setInitialSolution(PetscScalar facInc){
    setArrayPresetVals(facInc);
    PetscCall(VecZeroEntries(m_uIncInitial));
    PetscCall(VecSetValues(m_uIncInitial,m_mConstrainedDof,m_arrayConstrainedRows,m_arrayPresetVals,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(m_uIncInitial));
    PetscCall(VecAssemblyEnd(m_uIncInitial)); 
    
    return 0;
}
PetscErrorCode BCsSysStructured2d::applyResidualBoundaryCondition(Vec *residualPtr){
    PetscCall(VecSetValues(*residualPtr,m_mConstrainedDof,m_arrayConstrainedRows,m_arrayZero,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(*residualPtr));
    PetscCall(VecAssemblyEnd(*residualPtr));
    return 0;
}
PetscErrorCode BCsSysStructured2d::applyJacobianBoundaryCondition(Mat *AMatrixPtr){
    if(m_drclt_method==DirichletMethod::SETUNIT){
        const PetscScalar pivot=1.0;
        PetscCall(MatZeroRowsColumns(*AMatrixPtr,m_mConstrainedDof,m_arrayConstrainedRows,pivot,NULL,NULL));
    }
    else if(m_drclt_method==DirichletMethod::SETLARGE){
        for(PetscInt dofI=0;dofI<m_mConstrainedDof;++dofI){
            PetscCall(MatSetValue(*AMatrixPtr,m_arrayConstrainedRows[dofI],m_arrayConstrainedRows[dofI],m_penalty,INSERT_VALUES));
        }
        // delete?
        PetscCall(MatAssemblyBegin(*AMatrixPtr,MAT_FINAL_ASSEMBLY));
        PetscCall(MatAssemblyEnd(*AMatrixPtr,MAT_FINAL_ASSEMBLY));
    }
    return 0;
}
// PetscErrorCode BCsSysStructured2d::applyBoundaryCondition(PetscScalar facInc,Vec *residualPtr, Mat *AMatrixPtr){
//     PetscInt mConstrainedDof=m_bcValsMap.size();    /**< num of constrained dof in this rank*/
//     PetscScalar pivot=1.0;
//     PetscCall(VecZeroEntries(m_u));   PetscCall(VecZeroEntries(m_b));
//     PetscInt *arrayConstrainedRows=new PetscInt[mConstrainedDof];
//     PetscScalar *arrayPresetVals=new PetscScalar[mConstrainedDof];
//     PetscScalar *arrayZero=nullptr;
//     if(facInc!=0.0)
//         arrayZero=new PetscScalar[mConstrainedDof];
//     /**fill arrayConstrainedRows and arrayPresetVals*/
//     /************************************************/
//     PetscInt mDof=0;
//     for(map<PetscInt,PetscScalar>::iterator itBCMap=m_bcValsMap.begin();
//     itBCMap!=m_bcValsMap.end();++itBCMap){//loop over every constrained dof in this rank
//         arrayConstrainedRows[mDof]=itBCMap->first;
//         if(facInc!=0.0) arrayZero[mDof]=0.0;
//         arrayPresetVals[mDof++]=-facInc*itBCMap->second;    /**< because residual=-RHS, so set to negative*/
//     }
//     /**set prescribed dof in Vec m_u*/
//     PetscCall(VecSetValues(m_u,mDof,arrayConstrainedRows,arrayPresetVals,INSERT_VALUES));
//     PetscCall(VecAssemblyBegin(m_u));
//     if(facInc!=0.0){
//         PetscCall(VecSetValues(*residualPtr,mDof,arrayConstrainedRows,arrayZero,INSERT_VALUES));       
//     }
//     else{
//         PetscCall(VecSetValues(*residualPtr,mDof,arrayConstrainedRows,arrayPresetVals,INSERT_VALUES));
//     }
//     PetscCall(VecAssemblyEnd(m_u));
//     PetscCall(MatZeroRowsColumns(*AMatrixPtr,mDof,arrayConstrainedRows,pivot,m_u,m_b));
//     if(facInc!=0.0)
//         VecAYPX(*residualPtr,1.0,m_b); 
//     delete[] arrayConstrainedRows;
//     delete[] arrayPresetVals;
//     if(facInc!=0.0) delete[] arrayZero;
//     return 0;
// }
PetscErrorCode BCsSysStructured2d::update_penalty(Mat *AMatrixPtr){
    if(m_drclt_method!=DirichletMethod::SETLARGE)return 0;
    PetscCall(MatGetRowMaxAbs(*AMatrixPtr,m_max_entry_vec,NULL));
    PetscCall(VecMax(m_max_entry_vec,NULL,&m_penalty));
    m_penalty=m_penalty*m_penalty_coef;
    return 0;
}