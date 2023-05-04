#pragma once
#include "nlohmann/json.hpp"
#include "MathUtils/MatrixXd.h"
#include "InputSystem/EnumDataType.h"
class Material{
public:
    /**
    * construction
    */
    Material(){}
    /**
    * destruction
    */
    virtual ~Material(){}
    /**
     * Init material's property by properties nlohmann::json
     * @param t_det_dx0dr > det(dx0/dr), x0 is elmt's coords in ref config
    */
    virtual void initProperty(nlohmann::json *t_propPtr)=0;
    /**
     * update qpoint's material status, including cauchy stress, deformation tensor, ...
     * @param incStrainPtr > ptr to deriv of inc strain (du/dX for small strain, Finc for large strain)
     * @param converged < if update iteration converged.
    */
    virtual void updateMaterialBydudx(void *t_incStrainPtr,bool *t_converged)=0;
    /**
     * get tangent modulus by Finc
     * @param incStrainPtr > ptr to deriv of inc strain (du/dX for small strain, Finc for large strain)
     * @param D < ptr to get the tangent modulus
    */
    virtual void getTangentModulus(void *t_incStrainPtr,void *t_D)=0;
    /**
     * get spatial tangent modulus by Finc for large strain
     * @param incStrainPtr > ptr to deriv of inc strain (duInc/dX for small strain, Finc for large strain)
     * @param a < ptr to get the spatial tangent modulus a_4*4
    */
    virtual void getSpatialTangentModulus(void *t_incStrainPtr,MatrixXd *t_a)=0; 
    /**
     * Get material variable of elmtVarType
     * @param elmtVarType > required elemnt variable's type
     * @param elmtVarPtr < ptr to get the elemnt variable (need to preallocate)
    */
    virtual void getMatVariable(ElementVariableType elmtVarType,void *elmtVarPtr)=0;
    /**
     * Get material variable of elmtVarType in PetscScalar array form
     * tensor of rank 2's vector order: 11 22 33 12 13 23
     * @param elmtVarType > required elemnt variable's type
     * @param elmtVarPtr < ptr to store the elemnt variable (need to preallocate)
    */
    virtual void getMatVariableArray(ElementVariableType elmtVarType,PetscScalar *elmtVarPtr)=0;
};