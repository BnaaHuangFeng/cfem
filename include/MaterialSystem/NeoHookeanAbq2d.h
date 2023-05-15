#pragma once
#include "Material2D.h"
#include "MaterialSystem/ElasticConst.h"
#include "nlohmann/json.hpp"
/**
 * this class implement the Neohooken material of abaqus version
 * elastic energy phi = C10 (I1_ - 3) + 1/D1 (J-1)^2
 *                    = G/2 (I1_ - 3) + K/2 (J-1)^2
*/
class NeoHookeanAbq2d:public Material2D{
    private:
    bool m_ifPropInit;  /**< if the material inited*/
    private:
    /**
     * check if the m_nLarge is true, if not, print error and exit
    */
    PetscErrorCode checkIfLargeStrain();
    public:
    NeoHookeanAbq2d():Material2D(false,0.0),m_ifPropInit(false),
                        m_B(ViogtRank2Tensor2D::InitMethod::IDENTITY),
                        m_B0(ViogtRank2Tensor2D::InitMethod::IDENTITY),
                        m_K(0.0),m_G(0.0),m_planeState(false){checkIfLargeStrain();}

    NeoHookeanAbq2d(bool ifLarge,double t_det_dx0dr):Material2D(ifLarge,t_det_dx0dr),m_ifPropInit(false),
                        m_B(ViogtRank2Tensor2D::InitMethod::IDENTITY),
                        m_B0(ViogtRank2Tensor2D::InitMethod::IDENTITY),
                        m_K(0.0),m_G(0.0),m_planeState(false){checkIfLargeStrain();}

    NeoHookeanAbq2d(bool ifLarge,double t_det_dx0dr,nlohmann::json *t_propPtr);
    /**
     * Init material's property by properties nlohmann::json
     * @param t_det_dx0dr > det(dx0/dr), x0 is elmt's coords in ref config
    */
    virtual void initProperty(nlohmann::json *t_propPtr);
    virtual ~NeoHookeanAbq2d(){};
    /**
     * update qpoint's material status, including cauchy stress, deformation tensor, ...
     * @param incStrainPtr > ptr to deriv of inc strain (du/dX for small strain, Finc for large strain)
     * @param converged < if update iteration converged.
    */
    virtual void updateMaterialBydudx(void *t_incStrainPtr,bool *t_converged);
    virtual void updateConvergence();
    /**
     * get tangent modulus by Finc
     * @param incStrainPtr > ptr to deriv of inc strain (du/dx for small strain, Finc for large strain)
     * @param D < ptr to get the tangent modulus
     * @param nLarge > large strain flag, 0 for small strain, 1 for large stran
    */
    virtual void getTangentModulus(void *t_incStrainPtr,void *t_D);
    /**
     * get spatial tangent modulus by Finc for large strain refer to (14.99) of CMFP (page 598)
     * @param incStrainPtr > ptr to deriv of inc strain (du/dX for small strain, Finc for large strain)
     * @param a < ptr to get the spatial tangent modulus a_4*4
    */
    virtual void getSpatialTangentModulus(void *t_incStrainPtr,MatrixXd *t_a);
    /**
     * Get material variable of elmtVarType
     * @param elmtVarType > required elemnt variable's type
     * @param elmtVarPtr < ptr to store the elemnt variable (need to preallocate)
    */
    virtual void getMatVariable(ElementVariableType elmtVarType,void *elmtVarPtr);
    /**
     * Get material variable of elmtVarType in PetscScalar array form
     * tensor of rank 2's vector order: 11 22 33 12 13 23
     * @param elmtVarType > required elemnt variable's type
     * @param elmtVarPtr < ptr to store the elemnt variable (need to preallocate)
    */
    virtual void getMatVariableArray(ElementVariableType elmtVarType,PetscScalar *elmtVarPtr);
    virtual double getLame(){
        return ElasticConst::getLameByK_G(m_K,m_G);
    }
    virtual double getG(){
        return m_G;
    };
    public:
    ViogtRank2Tensor2D m_B;             /**< left cauchy-green tensor B=F*F^T*/
    ViogtRank2Tensor2D m_B0;            /**< last converged left cauchy-green tensor B0*/
    double             m_T33;           /**< T is kirichorff stress*/
    double m_K,m_G;                     /**< material props*/
    bool m_planeState;                  /**< false for plane strain, true for plane stress*/
};