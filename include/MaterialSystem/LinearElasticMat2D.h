#pragma once
#include "Material2D.h"
#include "nlohmann/json.hpp"
class LinearElasticMat2D:public Material2D{
    private:
    bool m_ifPropInit;  /**< if the material inited*/
    public:
    LinearElasticMat2D():Material2D(false,0.0),m_ifPropInit(false),
                        m_strain(),m_lame(0.0),m_G(0.0),m_planeState(false){}

    LinearElasticMat2D(bool ifLarge,double t_det_dx0dr):Material2D(ifLarge,t_det_dx0dr),m_ifPropInit(false),
                                m_strain(),m_lame(0.0),m_G(0.0),m_planeState(false){}

    LinearElasticMat2D(bool ifLarge,double t_det_dx0dr,nlohmann::json *t_propPtr);
    /**
     * Init material's property by properties nlohmann::json
     * @param t_det_dx0dr > det(dx0/dr), x0 is elmt's coords in ref config
    */
    virtual void initProperty(nlohmann::json *t_propPtr);
    virtual ~LinearElasticMat2D(){};
    /**
     * update qpoint's material status, including cauchy stress, deformation tensor, ...
     * @param incStrainPtr > ptr to deriv of inc strain (du/dX for small strain, Finc for large strain)
     * @param converged < if update iteration converged.
    */
    virtual void updateMaterialBydudx(void *t_incStrainPtr,bool *t_converged);
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
    virtual void getElementVariable(ElementVariableType elmtVarType,void *elmtVarPtr);
    public:
    ViogtRank2Tensor2D m_strain;        /**< ln(V)=ln(B)/2 (Eulerain logarithmic strain)*/
    ViogtRank2Tensor2D m_strain0;       /**< last converged ln(V)=ln(B)/2 (Eulerain logarithmic strain)*/
    double m_lame,m_G;                  /**< material props*/
    bool m_planeState;                  /**< false for plane strain, true for plane stress*/
};