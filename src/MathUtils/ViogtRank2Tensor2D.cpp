#include "MathUtils/Vector2d.h"
#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank2Tensor3d.h"
#include "MathUtils/ViogtRank2Tensor2D.h"
#include "MathUtils/ViogtRank4Tensor2D.h"
#include "MathUtils/Rank2Tensor2d.h"
#include "MathUtils/MatrixXd.h"
#include "MathUtils/TensorConst2D.h"
#include <cmath>
ViogtRank2Tensor2D::ViogtRank2Tensor2D(InitMethod initmethod):Vector3d(Vector3d::InitMethod::ZERO){
    switch (initmethod)
    {
    case InitMethod::ZERO:
        break;
    case InitMethod::IDENTITY:
        (*this)(0)=1.0;
        (*this)(1)=1.0;
        (*this)(2)=0.0;
        break;
    case InitMethod::RANDOM:
        for(int iViogt=0;iViogt<NViogt;iViogt++)
            (*this)(iViogt)=static_cast<double>(1.0*rand()/RAND_MAX);
        break;
    default:
        break;
    }
}
void ViogtRank2Tensor2D::setFromRank2Tensor2D(const Rank2Tensor2d &R){
    static const double small=1E-5;
    double R01=R(0,1), R10=R(1,0);
    double differ=abs(R01-R10);
    double maxval=std::max(abs(R10),abs(R01));
    if(maxval!=0.0)
        differ/=maxval;
    if(differ>small){
        MessagePrinter::printRankError("can not set ViogtRank2Tensor2D by a unsymmetric Rank2Tensor2d");
        MessagePrinter::exitcfem();
    }
    m_vals[0]=R(0,0);
    m_vals[1]=R(1,1);
    m_vals[2]=R(0,1);
}
double& ViogtRank2Tensor2D::operator()(const int i,const int j){
    if(i<0||i>=2||j<0||j>=2){
    MessagePrinter::printErrorTxt(to_string(i)+" is out of range for Vector2");
    MessagePrinter::exitcfem();
}
    if(i==0&&j==0)return Vector3d::operator()(0);
    else if (i!=j)return Vector3d::operator()(2);
    else return Vector3d::operator()(1);
};

double ViogtRank2Tensor2D::operator()(const int i,const int j)const{
    if(i<0||i>=2||j<0||j>=2){
    MessagePrinter::printErrorTxt(to_string(i)+" is out of range for Vector2");
    MessagePrinter::exitcfem();
}
    if(i==0&&j==0)return Vector3d::operator()(0);
    else if (i!=j)return Vector3d::operator()(2);
    else return Vector3d::operator()(1);
};

Vector2d ViogtRank2Tensor2D::operator*(const Vector2d &a)const{
    Vector2d tmp;
    tmp(0)=Vector3d::operator()(0)*a(0)+Vector3d::operator()(2)*a(1);
    tmp(1)=Vector3d::operator()(2)*a(0)+Vector3d::operator()(1)*a(1);
    return tmp;
}
Vector2d operator*(const Vector2d &a,const ViogtRank2Tensor2D &b){
    return b*a;
}
Rank2Tensor2d ViogtRank2Tensor2D::operator*(const Rank2Tensor2d &R){
    Rank2Tensor2d tmp(Rank2Tensor2d::InitMethod::ZERO);
    for(int i=0;i<this->dim;i++){
        for(int j=0;j<this->dim;j++){
            for(int k=0;k<this->dim;k++){
                tmp(i,j)+=(*this)(i,k)*R(k,j);
            }
        }
    }
    return tmp;
}
ViogtRank2Tensor2D ViogtRank2Tensor2D::operator/(const double R){
    ViogtRank2Tensor2D tmp;
    tmp(0)=(*this)(0)/R;
    tmp(1)=(*this)(1)/R;
    tmp(2)=(*this)(2)/R;
    return tmp;
}
/**
 * operation %: tensor product:a_ijkl=L_ij*R_Kl
*/
ViogtRank4Tensor2D ViogtRank2Tensor2D::operator%(const ViogtRank2Tensor2D &R){
    ViogtRank4Tensor2D tmp(ViogtRank4Tensor2D::InitMethod::ZERO);
    for(int viogtI=0;viogtI<NViogt;viogtI++){
        for(int viogtJ=0;viogtJ<NViogt;viogtJ++){
            tmp(viogtI,viogtJ)=(*this)(viogtI)*R(viogtJ);
        }
    }
    return tmp;
}
Rank2Tensor2d ViogtRank2Tensor2D::toRank2Tensor2d()const{
    Rank2Tensor2d tmp(Rank2Tensor2d::InitMethod::ZERO);
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            tmp(i,j)=(*this)(i,j);
        }
    }
    return tmp;
}
Rank2Tensor3d ViogtRank2Tensor2D::toRank2Tensor3d()const{
    Rank2Tensor3d tmp(Rank2Tensor3d::InitMethod::ZERO);
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            tmp(i,j)=(*this)(i,j);
        }
    }
    return tmp;
}

//**************************************************************
//*** For eigen value and eigen vectors and other 
//*** stress and strain decomposition related functions
//**************************************************************
void ViogtRank2Tensor2D::calcEigenValueAndEigenVectors(double eigvalPtr[2],Vector2d eigvecPtr[2]){
    Eigen::Matrix2d _M;

    _M<<(*this)(0,0),(*this)(0,1),
        (*this)(1,0),(*this)(1,1);
    Eigen::EigenSolver<Eigen::Matrix2d> m_eigen_solver;                        /**< solver for eigen solve*/
    m_eigen_solver.compute(_M);
    //
    eigvalPtr[0]=m_eigen_solver.eigenvalues()(0).real();
    eigvalPtr[1]=m_eigen_solver.eigenvalues()(1).real();
    for(int i=0;i<dim;i++){
        eigvecPtr[i](0)=m_eigen_solver.eigenvectors()(0,i).real();
        eigvecPtr[i](1)=m_eigen_solver.eigenvectors()(1,i).real();
    }
}
void ViogtRank2Tensor2D::spectralDecomposition(double eigvalPtr[2],ViogtRank2Tensor2D eigProjPtr[2],bool &repeated){
    static const double small=1E-5;
    Vector2d eigenvec[2];
    calcEigenValueAndEigenVectors(eigvalPtr,eigenvec);
    double differ=abs(eigvalPtr[0]-eigvalPtr[1]);
    double maxEigen=std::max(abs(eigvalPtr[0]),abs(eigvalPtr[1]));
    if(maxEigen!=0.0)
        differ=differ/maxEigen;
    repeated=differ<small;
    for(int eigI=0;eigI<dim;eigI++){
        eigProjPtr[eigI](0)=eigenvec[eigI](0)*eigenvec[eigI](0);
        eigProjPtr[eigI](1)=eigenvec[eigI](1)*eigenvec[eigI](1);
        eigProjPtr[eigI](2)=eigenvec[eigI](0)*eigenvec[eigI](1);
    }
}
ViogtRank2Tensor2D ViogtRank2Tensor2D::isotropicFunc(double (*func)(double)){
    double eigvalPtr[2];
    ViogtRank2Tensor2D eigProjPtr[2];
    bool repeated;
    spectralDecomposition(eigvalPtr,eigProjPtr,repeated);
    double y0=(*func)(eigvalPtr[0]);
    double y1=(*func)(eigvalPtr[1]);
    return eigProjPtr[0]*y0+eigProjPtr[1]*y1;
}
ViogtRank4Tensor2D ViogtRank2Tensor2D::iostropicFuncDeriv(double (*func)(double),double (*funcDeriv)(double)){
    double aEigVal[2]={0.0,0.0};
    double aEigValFunc[2]={0.0,0.0};
    double aEigValFuncDeri[2]={0.0,0.0};
    ViogtRank2Tensor2D aEigProj[2];
    bool ifRepeated=false;
    this->spectralDecomposition(aEigVal,aEigProj,ifRepeated);
    /** Cal eigval and eigvalderiv of y*/
    for(int dimI=0;dimI<2;dimI++){
        aEigValFunc[dimI]=(*func)(aEigVal[dimI]);
        aEigValFuncDeri[dimI]=(*funcDeriv)(aEigVal[dimI]);
    }
    ViogtRank4Tensor2D dydx(ViogtRank4Tensor2D::InitMethod::ZERO);
    if(ifRepeated){// for this tensor's eigenvalue are repeated condition
        dydx=TensorConst2D::IISym*aEigValFuncDeri[0];   // aEigValFuncDeri[0]=aEigValFuncDeri[1] in this conditon
    }
    else{ // for this tensor's eigenvalue are different condition
        ViogtRank4Tensor2D E1XE1=aEigProj[0]%aEigProj[0];
        ViogtRank4Tensor2D E2XE2=aEigProj[1]%aEigProj[1];
        dydx=(TensorConst2D::IISym-E1XE1-E2XE2)*((aEigValFunc[0]-aEigValFunc[1])/(aEigVal[0]-aEigVal[1]))+
            E1XE1*aEigValFuncDeri[0]+E2XE2*aEigValFuncDeri[1];
    }
    return dydx;
}
MatrixXd ViogtRank2Tensor2D::ijkl(const ViogtRank2Tensor2D &R)const{
    const int nMatrix=4;
    int i,j,k,l;
    MatrixXd tmp(nMatrix,nMatrix,0.0);
    for(int matI=0;matI<nMatrix;matI++){
        for(int matJ=0;matJ<nMatrix;matJ++){
            i=matInd2ij[matI][0]; j=matInd2ij[matI][1];   // matI -> i,j
            k=matInd2ij[matJ][0]; l=matInd2ij[matJ][1];   // matJ -> k,l
            tmp(matI,matJ)=(*this)(i,j)*R(k,l);
        }
    }
    return tmp;    
}  
MatrixXd ViogtRank2Tensor2D::iljk(const ViogtRank2Tensor2D &R)const{
    const int nMatrix=4;
    int i,j,k,l;
    MatrixXd tmp(nMatrix,nMatrix,0.0);
    for(int matI=0;matI<nMatrix;matI++){
        for(int matJ=0;matJ<nMatrix;matJ++){
            i=matInd2ij[matI][0]; j=matInd2ij[matI][1];   // matI -> i,j
            k=matInd2ij[matJ][0]; l=matInd2ij[matJ][1];   // matJ -> k,l
            tmp(matI,matJ)=(*this)(i,l)*R(j,k);
        }
    }
    return tmp;
}
MatrixXd ViogtRank2Tensor2D::ikjl(const ViogtRank2Tensor2D &R)const{
    const int nMatrix=4;
    int i,j,k,l;
    MatrixXd tmp(nMatrix,nMatrix,0.0);
    for(int matI=0;matI<nMatrix;matI++){
        for(int matJ=0;matJ<nMatrix;matJ++){
            i=matInd2ij[matI][0]; j=matInd2ij[matI][1];   // matI -> i,j
            k=matInd2ij[matJ][0]; l=matInd2ij[matJ][1];   // matJ -> k,l
            tmp(matI,matJ)=(*this)(i,k)*R(j,l);
        }
    }
    return tmp;
}