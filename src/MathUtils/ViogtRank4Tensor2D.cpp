# include "MathUtils/ViogtRank4Tensor2D.h"
# include "MathUtils/ViogtRank2Tensor2D.h"
# include "MathUtils/Rank4Tensor3d.h"
# include "MathUtils/Rank2Tensor2d.h"
ViogtRank4Tensor2D::ViogtRank4Tensor2D():MatrixXd(3,3){};
ViogtRank4Tensor2D::ViogtRank4Tensor2D(const double &val):MatrixXd(3,3,val){};
ViogtRank4Tensor2D::ViogtRank4Tensor2D(const double *vals):MatrixXd(3,3,vals){};
ViogtRank4Tensor2D::ViogtRank4Tensor2D(const MatrixXd &matrix):MatrixXd(matrix){};
ViogtRank4Tensor2D::ViogtRank4Tensor2D(InitMethod initmethod):MatrixXd(3,3,0.0){
    switch (initmethod)
    {
    case InitMethod::ZERO:
        break;
    case InitMethod::IDENTITY:
        for(int indij=0;indij<NViogt;indij++){
            for(int indkl=0;indkl<NViogt;indkl++){
                (*this)(indij,indkl)=indij==indkl;
            }
        }
        break;
    case InitMethod::RANDOM:
        for(int indij=0;indij<NViogt;indij++){
            for(int indkl=0;indkl<NViogt;indkl++){
                (*this)(indij,indkl)=static_cast<double>(1.0*rand()/RAND_MAX);
            }
        }
        break;
    default:
        break;
    }   
}

double & ViogtRank4Tensor2D::operator()(const int indij,const int indkl){
    return MatrixXd::operator()(indij,indkl);
}

double ViogtRank4Tensor2D::operator()(const int indij,const int indkl)const{
    return MatrixXd::operator()(indij,indkl);
}
double & ViogtRank4Tensor2D::operator()(const int i,const int j,const int k,const int l){
    return MatrixXd::operator()(ij2ind[i][j],ij2ind[k][l]);
}
double ViogtRank4Tensor2D::operator()(const int i,const int j,const int k,const int l)const{
    return MatrixXd::operator()(ij2ind[i][j],ij2ind[k][l]);
};
ViogtRank4Tensor2D::ViogtRank4Tensor2D(Rank4Tensor3d rank4Tensor):MatrixXd(3,3,0.0){
    for(int indij=0;indij<NViogt;indij++){
        for(int indkl=0;indkl<NViogt;indkl++){
            int i,j,k,l;
            i=ind2ij[indij][0];
            j=ind2ij[indij][1];
            k=ind2ij[indkl][0];
            l=ind2ij[indkl][1];
            if(rank4Tensor(i,j,k,l)!=rank4Tensor(j,i,k,l)||rank4Tensor(i,j,k,l)!=rank4Tensor(i,j,l,k)){
                MessagePrinter::printErrorTxt("the input Rank4Tensor3d doesn't have the property of viogt minor symtric.");
                MessagePrinter::exitcfem();
            }
            else{
                (*this)(indij,indkl)=rank4Tensor(i,j,k,l);
            }
        }
    }
}
Rank4Tensor3d ViogtRank4Tensor2D::toRank4Tensor(){
    Rank4Tensor3d tmp(0.0);
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            for(int k=0;k<dim;k++){
                for(int l=0;l<dim;l++){
                    tmp(i,j,k,l)=(*this)(i,j,k,l);
                }
            }
        }
    }
    return tmp;
}
ViogtRank2Tensor2D ViogtRank4Tensor2D::operator*(const ViogtRank2Tensor2D &b)const{
    ViogtRank2Tensor2D tmp(0.0);
    int dim=this->getM();
    for(int indij=0;indij<dim;indij++){
        for(int indkl=0;indkl<dim;indkl++){
            tmp(indij)+=(*this)(indij,indkl)*b(indkl)*(indkl==2?2:1);
        }
    }
    return tmp;
};

ViogtRank2Tensor2D operator*(const ViogtRank2Tensor2D &b,const ViogtRank4Tensor2D &a){
    ViogtRank2Tensor2D tmp(0.0);
    int dim=a.getM();
    for(int indkl=0;indkl<dim;indkl++){
        for(int indij=0;indij<dim;indij++){
            tmp(indkl)+=b(indij)*(indij==2?2:1)*a(indij,indkl);
        }
    }
    return tmp;
}
ViogtRank2Tensor2D ViogtRank4Tensor2D::operator*(const Rank2Tensor2d &b)const{
    ViogtRank2Tensor2D tmp(0.0);
    for(int indij=0;indij<3;indij++){     // indij: viogt index sf 1
        for(int k=0;k<2;k++){
            for(int l=0;l<2;l++){
                int indkl=ij2ind[k][l];
                tmp(indij)+=(*this)(indij,indkl)*b(k,l);
            }
        }
    }
    return tmp;
};
ViogtRank2Tensor2D operator*(const Rank2Tensor2d &L,const ViogtRank4Tensor2D & R){
    ViogtRank2Tensor2D tmp(0.0);
    for(int indkl=0;indkl<3;indkl++){
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                int indij=ViogtRank4Tensor2D::ij2ind[i][j];
                tmp(indkl)+=L(i,j)*R(indij,indkl);
            }
        }
    }
    return tmp;
}
ViogtRank4Tensor2D ViogtRank4Tensor2D::operator*(const ViogtRank4Tensor2D &R){
    ViogtRank4Tensor2D tmp(0.0);
    int Nviogt=this->getM();
    for(int indij=0;indij<Nviogt;indij++){
        for(int indkl=0;indkl<Nviogt;indkl++){
            for(int indmn=0;indmn<Nviogt;indmn++){
                tmp(indij,indkl)+=(*this)(indij,indmn)*R(indmn,indkl)*(indmn==(Nviogt-1)?2:1);
            }
        }
    }
    return tmp;
}
ViogtRank4Tensor2D ViogtRank4Tensor2D::operator+(ViogtRank4Tensor2D &R){
    ViogtRank4Tensor2D tmp=MatrixXd::operator+(R);
    return tmp;
}
ViogtRank4Tensor2D ViogtRank4Tensor2D::operator+(ViogtRank4Tensor2D R){
    ViogtRank4Tensor2D tmp=MatrixXd::operator+(R);
    return tmp;
}
ViogtRank4Tensor2D ViogtRank4Tensor2D::operator-(const ViogtRank4Tensor2D &R){
    ViogtRank4Tensor2D tmp=MatrixXd::operator-(R);
    return tmp;
}
MatrixXd ViogtRank4Tensor2D::toFullMatrix()const{
    const int nInd=4;
    MatrixXd tmp(nInd,nInd,0.0);
    const int ind2viogtInd[4]={0,2,2,1};
    for(int indI=0;indI<nInd;indI++){
        for(int indJ=0;indJ<nInd;indJ++){
            MatrixXd(indI,indJ)=(*this)(ind2viogtInd[indI],ind2viogtInd[indJ]);
        }
    }
    return tmp;
}
void ViogtRank4Tensor2D::print() const{
    PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e***\n",(*this)(0,0),(*this)(0,1),(*this)(0,2));
    PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e***\n",(*this)(1,0),(*this)(1,1),(*this)(1,2));
    PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e***\n",(*this)(2,0),(*this)(2,1),(*this)(2,2));
};
