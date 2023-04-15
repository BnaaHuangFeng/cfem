# include "MathUtils/ViogtRank4Tensor2D.h"
# include "MathUtils/Rank4Tensor.h"
ViogtRank4Tensor2D::ViogtRank4Tensor2D(InitMethod initmethod):MatrixXd(3,3,0.0){
    switch (initmethod)
    {
    case InitMethod::ZERO:
        break;
    case InitMethod::IDENTITY:
        for(int indij=1;indij<=NViogt;indij++){
            for(int indkl=1;indkl<=NViogt;indkl++){
                (*this)(indij,indkl)=indij==indkl;
            }
        }
        break;
    case InitMethod::RANDOM:
        for(int indij=1;indij<=NViogt;indij++){
            for(int indkl=1;indkl<=NViogt;indkl++){
                (*this)(indij,indkl)=static_cast<double>(1.0*rand()/RAND_MAX);
            }
        }
        break;
    default:
        break;
    }   
}
ViogtRank4Tensor2D::ViogtRank4Tensor2D(Rank4Tensor rank4Tensor):MatrixXd(3,3,0.0){
    for(int indij=1;indij<=NViogt;indij++){
        for(int indkl=1;indkl<=NViogt;indkl++){
            int i,j,k,l;
            i=ind2ij[indij-1][0];
            j=ind2ij[indij-1][1];
            k=ind2ij[indkl-1][0];
            l=ind2ij[indkl-1][1];
            if(rank4Tensor(i,j,k,l)!=rank4Tensor(j,i,k,l)||rank4Tensor(i,j,k,l)!=rank4Tensor(i,j,l,k)){
                MessagePrinter::printErrorTxt("the input Rank4Tensor doesn't have the property of viogt minor symtric.");
                MessagePrinter::exitAsFem();
            }
            else{
                (*this)(indij,indkl)=rank4Tensor(i,j,k,l);
            }
        }
    }
}
Rank4Tensor ViogtRank4Tensor2D::toRank4Tensor(){
    Rank4Tensor tmp(0.0);
    for(int i=1;i<=dim;i++){
        for(int j=1;j<=dim;j++){
            for(int k=1;k<=dim;k++){
                for(int l=1;l<=dim;l++){
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
    for(int indij=1;indij<=dim;indij++){
        for(int indkl=1;indkl<=dim;indkl++){
            tmp(indij)+=(*this)(indij,indkl)*b(indkl)*(indkl==3?2:1);
        }
    }
    return tmp;
};

ViogtRank2Tensor2D operator*(const ViogtRank2Tensor2D &b,const ViogtRank4Tensor2D &a){
    ViogtRank2Tensor2D tmp(0.0);
    int dim=a.getM();
    for(int indkl=1;indkl<=dim;indkl++){
        for(int indij=1;indij<=dim;indij++){
            tmp(indkl)+=b(indij)*(indij==3?2:1)*a(indij,indkl);
        }
    }
    return tmp;
}
ViogtRank2Tensor2D ViogtRank4Tensor2D::operator*(const Rank2Tensor2d &b)const{
    ViogtRank2Tensor2D tmp(0.0);
    for(int indij=1;indij<=3;indij++){     // indij: viogt index sf 1
        for(int k=1;k<=2;k++){
            for(int l=1;l<=2;l++){
                int indkl=ij2ind[k-1][l-1];
                tmp(indij)+=(*this)(indij,indkl)*b(k,l);
            }
        }
    }
    return tmp;
};
ViogtRank2Tensor2D operator*(const Rank2Tensor2d &L,const ViogtRank4Tensor2D & R){
    ViogtRank2Tensor2D tmp(0.0);
    for(int indkl=1;indkl<=3;indkl++){
        for(int i=1;i<=2;i++){
            for(int j=1;j<=2;j++){
                int indij=ViogtRank4Tensor2D::ij2ind[i-1][j-1];
                tmp(indkl)+=L(i,j)*R(indij,indkl);
            }
        }
    }
    return tmp;
}
ViogtRank4Tensor2D ViogtRank4Tensor2D::operator*(const ViogtRank4Tensor2D &R){
    ViogtRank4Tensor2D tmp(0.0);
    int Nviogt=this->getM();
    for(int indij=1;indij<=Nviogt;indij++){
        for(int indkl=1;indkl<=Nviogt;indkl++){
            for(int indmn=1;indmn<=Nviogt;indmn++){
                tmp(indij,indkl)+=(*this)(indij,indmn)*R(indmn,indkl)*(indmn==Nviogt?2:1);
            }
        }
    }
    return tmp;
}

void ViogtRank4Tensor2D::print() const{
    PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e***\n",(*this)(1,1),(*this)(1,2),(*this)(1,3));
    PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e***\n",(*this)(2,1),(*this)(2,2),(*this)(2,3));
    PetscPrintf(PETSC_COMM_WORLD,"*** %14.6e ,%14.6e ,%14.6e***\n",(*this)(3,1),(*this)(3,2),(*this)(3,3));
};