#include"MathUtils/Vector2d.h"
#include"MathUtils/Rank2Tensor2d.h"
#include"MathUtils/ViogtRank2Tensor2D.h"
#include"MathUtils/ViogtRank4Tensor2D.h"
#include "Utils/MessagePrinter.h"
#include"petsc.h"
int main(int argc,char **argv){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&argc,&argv,NULL,NULL);if (ierr) return ierr;
    /**
     * ViogtRank2Tensor2D's verification
    */
    // verificationm of parse
    {
    MessagePrinter::printStars(MessageColor::BLUE);
    MessagePrinter::printTxt("verificationm of parse function of Vector2D, Rank2Tensor2d, ViogtRank2Tensor2D, ViogtRank4Tensor2D",MessageColor::BLUE);
    Vector2d v2d(Vector2d::InitMethod::RANDOM);
    Vector3d v3d=v2d.toVector3d();
    MessagePrinter::printTxt("verification of parse Vector2d 2 Vector3d (if double[] a deep copy)",MessageColor::YELLOW);
    MessagePrinter::printTxt("the original Vector2D:");
    v2d.print();
    MessagePrinter::printTxt("Vector2d parsed Vector3d:");
    v3d.print();

    MessagePrinter::printTxt("verification of parse Rank2Tensor2d 2 Rank2Tensor (if vector<double> a deep copy)",MessageColor::YELLOW);
    Rank2Tensor2d rank2Tensor2d(Rank2Tensor2d::InitMethod::RANDOM);
    Rank2Tensor rank2Tensor=rank2Tensor2d.toRank2Tensor3d();
    MessagePrinter::printTxt("the original Rank2Tensor2d:");
    rank2Tensor2d.print();
    MessagePrinter::printTxt("Rank2Tensor2d parsed Rank2Tensor3d:");
    rank2Tensor.print();

    MessagePrinter::printTxt("verification of parse ViogtRank2Tensor2D 2 Rank2Tensor2d/Rank2Tensor3d (if double[] a deep copy)",MessageColor::YELLOW);
    ViogtRank2Tensor2D VR22(ViogtRank2Tensor2D::InitMethod::RANDOM);
    Rank2Tensor2d R22=VR22.toRank2Tensor2d();
    Rank2Tensor R23=VR22.toRank2Tensor3d();
    MessagePrinter::printTxt("the original ViogtRank2Tensor2D:");
    VR22.print();
    MessagePrinter::printTxt("ViogtRank2Tensor2D parsed Rank2Tensor2d:");
    R22.print();
    MessagePrinter::printTxt("ViogtRank2Tensor2D parsed Rank2Tensor3d:");
    R23.print();

    MessagePrinter::printTxt("verification of parse ViogtRank4Tensor2D 2 Rank4Tensor (if vector<double> a deep copy)",MessageColor::YELLOW);
    ViogtRank4Tensor2D VR42(ViogtRank4Tensor2D::InitMethod::RANDOM);
    Rank4Tensor R43=VR42.toRank4Tensor();
    MessagePrinter::printTxt("the original ViogtRank4Tensor2D:");
    VR42.print();
    MessagePrinter::printTxt("ViogtRank4Tensor2D parsed Rank4Tensor:");
    ViogtRank4Tensor2D(R43).print();
    }

    // verification of L_ij*R_j, R_i*L_ij   // L: ViogtRank2Tensor2D, R: Vector2d
    {
    MessagePrinter::printStars(MessageColor::BLUE);
    MessagePrinter::printTxt("verification of L_ij*R_j R_i*L_ij (L: ViogtRank2Tensor2D, C,R: Vector2d)",MessageColor::BLUE);
    Vector2d C,R(Vector2d::InitMethod::RANDOM);
    ViogtRank2Tensor2D L(ViogtRank2Tensor2D::InitMethod::RANDOM);
    Rank2Tensor L2=L.toRank2Tensor3d();
    Vector3d R2=R.toVector3d();
    Vector3d C2=L2*R2;
    C=L*R;
    MessagePrinter::printTxt("ViogtRank2Tensor2D L_ij:");
    L.print();
    MessagePrinter::printTxt("Vector2d R_j:");
    R.print();
    MessagePrinter::printTxt("L_ij*R_j:");
    C.print();
    MessagePrinter::printTxt("L_ij*R_j (reference):");
    C2.print();
    MessagePrinter::printTxt("R_i*L_ij:");
    (R*L).print();
    MessagePrinter::printTxt("R_i*L_ij (reference):");
    (R2*L2).print();
    }

    // verification of L_ik*R_kj, R_ik*L_kj   // L: ViogtRank2Tensor2D, R: Rank2Tensor2d
    {
    MessagePrinter::printStars(MessageColor::BLUE);
    MessagePrinter::printTxt("verification of L_ik*R_kj, R_ik*L_kj (L: ViogtRank2Tensor2D, R: Vector2d)",MessageColor::BLUE);
    ViogtRank2Tensor2D L(ViogtRank2Tensor2D::InitMethod::RANDOM);
    Rank2Tensor2d R(Rank2Tensor2d::InitMethod::RANDOM);
    MessagePrinter::printTxt("L_ik*R_kj:");
    (L*R).print();
    MessagePrinter::printTxt("L_ik*R_kj (reference):");
    (L.toRank2Tensor3d()*R.toRank2Tensor3d()).print();
    MessagePrinter::printTxt("R_ik*L_kj:");
    (R*L).print();
    MessagePrinter::printTxt("R_ik*L_kj (reference):");
    (R.toRank2Tensor3d()*L.toRank2Tensor3d()).print();   
    }

    // verification of Rank2Tensor2d::calcEigenValueAndEigenVectors
    {
    MessagePrinter::printTxt("verification of Rank2Tensor2d::calcEigenValueAndEigenVectors");
    Rank2Tensor2d L=ViogtRank2Tensor2D(ViogtRank2Tensor2D::InitMethod::RANDOM).toRank2Tensor2d();
    double eigenval[2];
    Rank2Tensor2d eigenvec;
    L.calcEigenValueAndEigenVectors(eigenval,eigenvec);
    PetscPrintf(MPI_COMM_WORLD,"eigenvalue: %14.6e, %14.6e\n",eigenval[0],eigenval[1]);
    MessagePrinter::printTxt("L*eigenvec1=");
    (L*eigenvec.getIthCol(1)).print();
    MessagePrinter::printTxt("eigenvec1*eigenval1=");
    (eigenval[0]*eigenvec.getIthCol(1)).print();
    MessagePrinter::printTxt("L*eigenvec2=");
    (L*eigenvec.getIthCol(2)).print();
    MessagePrinter::printTxt("eigenvec2*eigenval2=");
    (eigenval[1]*eigenvec.getIthCol(2)).print();
    }

    // verification of L_ijkl*R_kl, R_ij*L_ijkl (L: ViogtRank4Tensor2D, R: ViogtRank2Tensor2D) 
    {
    MessagePrinter::printTxt("verification of L_ijkl*R_kl, R_ij*L_ijkl (L: ViogtRank4Tensor2D, R: ViogtRank2Tensor2D)",MessageColor::BLUE);
    ViogtRank4Tensor2D L(ViogtRank4Tensor2D::InitMethod::RANDOM);
    ViogtRank2Tensor2D R(ViogtRank2Tensor2D::InitMethod::RANDOM);
    MessagePrinter::printTxt("L_ijkl*R_kl:");
    (L*R).print();
    MessagePrinter::printTxt("L_ijkl*R_kl (reference):");
    (L.toRank4Tensor().doubledot(R.toRank2Tensor3d())).print();
    MessagePrinter::printTxt("R_ij*L_ijkl:");
    (R*L).print();
    MessagePrinter::printTxt("R_ij*L_ijkl (reference):");
    (R.toRank2Tensor3d().doubledot(L.toRank4Tensor())).print();
    }

    // verification of L_ijkl*R_kl, R_ij*L_ijkl (L: ViogtRank4Tensor2D, R: Rank2Tensor2D) 
    {
    MessagePrinter::printTxt("verification of L_ijkl*R_kl, R_ij*L_ijkl (L: ViogtRank4Tensor2D, R: Rank2Tensor2D)",MessageColor::BLUE);
    ViogtRank4Tensor2D L(ViogtRank4Tensor2D::InitMethod::RANDOM);
    Rank2Tensor2d R(ViogtRank2Tensor2D::InitMethod::RANDOM);
    MessagePrinter::printTxt("L_ijkl*R_kl:");
    (L*R).print();
    MessagePrinter::printTxt("L_ijkl*R_kl (reference):");
    (L.toRank4Tensor().doubledot(R.toRank2Tensor3d())).print();
    MessagePrinter::printTxt("R_ij*L_ijkl:");
    (R*L).print();
    MessagePrinter::printTxt("R_ij*L_ijkl (reference):");
    (R.toRank2Tensor3d().doubledot(L.toRank4Tensor())).print();
    }
    // verification of L_ijmn*R_mnkl (L,R: ViogtRank4Tensor2D)
    {
    MessagePrinter::printTxt("verification of L_ijmn*R_mnkl (L,R: ViogtRank4Tensor2D)",MessageColor::BLUE);
    ViogtRank4Tensor2D L(ViogtRank4Tensor2D::InitMethod::RANDOM);
    ViogtRank4Tensor2D R(ViogtRank4Tensor2D::InitMethod::RANDOM);
    MessagePrinter::printTxt("L_ijmn*R_mnkl:");
    (L*R).print();
    MessagePrinter::printTxt("L_ijmn*R_mnkl (reference):");
    ViogtRank4Tensor2D(L.toRank4Tensor().doubledot(R.toRank4Tensor())).print();
    }
    ierr=PetscFinalize();CHKERRQ(ierr);
    return 0;
}