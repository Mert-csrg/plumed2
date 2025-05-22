/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/Torsion.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR PUCKERING
/*
Calculate sugar pseudorotation coordinates.

This command can be used to calculate pseudorotations for the rings in sugars (puckers). It works for both
5-membered and 6-membered rings. Notice that there are two different implementations depending if
one passes 5 or 6 atoms in the ATOMS keyword. This input tells plumed to print the puckering phase angle of the
 second nucleotide of a RNA molecule on file COLVAR.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt65/AA.pdb
MOLINFO STRUCTURE=regtest/basic/rt65/AA.pdb MOLTYPE=rna
puck: PUCKERING ATOMS=@sugar-2
PRINT ARG=puck.phs FILE=COLVAR
```

For 5-membered rings the implementation is the one discussed in the first of the papers in the bibliography below.
This implementation is simple and can be used in RNA to distinguish C2'-endo and C3'-endo conformations.
Both the polar coordinates (phs and amp) and the Cartesian coordinates (Zx and Zy) are provided.
C2'-endo conformations have negative Zx, whereas C3'-endo conformations have positive Zy.
The notation is consistent with the notation in that first paper.
The five atoms should be provided as C4',O4',C1',C2',C3'.
Notice that this is the same order that can be obtained using the [MOLINFO](MOLINFO.md) syntax (see example below).

For 6-membered rings the implementation is the general Cremer-Pople one that is discussed in the second and third
papers in the bibliography.
This implementation provides both a triplet with Cartesian components (qx, qy, and qz)
and a triplet of polar components (amplitude, phi, and theta).
Applications of this particular implementation are yet to be published (paper in preparation).

> [!NOTE]
> The 6-membered ring implementation distributed with previous versions of PLUMED lead to
> qx and qy values that had an opposite sign with respect to those originally defined in the
> thid reference in the bibliograph below.  The bug was fixed in version 2.5.

*/
//+ENDPLUMEDOC

class Puckering : public Colvar {
    
public:
    explicit Puckering(const ActionOptions&);
    void calculate() override;
    static void registerKeywords(Keywords& keys);
    void calculate5m();
    void calculate6m();
    void calculate7m();
};

PLUMED_REGISTER_ACTION(Puckering,"PUCKERING")

void Puckering::registerKeywords(Keywords& keys) {
    Colvar::registerKeywords( keys );
    keys.remove("NOPBC");
    keys.add("atoms","ATOMS","the five, six or seven atoms of the sugar ring in the proper order");
    keys.addOutputComponent("phs","default","Pseudorotation phase (5 membered rings)");
    keys.addOutputComponent("amp","default","Pseudorotation amplitude (5 membered rings)");
    keys.addOutputComponent("Zx","default","Pseudorotation x Cartesian component (5 membered rings)");
    keys.addOutputComponent("Zy","default","Pseudorotation y Cartesian component (5 membered rings)");
    keys.addOutputComponent("phi","default","Pseudorotation phase (6 membered rings)");
    keys.addOutputComponent("theta","default","Theta angle (6 membered rings)");
    keys.addOutputComponent("amplitude","default","Pseudorotation amplitude (6 membered rings)");
    keys.addOutputComponent("qx","default","Cartesian component x (6 membered rings)");
    keys.addOutputComponent("qy","default","Cartesian component y (6 membered rings)");
    keys.addOutputComponent("qz","default","Cartesian component z (6 membered rings)");
    keys.addOutputComponent("q2","default","Amplitude 2 (7 membered rings)");
    keys.addOutputComponent("q3","default","Amplitude 3 (7 membered rings)");
    keys.addOutputComponent("phi2","default","Phase 2 (7 membered rings)");
    keys.addOutputComponent("phi3","default","Phase 3 (7 membered rings)");
    keys.addOutputComponent("x7m","default","Hyper-Cartesian component x (7 membered rings)");
    keys.addOutputComponent("y7m","default","Hyper-Cartesian component y (7 membered rings)");
    keys.addOutputComponent("z7m","default","Hyper-Cartesian component z (7 membered rings)");
    keys.addOutputComponent("w7m","default","Hyper-Cartesian component w (7 membered rings)");
    keys.addOutputComponent("A2","default","Cremer-Pople Summation A2 (7 membered rings)");
    keys.addOutputComponent("B2","default","Cremer-Pople Summation B2 (7 membered rings)");
    keys.addOutputComponent("A3","default","Cremer-Pople Summation A3 (7 membered rings)");
    keys.addOutputComponent("B3","default","Cremer-Pople Summation B3 (7 membered rings)");
    keys.addOutputComponent("tau2","default","Transformed amplitude tau2 (nonlinear amplitude for mode 2)");
    keys.addOutputComponent("tau3","default","Transformed amplitude tau3 (nonlinear amplitude for mode 3)");
    keys.addOutputComponent("theta7m","default","Angular deviation angle for 7 membered rings");
    keys.addOutputComponent("rho","default","Total puckering amplitude for 7 membered rings");
    keys.addDOI("10.1021/ct401013s");
    keys.addDOI("10.1021/ja00839a011");
    keys.addDOI("10.1021/ja068411o");
}

Puckering::Puckering(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
    std::vector<AtomNumber> atoms;
    parseAtomList("ATOMS",atoms);
    if(atoms.size()!=5 && atoms.size()!=6 && atoms.size()!=7) error("only for 5, 6 or 7-membered rings");
    checkRead();
    
    if(atoms.size()==5) {
        log.printf("  between atoms %d %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial(),atoms[4].serial());
    } else if(atoms.size()==6) {
        log.printf("  between atoms %d %d %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial(),atoms[4].serial(),atoms[5].serial());
    } else if(atoms.size()==7) {
        log.printf("  between atoms %d %d %d %d %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial(),atoms[4].serial(),atoms[5].serial(),atoms[6].serial());
    }    else error("ATOMS should specify 5, 6 or 7 atoms");
    
    if(atoms.size()==5) {
      addComponentWithDerivatives("phs");
      componentIsPeriodic("phs","-pi","pi");
      addComponentWithDerivatives("amp");
      componentIsNotPeriodic("amp");
      addComponentWithDerivatives("Zx");
      componentIsNotPeriodic("Zx");
      addComponentWithDerivatives("Zy");
      componentIsNotPeriodic("Zy");
    } else if(atoms.size()==6) {
      addComponentWithDerivatives("qx");
      componentIsNotPeriodic("qx");
      addComponentWithDerivatives("qy");
      componentIsNotPeriodic("qy");
      addComponentWithDerivatives("qz");
      componentIsNotPeriodic("qz");
      addComponentWithDerivatives("phi");
      componentIsPeriodic("phi","0","2pi");
      addComponentWithDerivatives("theta");
      componentIsNotPeriodic("theta");
      addComponentWithDerivatives("amplitude");
      componentIsNotPeriodic("amplitude");
    } else if(atoms.size()==7) {
      addComponentWithDerivatives("phi2");
      componentIsPeriodic("phi2","0","2pi");
      addComponentWithDerivatives("phi3");
      componentIsPeriodic("phi3","0","2pi");
      addComponentWithDerivatives("q2");
      componentIsNotPeriodic("q2");
      addComponentWithDerivatives("q3");
      componentIsNotPeriodic("q3");
      addComponentWithDerivatives("x7m");
      componentIsNotPeriodic("x7m");
      addComponentWithDerivatives("y7m");
      componentIsNotPeriodic("y7m");
      addComponentWithDerivatives("z7m");
      componentIsNotPeriodic("z7m");
      addComponentWithDerivatives("w7m");
      componentIsNotPeriodic("w7m");
      addComponentWithDerivatives("A2");
      componentIsNotPeriodic("A2");
      addComponentWithDerivatives("B2");
      componentIsNotPeriodic("B2");
      addComponentWithDerivatives("A3");
      componentIsNotPeriodic("A3");
      addComponentWithDerivatives("B3");
      componentIsNotPeriodic("B3");
      addComponentWithDerivatives("tau2");
      componentIsNotPeriodic("tau2");
      addComponentWithDerivatives("tau3");
      componentIsNotPeriodic("tau3");
      addComponentWithDerivatives("theta7m");
      componentIsNotPeriodic("theta7m");
      addComponentWithDerivatives("rho");
      componentIsNotPeriodic("rho");
        
    }
    
    log<<"  Bibliography ";
    if(atoms.size()==5) log<<plumed.cite("Huang, Giese, Lee, York, J. Chem. Theory Comput. 10, 1538 (2014)");
    if(atoms.size()==6) log<<plumed.cite("Cremer and Pople, J. Am. Chem. Soc. 97, 1354 (1975)");
    if(atoms.size()==7) log<<plumed.cite("Boessenkool and Boeyens, J. Cryst. Mol. Struct., 10, 11–18 (1980)");
    
    log<<"\n";
    
    requestAtoms(atoms);
}

// calculator
void Puckering::calculate() {
    makeWhole();
    if (getNumberOfAtoms()==5) {
        calculate5m();
    } else if(getNumberOfAtoms()==6) {
        calculate6m();
    } else
        calculate7m();
}

void Puckering::calculate5m() {

  Vector d0,d1,d2,d3,d4,d5;

  d0=delta(getPosition(2),getPosition(1));
  d1=delta(getPosition(3),getPosition(2));
  d2=delta(getPosition(4),getPosition(3));
  d3=delta(getPosition(4),getPosition(3));
  d4=delta(getPosition(0),getPosition(4));
  d5=delta(getPosition(1),getPosition(0));

  Vector dd0,dd1,dd2,dd3,dd4,dd5;

  PLMD::Torsion t;

  double v1=t.compute(d0,d1,d2,dd0,dd1,dd2);
  double v3=t.compute(d3,d4,d5,dd3,dd4,dd5);

  double Zx=(v1+v3)/(2.0*std::cos(4.0*pi/5.0));
  double Zy=(v1-v3)/(2.0*std::sin(4.0*pi/5.0));
  double phase=std::atan2(Zy,Zx);
  double amplitude=std::sqrt(Zx*Zx+Zy*Zy);

  Vector dZx_dR[5];
  Vector dZy_dR[5];

  dZx_dR[0]=(dd5-dd4);
  dZx_dR[1]=(dd0-dd5);
  dZx_dR[2]=(dd1-dd0);
  dZx_dR[3]=(dd2+dd3-dd1);
  dZx_dR[4]=(dd4-dd3-dd2);

  dZy_dR[0]=(dd4-dd5);
  dZy_dR[1]=(dd0+dd5);
  dZy_dR[2]=(dd1-dd0);
  dZy_dR[3]=(dd2-dd3-dd1);
  dZy_dR[4]=(dd3-dd4-dd2);

  for(unsigned j=0; j<5; j++) {
    dZx_dR[j]*=(1.0/(2.0*std::cos(4.0*pi/5.0)));
  }
  for(unsigned j=0; j<5; j++) {
    dZy_dR[j]*=(1.0/(2.0*std::sin(4.0*pi/5.0)));
  }

  Vector dphase_dR[5];
  for(unsigned j=0; j<5; j++) {
    dphase_dR[j]=(1.0/(Zx*Zx+Zy*Zy))*(-Zy*dZx_dR[j] + Zx*dZy_dR[j]);
  }

  Vector damplitude_dR[5];
  for(unsigned j=0; j<5; j++) {
    damplitude_dR[j]=(1.0/amplitude)*(Zx*dZx_dR[j] + Zy*dZy_dR[j]);
  }

  Value* vzx=getPntrToComponent("Zx");
  vzx->set(Zx);
  setAtomsDerivatives (vzx,0, dZx_dR[0]);
  setAtomsDerivatives (vzx,1, dZx_dR[1]);
  setAtomsDerivatives (vzx,2, dZx_dR[2]);
  setAtomsDerivatives (vzx,3, dZx_dR[3]);
  setAtomsDerivatives (vzx,4, dZx_dR[4]);
  setBoxDerivativesNoPbc(vzx);

  Value* vzy=getPntrToComponent("Zy");
  vzy->set(Zy);
  setAtomsDerivatives (vzy,0, dZy_dR[0]);
  setAtomsDerivatives (vzy,1, dZy_dR[1]);
  setAtomsDerivatives (vzy,2, dZy_dR[2]);
  setAtomsDerivatives (vzy,3, dZy_dR[3]);
  setAtomsDerivatives (vzy,4, dZy_dR[4]);
  setBoxDerivativesNoPbc(vzy);


  Value* vph=getPntrToComponent("phs");
  vph->set(phase);
  setAtomsDerivatives (vph,0, dphase_dR[0]);
  setAtomsDerivatives (vph,1, dphase_dR[1]);
  setAtomsDerivatives (vph,2, dphase_dR[2]);
  setAtomsDerivatives (vph,3, dphase_dR[3]);
  setAtomsDerivatives (vph,4, dphase_dR[4]);
  setBoxDerivativesNoPbc(vph);

  Value* vam=getPntrToComponent("amp");
  vam->set(amplitude);
  setAtomsDerivatives (vam,0, damplitude_dR[0]);
  setAtomsDerivatives (vam,1, damplitude_dR[1]);
  setAtomsDerivatives (vam,2, damplitude_dR[2]);
  setAtomsDerivatives (vam,3, damplitude_dR[3]);
  setAtomsDerivatives (vam,4, damplitude_dR[4]);
  setBoxDerivativesNoPbc(vam);


}

void Puckering::calculate6m() {

  std::vector<Vector> r(6);
  for(unsigned i=0; i<6; i++) {
    r[i]=getPosition(i);
  }

  std::vector<Vector> R(6);
  Vector center;
  for(unsigned j=0; j<6; j++) {
    center+=r[j]/6.0;
  }
  for(unsigned j=0; j<6; j++) {
    R[j]=(r[j]-center);
  }

  Vector Rp,Rpp;
  for(unsigned j=0; j<6; j++) {
    Rp +=R[j]*std::sin(2.0/6.0*pi*j);
  }
  for(unsigned j=0; j<6; j++) {
    Rpp+=R[j]*std::cos(2.0/6.0*pi*j);
  }

  Vector n=crossProduct(Rp,Rpp);
  Vector nhat=n/modulo(n);

  Tensor dn_dRp=dcrossDv1(Rp,Rpp);
  Tensor dn_dRpp=dcrossDv2(Rp,Rpp);

  Tensor dnhat_dn= 1.0/modulo(n)*( Tensor::identity() - extProduct(nhat,nhat));
  Tensor dnhat_dRp=matmul(dnhat_dn,dn_dRp);
  Tensor dnhat_dRpp=matmul(dnhat_dn,dn_dRpp);

  std::vector<double> z(6);
  for(unsigned j=0; j<6; j++) {
    z[j]=dotProduct(R[j],nhat);
  }

  std::vector<std::vector<Vector> > dz_dR(6);
  for(unsigned j=0; j<6; j++) {
    dz_dR[j].resize(6);
  }

  for(unsigned i=0; i<6; i++)
    for(unsigned j=0; j<6; j++) {
      if(i==j) {
        dz_dR[i][j]+=nhat;
      }
      dz_dR[i][j]+=matmul(R[i],dnhat_dRp)*std::sin(2.0/6.0*pi*j);
      dz_dR[i][j]+=matmul(R[i],dnhat_dRpp)*std::cos(2.0/6.0*pi*j);
    }

  double B=0.0;
  for(unsigned j=0; j<6; j++) {
    B+=z[j]*std::cos(4.0/6.0*pi*j);
  }

  std::vector<Vector> dB_dR(6);
  for(unsigned i=0; i<6; i++)
    for(unsigned j=0; j<6; j++) {
      dB_dR[i]+=dz_dR[j][i]*std::cos(4.0/6.0*pi*j);
    }
  Vector Bsum;
  for(unsigned j=0; j<6; j++) {
    Bsum+=dB_dR[j];
  }
  for(unsigned j=0; j<6; j++) {
    dB_dR[j]-=Bsum/6.0;
  };

  double A=0.0;
  for(unsigned j=0; j<6; j++) {
    A+=z[j]*std::sin(4.0/6.0*pi*j);
  }

  std::vector<Vector> dA_dR(6);
  for(unsigned i=0; i<6; i++)
    for(unsigned j=0; j<6; j++) {
      dA_dR[i]+=dz_dR[j][i]*std::sin(4.0/6.0*pi*j);
    }
  Vector Asum;
  for(unsigned j=0; j<6; j++) {
    Asum+=dA_dR[j];
  }
  for(unsigned j=0; j<6; j++) {
    dA_dR[j]-=Asum/6.0;
  };

  double C=0.0;
  for(unsigned j=0; j<6; j++) {
    C+=z[j]*Tools::fastpow(-1.0,(j));
  }

  std::vector<Vector> dC_dR(6);
  for(unsigned i=0; i<6; i++)
    for(unsigned j=0; j<6; j++) {
      dC_dR[i]+=dz_dR[j][i]*Tools::fastpow(-1.0,(j));
    }

  Vector Csum;
  for(unsigned j=0; j<6; j++) {
    Csum+=dC_dR[j];
  }
  for(unsigned j=0; j<6; j++) {
    dC_dR[j]-=Csum/6.0;
  };


// qx
  double qx = B/std::sqrt(3);

// qx derivaties
  std::vector<Vector> dqx_dR(6);
  for(unsigned j=0; j<6; j++) {
    dqx_dR[j]=-1/std::sqrt(3) * dA_dR[j];
  }

  Value* vqx=getPntrToComponent("qx");
  vqx->set(qx);
  setAtomsDerivatives (vqx,0, dqx_dR[0] );
  setAtomsDerivatives (vqx,1, dqx_dR[1] );
  setAtomsDerivatives (vqx,2, dqx_dR[2] );
  setAtomsDerivatives (vqx,3, dqx_dR[3] );
  setAtomsDerivatives (vqx,4, dqx_dR[4] );
  setAtomsDerivatives (vqx,5, dqx_dR[5] );
  setBoxDerivativesNoPbc(vqx);

// qy
  double qy = -A/std::sqrt(3);

// qy derivatives
  std::vector<Vector> dqy_dR(6);
  for(unsigned j=0; j<6; j++) {
    dqy_dR[j]=1/std::sqrt(3) * dB_dR[j];
  }

  Value* vqy=getPntrToComponent("qy");
  vqy->set(qy);
  setAtomsDerivatives (vqy,0, dqy_dR[0] );
  setAtomsDerivatives (vqy,1, dqy_dR[1] );
  setAtomsDerivatives (vqy,2, dqy_dR[2] );
  setAtomsDerivatives (vqy,3, dqy_dR[3] );
  setAtomsDerivatives (vqy,4, dqy_dR[4] );
  setAtomsDerivatives (vqy,5, dqy_dR[5] );
  setBoxDerivativesNoPbc(vqy);

// qz
  double qz = C/std::sqrt(6);

// qz derivatives
  std::vector<Vector> dqz_dR(6);
  for(unsigned j=0; j<6; j++) {
    dqz_dR[j]=1/std::sqrt(6) * dC_dR[j];
  }

  Value* vqz=getPntrToComponent("qz");
  vqz->set(qz);
  setAtomsDerivatives (vqz,0, dqz_dR[0] );
  setAtomsDerivatives (vqz,1, dqz_dR[1] );
  setAtomsDerivatives (vqz,2, dqz_dR[2] );
  setAtomsDerivatives (vqz,3, dqz_dR[3] );
  setAtomsDerivatives (vqz,4, dqz_dR[4] );
  setAtomsDerivatives (vqz,5, dqz_dR[5] );
  setBoxDerivativesNoPbc(vqz);


// PHASE
  double phi=std::atan2(-A,B);

// PHASE DERIVATIVES
  std::vector<Vector> dphi_dR(6);
  for(unsigned j=0; j<6; j++) {
    dphi_dR[j]=1.0/(A*A+B*B) * (-B*dA_dR[j] + A*dB_dR[j]);
  }

  Value* vphi=getPntrToComponent("phi");
  vphi->set(phi);
  setAtomsDerivatives (vphi,0, dphi_dR[0] );
  setAtomsDerivatives (vphi,1, dphi_dR[1] );
  setAtomsDerivatives (vphi,2, dphi_dR[2] );
  setAtomsDerivatives (vphi,3, dphi_dR[3] );
  setAtomsDerivatives (vphi,4, dphi_dR[4] );
  setAtomsDerivatives (vphi,5, dphi_dR[5] );
  setBoxDerivativesNoPbc(vphi);

//  AMPLITUDE
  double amplitude=std::sqrt((2*(A*A+B*B)+C*C)/6);

//  AMPLITUDE DERIVATIES
  std::vector<Vector> damplitude_dR(6);
  for (unsigned j=0; j<6; j++) {
    damplitude_dR[j]=0.5*std::sqrt(2.0/6.0)/(std::sqrt(A*A+B*B+0.5*C*C)) * (2*A*dA_dR[j] + 2*B*dB_dR[j] + C*dC_dR[j]);
  }

  Value* vamplitude=getPntrToComponent("amplitude");
  vamplitude->set(amplitude);
  setAtomsDerivatives (vamplitude,0, damplitude_dR[0] );
  setAtomsDerivatives (vamplitude,1, damplitude_dR[1] );
  setAtomsDerivatives (vamplitude,2, damplitude_dR[2] );
  setAtomsDerivatives (vamplitude,3, damplitude_dR[3] );
  setAtomsDerivatives (vamplitude,4, damplitude_dR[4] );
  setAtomsDerivatives (vamplitude,5, damplitude_dR[5] );
  setBoxDerivativesNoPbc(vamplitude);

//  THETA
  double theta=std::acos( C / std::sqrt(2.*(A*A+B*B) +C*C ) );

//  THETA DERIVATIVES
  std::vector<Vector> dtheta_dR(6);
  for(unsigned j=0; j<6; j++) {
    dtheta_dR[j]=1.0/(3.0*std::sqrt(2)*amplitude*amplitude) * (C/(std::sqrt(A*A+B*B)) * (A*dA_dR[j] + B*dB_dR[j]) - std::sqrt(A*A+B*B)*dC_dR[j]);
  }
  Value* vtheta=getPntrToComponent("theta");
  vtheta->set(theta);
  setAtomsDerivatives (vtheta,0, dtheta_dR[0] );
  setAtomsDerivatives (vtheta,1, dtheta_dR[1] );
  setAtomsDerivatives (vtheta,2, dtheta_dR[2] );
  setAtomsDerivatives (vtheta,3, dtheta_dR[3] );
  setAtomsDerivatives (vtheta,4, dtheta_dR[4] );
  setAtomsDerivatives (vtheta,5, dtheta_dR[5] );
  setBoxDerivativesNoPbc(vtheta);
}

void Puckering::calculate7m() {
    // Collect positions
    std::vector<Vector> r(7);
    for (unsigned i = 0; i < 7; ++i) r[i] = getPosition(i);
    
    // Center coordinates
    std::vector<Vector> R(7);
    Vector center(0.0, 0.0, 0.0);
    for (unsigned j = 0; j < 7; ++j) center += r[j] / 7.0;
    for (unsigned j = 0; j < 7; ++j) R[j] = r[j] - center;
    
    // Compute projection vectors
    Vector Rp(0.0, 0.0, 0.0);
    Vector Rpp(0.0, 0.0, 0.0);
    for (unsigned j = 0; j < 7; ++j) {
        Rp  += R[j] * std::sin(2.0/7.0 * pi * j);
        Rpp += R[j] * std::cos(2.0/7.0 * pi * j);
    }
    
    // Normal to ring plane
    Vector n    = crossProduct(Rp, Rpp);
    Vector nhat = n / modulo(n);
    
    // Derivatives of normal
    Tensor dn_dRp  = dcrossDv1(Rp, Rpp);
    Tensor dn_dRpp = dcrossDv2(Rp, Rpp);
    Tensor dnhat_dn    = 1.0/modulo(n) * (Tensor::identity() - extProduct(nhat, nhat));
    Tensor dnhat_dRp   = matmul(dnhat_dn, dn_dRp);
    Tensor dnhat_dRpp  = matmul(dnhat_dn, dn_dRpp);
    
    // Distances from plane
    std::vector<double> z(7);
    for (unsigned j = 0; j < 7; ++j) z[j] = dotProduct(R[j], nhat);
    
    // Derivative of z
    std::vector<std::vector<Vector>> dz_dR(7, std::vector<Vector>(7));
    for (unsigned i = 0; i < 7; ++i) {
        for (unsigned j = 0; j < 7; ++j) {
            if (i == j) dz_dR[i][j] += nhat;
            dz_dR[i][j] += matmul(R[i], dnhat_dRp) * std::sin(2.0/7.0 * pi * j);
            dz_dR[i][j] += matmul(R[i], dnhat_dRpp) * std::cos(2.0/7.0 * pi * j);
        }
    }
    
    // Cremer-Pople A2 and B2
    double A2 = 0.0, B2 = 0.0;
    for (unsigned j = 0; j < 7; ++j) {
        A2 += z[j] * std::cos(4.0/7.0 * pi * j);
        B2 += z[j] * std::sin(4.0/7.0 * pi * j);
    }
    Value* vA2 = getPntrToComponent("A2"); vA2->set(A2);
    Value* vB2 = getPntrToComponent("B2"); vB2->set(B2);
    
    // Derivatives A2, B2
    std::vector<Vector> dA2_dR(7), dB2_dR(7);
    for (unsigned i = 0; i < 7; ++i) {
        for (unsigned j = 0; j < 7; ++j) {
            dA2_dR[i] += dz_dR[j][i] * std::cos(4.0/7.0 * pi * j);
            dB2_dR[i] += dz_dR[j][i] * std::sin(4.0/7.0 * pi * j);
        }
    }
    Vector A2sum(0.0, 0.0, 0.0), B2sum(0.0, 0.0, 0.0);
    for (unsigned j = 0; j < 7; ++j) {
        A2sum += dA2_dR[j];
        B2sum += dB2_dR[j];
    }
    for (unsigned j = 0; j < 7; ++j) {
        dA2_dR[j] -= A2sum/7.0;
        dB2_dR[j] -= B2sum/7.0;
        setAtomsDerivatives(vA2, j, dA2_dR[j]);
        setAtomsDerivatives(vB2, j, dB2_dR[j]);
    }
    setBoxDerivativesNoPbc(vA2);
    setBoxDerivativesNoPbc(vB2);
    
    // Cremer-Pople A3 and B3
    double A3 = 0.0, B3 = 0.0;
    for (unsigned j = 0; j < 7; ++j) {
        A3 += z[j] * std::cos(6.0/7.0 * pi * j);
        B3 += z[j] * std::sin(6.0/7.0 * pi * j);
    }
    Value* vA3 = getPntrToComponent("A3"); vA3->set(A3);
    Value* vB3 = getPntrToComponent("B3"); vB3->set(B3);
    
    // Derivatives A3, B3
    std::vector<Vector> dA3_dR(7), dB3_dR(7);
    for (unsigned i = 0; i < 7; ++i) {
        for (unsigned j = 0; j < 7; ++j) {
            dA3_dR[i] += dz_dR[j][i] * std::cos(6.0/7.0 * pi * j);
            dB3_dR[i] += dz_dR[j][i] * std::sin(6.0/7.0 * pi * j);
        }
    }
    Vector A3sum(0.0, 0.0, 0.0), B3sum(0.0, 0.0, 0.0);
    for (unsigned j = 0; j < 7; ++j) {
        A3sum += dA3_dR[j];
        B3sum += dB3_dR[j];
    }
    for (unsigned j = 0; j < 7; ++j) {
        dA3_dR[j] -= A3sum/7.0;
        dB3_dR[j] -= B3sum/7.0;
        setAtomsDerivatives(vA3, j, dA3_dR[j]);
        setAtomsDerivatives(vB3, j, dB3_dR[j]);
    }
    setBoxDerivativesNoPbc(vA3);
    setBoxDerivativesNoPbc(vB3);
    
    //  AMPLITUDE 2
    double q2 = std::sqrt((2.0/7.0) * (A2*A2 + B2*B2));
    
    //  AMPLITUDE 2 DERIVATIVES
    std::vector<Vector> dq2_dR(7);
    for (unsigned j = 0; j < 7; ++j) {
        dq2_dR[j] = std::sqrt(2.0/7.0) * (1.0/std::sqrt(A2*A2 + B2*B2)) * (A2*dA2_dR[j] + B2*dB2_dR[j]);
    }
    
    Value* vq2 = getPntrToComponent("q2");
    vq2->set(q2);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vq2, j, dq2_dR[j]);
    setBoxDerivativesNoPbc(vq2);
    
    //  AMPLITUDE 3
    double q3 = std::sqrt((2.0/7.0) * (A3*A3 + B3*B3));
    
    //  AMPLITUDE 3 DERIVATIVES
    std::vector<Vector> dq3_dR(7);
    for (unsigned j = 0; j < 7; ++j) {
        dq3_dR[j] = std::sqrt(2.0/7.0) * (1.0/std::sqrt(A3*A3 + B3*B3)) * (A3*dA3_dR[j] + B3*dB3_dR[j]);
    }
    
    Value* vq3 = getPntrToComponent("q3");
    vq3->set(q3);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vq3, j, dq3_dR[j]);
    setBoxDerivativesNoPbc(vq3);
    
    // PHASE 2
    double phi2 = std::atan2(-B2, A2);
    
    // PHASE 2 DERIVATIVES
    std::vector<Vector> dphi2_dR(7);
    for (unsigned j = 0; j < 7; ++j) {
        dphi2_dR[j] = (1.0/(A2*A2 + B2*B2)) * (B2*dA2_dR[j] - A2*dB2_dR[j]);
    }
    
    Value* vphi2 = getPntrToComponent("phi2");
    vphi2->set(phi2);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vphi2, j, dphi2_dR[j]);
    setBoxDerivativesNoPbc(vphi2);
    
    // PHASE 3
    double phi3 = std::atan2(-B3, A3);
    
    // PHASE 3 DERIVATIVES
    std::vector<Vector> dphi3_dR(7);
    for (unsigned j = 0; j < 7; ++j) {
        dphi3_dR[j] = (1.0/(A3*A3 + B3*B3)) * (B3*dA3_dR[j] - A3*dB3_dR[j]);
    }
    
    Value* vphi3 = getPntrToComponent("phi3");
    vphi3->set(phi3);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vphi3, j, dphi3_dR[j]);
    setBoxDerivativesNoPbc(vphi3);
    
    // GENERAL AMPLITUDE rho
    double rho = std::sqrt((2.0/7.0) * (A2*A2 + B2*B2 + A3*A3 + B3*B3));
    
    
    // rho DERIVATIVES
    std::vector<Vector> drho_dR(7);
    for (unsigned j = 0; j < 7; ++j) {
        drho_dR[j] = (2.0/7.0)/rho * ( A2 * dA2_dR[j] + B2 * dB2_dR[j] + A3 * dA3_dR[j] + B3 * dB3_dR[j] );
    }
    
    Value* vrho = getPntrToComponent("rho");
    vrho->set(rho);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vrho, j, drho_dR[j]);
    setBoxDerivativesNoPbc(vrho);
    
    
    // ANGULAR ROTATION or BOATNESS MEASURE by Epsilon
    double theta7m = std::atan2(q3, q2);
    
    // Epsilon DERIVATIVES
    std::vector<Vector> dtheta7m_dR(7);
    for (unsigned j = 0; j < 7; ++j) {
        dtheta7m_dR[j] = (1.0/(std::sqrt(A2*A2 + B2*B2) * std::sqrt(A3*A3 + B3*B3) * (A2*A2 + B2*B2 + A3*A3 + B3*B3))) *
        (A3 * dA3_dR[j] * (A2*A2 + B2*B2) - A2 * dA2_dR[j] * (A3*A3 + B3*B3) + B3 * dB3_dR[j] * (A2*A2 + B2*B2) - B2 * dB2_dR[j] * (A3*A3 + B3*B3));
    }
    Value* vtheta7m = getPntrToComponent("theta7m"); vtheta7m->set(theta7m);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vtheta7m, j, dtheta7m_dR[j]);
    setBoxDerivativesNoPbc(vtheta7m);
    
    // Nonlinear Amplitude tau2
    double tau2_squared_base = (2.0/7.0) * (A2*A2 + B2*B2);
    double tau2 = std::sqrt(tau2_squared_base) + std::sqrt(tau2_squared_base + 1.0);
    
    // Derivatives
    std::vector<Vector> dtau2_dR(7);
    double dtau2_factor = 0.5 / std::sqrt(tau2_squared_base) + 0.5 / std::sqrt(tau2_squared_base + 1.0);
    for (unsigned j = 0; j < 7; ++j) {
        Vector dtau2_dx = (4.0/7.0) * (A2 * dA2_dR[j] + B2 * dB2_dR[j]);
        dtau2_dR[j] = dtau2_factor * dtau2_dx;
    }
    Value* vtau2 = getPntrToComponent("tau2"); vtau2->set(tau2);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vtau2, j, dtau2_dR[j]);
    setBoxDerivativesNoPbc(vtau2);
    
    // Nonlinear Amplitude tau3
    double tau3_squared_base = (2.0/7.0) * (A3*A3 + B3*B3);
    double tau3 = std::sqrt(tau3_squared_base) + std::sqrt(tau3_squared_base + 1.0);
    
    // Derivatives
    std::vector<Vector> dtau3_dR(7);
    double dtau3_factor = 0.5 / std::sqrt(tau3_squared_base) + 0.5 / std::sqrt(tau3_squared_base + 1.0);
    for (unsigned j = 0; j < 7; ++j) {
        Vector dtau3_dx = (4.0/7.0) * (A3 * dA3_dR[j] + B3 * dB3_dR[j]);
        dtau3_dR[j] = dtau3_factor * dtau3_dx;
    }
    Value* vtau3 = getPntrToComponent("tau3"); vtau3->set(tau3);
    for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(vtau3, j, dtau3_dR[j]);
    setBoxDerivativesNoPbc(vtau3);
    
    // qx, qy, qz, qw
    double x7m = A2 * std::sqrt(2.0/7.0);
    double y7m = -B2 * std::sqrt(2.0/7.0);
    double z7m = -B3 * std::sqrt(2.0/7.0);
    double w7m = A3 * std::sqrt(2.0/7.0);
    
    std::vector<std::string> labels7m = {"x7m", "y7m", "z7m", "w7m"};
    std::vector<double> qvals7m = {x7m, y7m, z7m, w7m};
    std::vector<std::vector<Vector>> dqvals_dR7m(4, std::vector<Vector>(7));
    
    for (unsigned j = 0; j < 7; ++j) {
        dqvals_dR7m[0][j] = dA2_dR[j] * std::sqrt(2.0/7.0);
        dqvals_dR7m[1][j] = -dB2_dR[j] * std::sqrt(2.0/7.0);
        dqvals_dR7m[2][j] = -dB3_dR[j] * std::sqrt(2.0/7.0);
        dqvals_dR7m[3][j] = dA3_dR[j] * std::sqrt(2.0/7.0);
    }
    
    for (int idx = 0; idx < 4; ++idx) {
        Value* v = getPntrToComponent(labels7m[idx]);
        v->set(qvals7m[idx]);
        for (unsigned j = 0; j < 7; ++j) setAtomsDerivatives(v, j, dqvals_dR7m[idx][j]);
        setBoxDerivativesNoPbc(v);
    }
}
}
}
