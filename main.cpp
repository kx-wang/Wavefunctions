#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <eigen-3.3.7/Eigen/Core> // need this for MatrixXd to be recognized  
#include <eigen-3.3.7/Eigen/Eigen> 
//#include <eigen-3.3.7/Eigen/Dense> 

using namespace std;
using namespace Eigen;

int main()
{
    // constants
    double hbar = 1.054571628e-34; // planck's constant
    double m0 = 9.109389e-31;      // electron effective mass
    double eV = 1.602177e-19;      // electron volt

    // define the design parameters
    double Al_frac = 0.2;
    double elec_field = 0;

    // assign the barrier and well widths
    vector <double> bw;
    vector <double> bw_sum;
    vector <int> isBarrier; 
    int bw_length = 3;

    double tot_length = 0.0;

    bw.assign(3, 0); // assign zeros to this vector (initialization of the size)
    isBarrier.assign(3, 0);

    bw[0] = 80; isBarrier[0] = 1; 
    bw[1] = 160; isBarrier[1] = 0;
    bw[2] = 80; isBarrier[2] = 1; 
    //bw[2] = 38; isBarrier[2] = 1;
   // bw[3] = 83.7; isBarrier[3] = 0;
   // bw[4] = 38; isBarrier[4] = 1;
   // bw[5] = 160; isBarrier[5] = 0;
   // bw[6] = 100; isBarrier[6] = 1;

    for (int i = 0; i < bw_length; i++) {
        bw[i] = bw[i] * 1e-10;
        tot_length = tot_length + bw[i];
    }

    // z coordinate
    double del_z = 2e-10;
    vector <double> z;
    int Npts;

    Npts = (int) floor((tot_length - del_z) / del_z)+1;
    cout << Npts << endl; 

    cout << Npts << endl; 
    z.assign(Npts, 0);
    z[0] = del_z;

    for (int i = 1; i < Npts ; i++) {
        z[i] = z[i - 1] + del_z;
    }

    // define the Al fraction profile, effective mass profile, and the potential profile
    vector <double> m_z;
    vector <double> mL, mR, mM;
    vector <double> ai, bi, ci;
    vector <double> Vcrystal;
    vector <double> Al_frac_z;
    MatrixXd Hamiltonian(Npts,Npts);  

    double zeroV;
    double left = 0;
    double right = bw[0];

    Vcrystal.assign(Npts, 0);
    m_z.assign(Npts, 0);
    mL.assign(Npts-2, 0);
    mR.assign(Npts-2, 0);
    mM.assign(Npts-2, 0);
    ai.assign(Npts-2, 0);
    bi.assign(Npts-2, 0);
    ci.assign(Npts-2, 0);
    Al_frac_z.assign(Npts, 0);

    for (int j = 0; j <= bw_length-2; j++)
    {
        if (isBarrier[j] == 1)
        {
            for (int i = 0; i < Npts; i++) {
                if ( (z[i] > left)& (z[i] < right) ) {
                    Al_frac_z[i] = Al_frac;                   
                }
            }
        }
        left = left + bw[j];
        right = right + bw[j+1];
       // cout << "new left=" << bw[j] << ", new right=" << bw[j + 1] << endl;
        
    } 
    for (int i = 0; i < Npts; i++) {
        if ((z[i] > left)& (z[i] < right)) {
            Al_frac_z[i] = Al_frac;
        }
    }
 
   
    for (int i = 0; i < Npts; i++) {
        m_z[i] = (0.067 + 0.083 * Al_frac_z[i]) * m0;
        Vcrystal[i] = 0.65 * (1.36 + 0.22 * Al_frac_z[i]) * Al_frac_z[i] * eV; // conduction band offset
        Vcrystal[i] = Vcrystal[i] - eV * elec_field * (z[i] - z[0]);           // apply the electric field
        
    }
   
    // Make it so that the potential profile starts at zero
    // find the minimum
    zeroV = Vcrystal[0];
    for (int i = 0; i < Npts; i++) {
        if (zeroV > Vcrystal[i]) {
            zeroV = Vcrystal[i];
        }
    }
    // add that to the potential profile 
    for (int i = 0; i < Npts; i++) {
        Vcrystal[i] = Vcrystal[i] - zeroV;       
    }
    Vcrystal[Npts - 1] = Vcrystal[0];
    

    // construct the Hamiltonian which will be solved for 
    for (int i = 0; i <= Npts-3; i++) {
        mL[i] = m_z[i] + m_z[i+1]; // m_z(1:end - 2) + m_z(2:end - 1)
        mR[i] = m_z[i+1] + m_z[i+2]; // (m_z(2:end - 1) + m_z(3:end));
        mM[i] = 1/((1/mL[i]) + (1/mR[i])); // 1. / (1. / mL + 1. / mR);
    }
    
    ofstream myfile;
    myfile.open("myfile.txt");
    for (int i = 0; i <= Npts - 3; i++) {
        ai[i] = -pow(hbar,2) / (mL[i] * pow(del_z, 2));
        ci[i] = -pow(hbar, 2) / (mR[i] * pow(del_z, 2));
        bi[i] = (pow(hbar, 2) / (mM[i] * pow(del_z, 2))) + Vcrystal[i + 1];
        myfile << ai[i] << " " << bi[i] << " " << ci[i] << endl;
    }
    myfile.close();
        
    for (int i = 0; i < Npts - 2; i++)
    {
        for (int j = 0; j < Npts - 2; j++)
        {
            Hamiltonian(i, j) = 0;
        }
    }
    for (int i = 0; i <= Npts - 3; i++)
    {
        Hamiltonian(i, i) = bi[i]; // on the diagonal 
    } 
    for (int i = 1; i <= Npts - 3; i++)
    {
        Hamiltonian(i, i - 1) = ai[i]; // lower diagonal  
    }
    
    for (int i = 0; i <= Npts - 4; i++)
    {
        Hamiltonian(i, i + 1) = ci[i];  // upper diagonal    
    }

    ofstream myfileH;
    myfileH.open("fileH.txt");
    for (int i = 0; i < Npts - 2; i++)
    {
        for (int j = 0; j < Npts - 2; j++)
        {
            myfileH << Hamiltonian(i,j) << " " ;
        }
        myfileH << endl; 
    }
    myfileH.close(); 

    
    cout << "calling eigen solver" << endl;
    EigenSolver<MatrixX2d> eig;
    eig.compute(Hamiltonian);
    VectorXd eigen_vals = eig.eigenvalues().real();
    VectorXd eigen_energies; 
    cout << eigen_vals << endl; 

    double maxV = 0.65 * (1.36 + 0.22 * Al_frac) * Al_frac * eV;
    double minV = 0; 

    for (int i = 0; i < eigen_vals.size(); i++)
    {
        if ((eigen_vals[i] < maxV) & (eigen_vals[i] > minV))
        {
            eigen_energies[i] = eigen_vals[i];
        }
    }

    cout << eigen_energies;

}
