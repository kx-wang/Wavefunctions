#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <eigen-3.3.7/Eigen/Core> // need this for MatrixXd to be recognized
#include <eigen-3.3.7/Eigen/Eigen>

using namespace std;
using namespace Eigen;

/* constants, global variables */
double hbar = 1.054571628e-34; // Planck's constant [J.s]
double m0 = 9.109389e-31;      // electron effective mass [kg]
double eV = 1.602177e-19;      // electron volt [J]

double GetMin(vector <double> aVector)
{
    /* returns the minimum value in a vector */
    double minV = 200;
    for (int i=0; i<aVector.size(); i++){
        if (minV > aVector[i])
        {
            minV = aVector[i];
        }
    }
    return minV;
}

double GetMax(vector <double> aVector)
{
    /* returns the maximum value in a vector */
    double maxV = 0;
    for (int i=0; i<aVector.size(); i++){
        if (maxV < aVector[i])
        {
            maxV = aVector[i];
        }
    }
    return maxV;
}

MatrixXd ConstructHamiltonian(vector <double> m_z, vector <double> Vcrystal, double del_z, int Npts)
{
    /* returns the Hamiltonian that will be solved */
    MatrixXd Hamiltonian(Npts-2,Npts-2);
    vector <double> mL, mR, mM;
    vector <double> ai, bi, ci;
    mL.assign(Npts-2, 0);
    mR.assign(Npts-2, 0);
    mM.assign(Npts-2, 0);
    ai.assign(Npts-2, 0);
    bi.assign(Npts-2, 0);
    ci.assign(Npts-2, 0);

    for (int i = 0; i <= Npts-3; i++) {
        mL[i] = m_z[i] + m_z[i+1];
        mR[i] = m_z[i+1] + m_z[i+2];
        mM[i] = 1/((1/mL[i]) + (1/mR[i]));
    }

    for (int i = 0; i <= Npts - 3; i++) {
        ai[i] = -pow(hbar,2) / (mL[i] * pow(del_z, 2));
        ci[i] = -pow(hbar, 2) / (mR[i] * pow(del_z, 2));
        bi[i] = (pow(hbar, 2) / (mM[i] * pow(del_z, 2))) + Vcrystal[i + 1];
    }

    for (int i = 0; i < Npts - 2; i++) // assigns zero to all elements of matrix Hamiltonian
    {
        for (int j = 0; j < Npts - 2; j++)
        {
            Hamiltonian(i, j) = 0;
        }
    }
    for (int i = 0; i <= Npts - 3; i++)
    {
        Hamiltonian(i, i) = bi[i]; // place on diagonal
    }
    for (int i = 1; i <= Npts - 3; i++)
    {
        Hamiltonian(i, i - 1) = ai[i]; // place on lower diagonal
    }

    for (int i = 0; i <= Npts - 4; i++)
    {
        Hamiltonian(i, i + 1) = ci[i];  // place on upper diagonal
    }

    return Hamiltonian;
}



int main()
{
    /*  define the design parameters */
    double Al_frac = 0.2;
    double elec_field = 0;

    // assign the barrier and well widths
    vector <double> bw;
    vector <double> bw_sum;
    vector <int> isBarrier;
    int bw_length = 3;
    double tot_length = 0.0;

    bw.assign(bw_length, 0); // assign 3 zeros to this vector
    isBarrier.assign(bw_length, 0);

    bw[0] = 80; isBarrier[0] = 1;
    bw[1] = 160; isBarrier[1] = 0;
    bw[2] = 80; isBarrier[2] = 1;

    for (int i = 0; i < bw_length; i++) {
        bw[i] = bw[i] * 1e-10;
        tot_length = tot_length + bw[i];
    }

    // define z coordinate
    double del_z = 2e-10;
    vector <double> z;
    int Npts = (int) floor((tot_length - del_z) / del_z)+1;
    z.assign(Npts, 0);
    z[0] = del_z;

    for (int i = 1; i < Npts ; i++) {
        z[i] = z[i - 1] + del_z;
    }

    // define the Al fraction profile, effective mass profile, and the potential profile
    vector <double> m_z;
    vector <double> Vcrystal;
    vector <double> Al_frac_z;

    double left = 0;
    double right = bw[0];

    Vcrystal.assign(Npts, 0);
    m_z.assign(Npts, 0);


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
    }
    for (int i = 0; i < Npts; i++) {
        if ((z[i] > left)& (z[i] < right)) {
            Al_frac_z[i] = Al_frac;
        }
    }

    // using Al_frac_z, calculate the potential profile and effective mass as function of z
    for (int i = 0; i < Npts; i++) {
        m_z[i] = (0.067 + 0.083 * Al_frac_z[i]) * m0; // effective mass
        Vcrystal[i] = 0.65 * (1.36 + 0.22 * Al_frac_z[i]) * Al_frac_z[i] * eV; // conduction band offset
        Vcrystal[i] = Vcrystal[i] - eV * elec_field * (z[i] - z[0]);           // apply the electric field

    }

    // Make it so that the potential profile starts at zero */
    double zeroV = GetMin(Vcrystal);
    for (int i = 0; i < Npts; i++) {
        Vcrystal[i] = Vcrystal[i] - zeroV;
    }
    Vcrystal[Npts - 1] = Vcrystal[0]; // for easier convergence, make the last point equal to the first


    
    /* construct the Hamiltonian which will be solved for */
    MatrixXd Hamiltonian = ConstructHamiltonian(m_z, Vcrystal, del_z, Npts);

    /* Call on the eigen solver*/
    cout << "calling eigen solver" << endl;
    EigenSolver<MatrixXd> eig;
    eig.compute(Hamiltonian);
    VectorXd eigen_vals = eig.eigenvalues().real();
    VectorXd eigen_energies(eigen_vals.size());



    /* determine which eigenvalues are the eigen energies of interest, print results */
    double maxV = GetMax(Vcrystal); // get maximum of Vcrystal
    double minV = 0; // minimum of Vcrystal is zero, since we have zeroed it above

    cout << "these are the eigen energies in [eV]" << endl;
    for (int i = 0; i < eigen_vals.size(); i++)
    {
        if ((eigen_vals(i) < maxV) & (eigen_vals(i) > minV)) // only accept eigen values that are within the potential well
        {
            eigen_energies(i) = eigen_vals(i);
            cout << eigen_energies(i)/eV << ", ";
        }
    }
    cout << endl;
    cout << "The maximum of the potential profile is: " << maxV/eV << "eV" << endl;


}
