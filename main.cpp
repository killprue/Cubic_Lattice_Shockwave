#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;


void acceleration(double &a1x,double &a1y,double &a1z,double &a2x,double &a2y,double &a2z,double &potential,double m,double k_spring,double r0,double x1,double x2,double y1,double y2,double z1,double z2,bool instantiate);

int main() {
    int Xn = 10;
    int Yn = 10;
    int Zn = 10;
    int nt = 300000;
    double m = 1.0;
    double k_spring = 1.0;
    double r0 = 1.0;
    double h=0.0025;
    
    double a1x = 0;
    double a2x = 0;
    double a1y = 0;
    double a2y = 0;
    double a1z = 0;
    double a2z = 0;
    
    double*** ax = new double**[Xn * Yn * Zn];
    double*** ay = new double**[Xn * Yn * Zn];
    double*** az = new double**[Xn * Yn * Zn];
    
    double**** XS = new double***[Xn * Yn * Zn * nt];
    double**** YS = new double***[Xn * Yn * Zn * nt];
    double**** ZS = new double***[Xn * Yn * Zn * nt];
    
    double* Momentum = new double[nt];
    double* Kinetic = new double[nt];
    double* Potential = new double[nt];

    
    for (int i = 0; i < Xn; i++) {
        ax[i] = new double*[Yn];
        ay[i] = new double*[Yn];
        az[i] = new double*[Yn];
        
        XS[i] = new double**[Yn];
        YS[i] = new double**[Yn];
        ZS[i] = new double**[Yn];
        for(int j = 0; j < Yn; j++){
            ax[i][j] = new double[Zn];
            ay[i][j] = new double[Zn];
            az[i][j] = new double[Zn];
            
            XS[i][j] = new double*[Zn];
            YS[i][j] = new double*[Zn];
            ZS[i][j] = new double*[Zn];
            for(int k = 0; k < Zn; k++){
                ax[i][j][k] = 0;
                ay[i][j][k] = 0;
                az[i][j][k] = 0;
                
                XS[i][j][k] = new double[nt];
                YS[i][j][k] = new double[nt];
                ZS[i][j][k] = new double[nt];
                for(int t = 0; t < nt; t++){
                    XS[i][j][k][t] = 0;
                    YS[i][j][k][t] = 0;
                    ZS[i][j][k][t] = 0;
                }
            }
        }
    }
    
    for(int i = 0; i < Xn; i++){
        for(int j = 0; j < Yn; j++){
            for(int k = 0; k < Zn; k++){
                XS[i][j][k][0] = i*r0;
                YS[i][j][k][0] = j*r0;
                ZS[i][j][k][0] = k*r0;
                
                XS[i][j][k][1] = i*r0;
                YS[i][j][k][1] = j*r0;
                ZS[i][j][k][1] = k*r0;
                
            }
        }
    }
    
    for(int it = 0; it < nt; it++){
        Potential[it] = 0;
        Kinetic[it] = 0;
        Momentum[it] = 0;
    }
    
    int i = 0;
    for(int j = 0; j < Yn; j++){
        for(int k = 0; k < Zn; k++){
            XS[i][j][k][1] = XS[i][j][k][0]+ 0*h + (0.5)* (1.0) * h*h;
            YS[i][j][k][1] = j*r0;
            ZS[i][j][k][1] = k*r0;
        }
    }
    
    double current_kinetic_energy = 0;
    double current_potential_energy = 0;
    double current_momentum = 0;
    double velocity = 0;
    for(i = 0; i < Xn; i++){
        for(int j = 0; j < Yn; j++){
            for(int k = 0; k < Zn; k++){
                velocity = (XS[i][j][k][1] - XS[i][j][k][0])/h;
                current_kinetic_energy = (0.5)*(m)*velocity*velocity;
                Kinetic[1] += current_kinetic_energy;
                current_momentum = m*velocity;
                Momentum[1] += current_momentum;
                
            }
        }
    }
    
    for(i = 0; i < Xn; i++){
        for(int j = 0; j < Yn; j++){
            for(int k = 0; k < Zn; k++){
                if(i != Xn-1){
                    current_potential_energy = 0;
                    acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy,m, k_spring, r0, XS[i][j][k][0], XS[i+1][j][k][0], YS[i][j][k][0], YS[i+1][j][k][0], ZS[i][j][k][0], ZS[i+1][j][k][0],true);
                    Potential[0] += current_potential_energy;
                    
                    current_potential_energy = 0;
                    acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy,m, k_spring, r0, XS[i][j][k][1], XS[i+1][j][k][1], YS[i][j][k][1], YS[i+1][j][k][1], ZS[i][j][k][1], ZS[i+1][j][k][1],true);
                    Potential[1] += current_potential_energy;
                    
            
                }
                if(j != Yn-1){
                    
                    current_potential_energy = 0;
                    acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy,m, k_spring, r0, XS[i][j][k][0], XS[i][j+1][k][0], YS[i][j][k][0], YS[i][j+1][k][0], ZS[i][j][k][0], ZS[i][j+1][k][0],true);
                    Potential[0] += current_potential_energy;
                    
                    current_potential_energy = 0;
                    acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy,m, k_spring, r0, XS[i][j][k][0], XS[i][j+1][k][0], YS[i][j][k][0], YS[i][j+1][k][0], ZS[i][j][k][0], ZS[i][j+1][k][0],true);
                    Potential[1] += current_potential_energy;
    
                
                }
                if(k != Zn-1){
                    
                    current_potential_energy = 0;
                    acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy, m, k_spring, r0, XS[i][j][k][0], XS[i][j][k+1][0], YS[i][j][k][0], YS[i][j][k+1][0], ZS[i][j][k][0], ZS[i][j][k+1][0],true);
                    Potential[0] += current_potential_energy;
                    
                    current_potential_energy = 0;
                    acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy, m, k_spring, r0, XS[i][j][k][1], XS[i][j][k+1][1], YS[i][j][k][1], YS[i][j][k+1][1], ZS[i][j][k][1], ZS[i][j][k+1][1],true);
                    Potential[1] += current_potential_energy;
                    

                }
                    
            }
        }
    }
    
    
    for(int it = 1; it < nt-1; it++){
        for(i = 0; i < Xn; i++){
            for(int j = 0; j < Yn; j++){
                for(int k = 0; k < Zn; k++){
                    ax[i][j][k] = 0;
                    ay[i][j][k] = 0;
                    az[i][j][k] = 0;
                }
            }
        }
    
        
        for(i = 0; i < Xn; i++){
            for(int j = 0; j < Yn; j++){
                for(int k = 0; k < Zn; k++){
                    if(i != Xn-1){
                        current_potential_energy = 0;
                        acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy,m, k_spring, r0, XS[i][j][k][it], XS[i+1][j][k][it], YS[i][j][k][it], YS[i+1][j][k][it], ZS[i][j][k][it], ZS[i+1][j][k][it],false);
                        Potential[it+1] += current_potential_energy;
                        
                        ax[i][j][k] = ax[i][j][k] + a1x;
                        ay[i][j][k] = ay[i][j][k] + a1y;
                        az[i][j][k] = az[i][j][k] + a1z;
                        
                        ax[i+1][j][k] = ax[i+1][j][k] + a2x;
                        ay[i+1][j][k] = ay[i+1][j][k] + a2y;
                        az[i+1][j][k] = az[i+1][j][k] + a2z;
                
                    }
                    if(j != Yn-1){
                        
                        current_potential_energy = 0;
                        acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy,m, k_spring, r0, XS[i][j][k][it], XS[i][j+1][k][it], YS[i][j][k][it], YS[i][j+1][k][it], ZS[i][j][k][it], ZS[i][j+1][k][it],false);
                        Potential[it+1] += current_potential_energy;
                        
                        ax[i][j][k] = ax[i][j][k] + a1x;
                        ay[i][j][k] = ay[i][j][k] + a1y;
                        az[i][j][k] = az[i][j][k] + a1z;
                        
                        ax[i][j+1][k] = ax[i][j+1][k] + a2x;
                        ay[i][j+1][k] = ay[i][j+1][k] + a2y;
                        az[i][j+1][k] = az[i][j+1][k] + a2z;
                    
                    }
                    if(k != Zn-1){
                        
                        current_potential_energy = 0;
                        acceleration(a1x, a1y, a1z, a2x, a2y, a2z,current_potential_energy, m, k_spring, r0, XS[i][j][k][it], XS[i][j][k+1][it], YS[i][j][k][it], YS[i][j][k+1][it], ZS[i][j][k][it], ZS[i][j][k+1][it],false);
                        Potential[it+1] += current_potential_energy;
                        
                        ax[i][j][k] = ax[i][j][k] + a1x;
                        ay[i][j][k] = ay[i][j][k] + a1y;
                        az[i][j][k] = az[i][j][k] + a1z;
                        
                        ax[i][j][k+1] = ax[i][j][k+1] + a2x;
                        ay[i][j][k+1] = ay[i][j][k+1] + a2y;
                        az[i][j][k+1] = az[i][j][k+1] + a2z;

                    }
                }
            }
        
        }
        
        for(i = 0; i < Xn; i++){
            for(int j = 0; j < Yn; j++){
                for(int k = 0; k < Zn; k++){
                    
                    XS[i][j][k][it+1] = 2*XS[i][j][k][it] - XS[i][j][k][it-1] + ax[i][j][k]*h*h;
                    YS[i][j][k][it+1] = 2*YS[i][j][k][it] - YS[i][j][k][it-1] + ay[i][j][k]*h*h;
                    ZS[i][j][k][it+1] = 2*ZS[i][j][k][it] - ZS[i][j][k][it-1] + az[i][j][k]*h*h;
                    
                    
                    velocity = (XS[i][j][k][it+1] - XS[i][j][k][it])/h;
                    current_kinetic_energy = (0.5)*(m)*velocity*velocity;
                    Kinetic[it+1] += current_kinetic_energy;
                    current_momentum = m*velocity;
                    Momentum[it+1] += current_momentum;
                    
                }
            }
        }
    }

    
    double*x1 = new double[nt];
    double*x2 = new double[nt];
    for(i = 0; i < nt; i++){
        x1[i] = 0;
        x2[i] = 0;
    }
    
    for(i = 0; i < nt; i++){
        x1[i] = XS[0][int(floor((Yn-1)/2))][int(floor((Zn-1)/2))][i];
        x2[i] = XS[Xn-1][int(floor((Yn-1)/2))][int(floor((Zn-1)/2))][i];
    }
    
    ofstream file;
    file.open("/Users/killianprue/machine_learning/nbodyModel.csv");
    file << "X1,X2,KE,PE,MV" << endl;
    for(i = 0; i < nt; i++){
        file << x1[i] << "," << x2[i] << "," << Kinetic[i] << "," << Potential[i] << "," << Momentum[i] << endl;
    }
    file.close();
    

    
    delete[] x1;
    delete[] x2;
    
    
    for (int i = 0; i < Xn; i++) {
        for(int j = 0; j < Yn; j++){
            for(int k = 0; k < Zn; k++){
                delete[] XS[i][j][k];
                delete[] YS[i][j][k];
                delete[] ZS[i][j][k];

            }
        }
    }
    
    for (int i = 0; i < Xn; i++) {
        for(int j = 0; j < Yn; j++){
            delete[] ax[i][j];
            delete[] ay[i][j];
            delete[] az[i][j];
            
            delete[] XS[i][j];
            delete[] YS[i][j];
            delete[] ZS[i][j];
        
        }
    }
    
    for (int i = 0; i < Xn; i++) {
        delete[] ax[i];
        delete[] ay[i];
        delete[] az[i];
        
        delete[] XS[i];
        delete[] YS[i];
        delete[] ZS[i];
        
    }
    
    delete[] ax;
    delete[] ay;
    delete[] az;
    
    delete[] XS;
    delete[] YS;
    delete[] ZS;
    
    delete[] Momentum;
    delete[] Kinetic;
    delete[] Potential;
    
}


void acceleration(double &a1x,double &a1y,double &a1z,double &a2x,double &a2y,double &a2z,double &potential,double m,double k_spring,double r0,double x1,double x2,double y1,double y2,double z1,double z2,bool instantiate){
    double magr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
    double f = -k_spring * (magr-r0);
    double amag = f/m;
    potential = (0.5)*k_spring * (magr-r0) * (magr-r0);


    if(instantiate == false){
        a1x = -((x2-x1)/magr)*amag;
        a1y = -((y2-y1)/magr)*amag;
        a1z = -((z2-z1)/magr)*amag;
        a2x = -a1x;
        a2y = -a1y;
        a2z = - a1z;
    }
}
