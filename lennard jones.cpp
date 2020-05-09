//
//  main.cpp
//  lennard jones
//
//  Created by Chen Xie on 2019/10/29.
//  Copyright © 2019 Chen Xie. All rights reserved.
//

#include <iostream>
#include <random>
#include <cmath>
#include <time.h>
using namespace std;
const int N = 108;                             // number of particles
const double bd = 0.8;                        // particle density
const double length = pow(N/bd,1.0/3);        //  length of the cube
const double L2 = length / 2;                //  half of the length of the cube
const double epsilon = 1.0 , sigma = 1.0;    // reduce units
const double M = 1 ,t = 0;                   //  particle mass , time
const int NC = 3;
const double TEMP_melt = 5.0;               // melting temperature
const double TEMP_target = 1.0;             // reducing temperature
//   periodic boundary condition
double pbc(double a){
    while( a > L2 || a <= -L2){
        if( a > L2){
            a -= length;
        }
        if( a <= -L2){
            a += length;
        }
    }
    return a;
}
//  acceleration of Lennard-Jones potential
double aLJ(double r, double L){
    double a=0;
    double r2 = L * L ;
    double rm2 = 1 / r2;
    double rm6 = rm2 * rm2 *rm2;
    a = 48 * epsilon * r * rm2 * rm6 *(rm6 - 0.5)/M;
    return a;
}
int main(int argc, const char * argv[]) {
    // insert code here...
    std::default_random_engine e;
    std::normal_distribution <double>  d1(0,1);
    std::uniform_real_distribution <double> d2(-L2,L2);
    double v[N][3];  //  three-dimensional velocities
    double r[N][3];  //  three-dimensional coordinates
    double a[N][3];  //  three-dimensional acceleration
    double epsilon = 1.0 , sigma = 1.0; // reduce units
    double H = 0, K = 0, V = 0;         //  total energy, kinetic energy, potential energy
    double M = 1, t = 0;                //  particle quality, time
    double L[N][N];                     //  interparticle distance
    double sumx = 0, sumy = 0, sumz = 0;//  sum of the speeds  of three dimension
    double kb = 1;                     //  Boltzmann constant
    double constant = 0;
    int i,j,k;
    double dt = 0.00002;                 // time step size
    double dx[N][N];     // distance of particle in x
    double dy[N][N];     // distance of particle in y
    double dz[N][N];     // distance of particle in z
    double Vsq,SUMVsq;
    double CELL = length/NC;
    double f = 0;
    double CELL2 = 0.5 * CELL;
    //  intialization
    for (i = 0;i < N;i++){
        for(j = 0;j < 3;j++){
            v[i][j] = 0;
            a[i][j] = 0;
            r[i][j] = 0;
        }
    }
    for (i = 0;i < N;i++){
        for(j = 0;j < N ;j++){
            L[i][j] = 0;
            dx[i][j] = 0;
            dy[i][j] = 0;
            dz[i][j] = 0;
        }
    }
    //  velocity ,acceleration assignment
    for(i = 0; i< N;i++){
        for(j = 0; j < 3;j++){
            v[i][j] = d1(e);
        //  r[i][j] = d2(e);
        }
    }
    //  Build the unit cell
     r[0][0] = r[0][1] = r[0][2] = 0.0;
     r[1][0] = r[1][1] = CELL2; r[1][2] = 0.0;
     r[2][1] = r[2][2] = CELL2; r[2][0] = 0.0;
     r[3][0] = r[3][1] = CELL2; r[3][1] = 0.0;
    // Build the lattice from the unit cell
     int m = 0;
     int IX, IY, IZ, IR;
     for (IZ=1; IZ<=NC; IZ++) {
         for (IY=1; IY<=NC; IY++) {
             for (IX=1; IX<=NC; IX++) {
                 for (IR=1; IR<=4; IR++) {
                     r[IR+m-1][0] = r[IR-1][0]+(CELL*(IX-1));
                     r[IR+m-1][1] = r[IR-1][1]+(CELL*(IY-1));
                     r[IR+m-1][2] = r[IR-1][2]+(CELL*(IZ-1));
                 }
                 m = m + 4;
             }
         }
     }
     for (i = 0; i < N; i++) {
             r[i][0] -= L2;
             r[i][1] -= L2;
             r[i][2] -= L2;
     }
    // let the total momentum of the system be 0
    for( i = 0; i < N;i++){
        sumx += v[i][0];
        sumy += v[i][1];
        sumz += v[i][2];
    }
    for( i = 0;i < N;i++){
        v[i][0] -= sumx/N;
        v[i][1] -= sumy/N;
        v[i][2] -= sumz/N;
    }
    // calculate the periodic boundary condition of distance
    for(i = 0; i < N;i++){
        for(j = 0; j < N;j++){
            dx[i][j] = pbc(r[i][0] - r[j][0]);
            dy[i][j] = pbc(r[i][1] - r[j][1]);
            dz[i][j] = pbc(r[i][2] - r[j][2]);
            L[i][j] = sqrt(dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j] + dz[i][j] * dz[i][j]);
        }
    }
    // calculate the acceleration in three dimensions
    for(i = 0;i < N;i++){
        for(j = 0;j < N;j++){
            if( (i != j)  && (L[i][j] <= L2)){
                a[i][0] += aLJ(dx[i][j] , L[i][j]);
                a[i][1] += aLJ(dy[i][j] , L[i][j]);
                a[i][2] += aLJ(dz[i][j] , L[i][j]);
            }
        }
    }
    //  scaling the velocity
    constant = 3.0 * N * TEMP_melt;
    SUMVsq = 0.0;
    for (i = 0; i < N; i++){
        Vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
        SUMVsq += Vsq;
    }
    f = sqrt(constant/SUMVsq);
    for (i = 0;i < N;i++){
        v[i][0] *= f;
        v[i][1] *= f;
        v[i][2] *= f;
    }
    //  Melt the lattice at reduced temp = 5.0
    for(k = 0;k < 2000;k++){
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                r[i][j] += dt * v[i][j] + 0.5 * a[i][j] * dt * dt;
                r[i][j] = pbc(r[i][j]);
            }
        }
        for(i = 0; i < N;i++){
            for(j = 0; j < N;j++){
                dx[i][j] = pbc(r[i][0] - r[j][0]);
                dy[i][j] = pbc(r[i][1] - r[j][1]);
                dz[i][j] = pbc(r[i][2] - r[j][2]);
                L[i][j] = sqrt(dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j] + dz[i][j] * dz[i][j]);
            }
        }
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                v[i][j] += 0.5 * dt * a[i][j];
            }
        }
        // initalize the acceleration
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                a[i][j] = 0;
            }
        }
        for(i = 0;i < N;i++){
             for(j = 0;j < N;j++){
                 if( (i != j)  && (L[i][j] <= L2)){
                  a[i][0] += aLJ(dx[i][j] , L[i][j]);
                  a[i][1] += aLJ(dy[i][j] , L[i][j]);
                  a[i][2] += aLJ(dz[i][j] , L[i][j]);
                 }
             }
         }
       for(i = 0;i < N;i++){
               for(j = 0;j < 3;j++){
                   v[i][j] += 0.5 * dt * a[i][j];
               }
           }
    }
    // reduce the temperature
    constant = 3.0 * N * TEMP_target;
    SUMVsq = 0.0;
    for (i = 0; i < N; i++){
        Vsq = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
        SUMVsq += Vsq;
    }
    f = sqrt(constant/SUMVsq);
    for (i = 0;i < N;i++){
        v[i][0] *= f;
        v[i][1] *= f;
        v[i][2] *= f;
    }
    //  Melt the lattice at reduced temp = 5.0
    for(k = 0;k < 2000;k++){
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                r[i][j] += dt * v[i][j] + 0.5 * a[i][j] * dt * dt;
                r[i][j] = pbc(r[i][j]);
            }
        }
        for(i = 0; i < N;i++){
            for(j = 0; j < N;j++){
                dx[i][j] = pbc(r[i][0] - r[j][0]);
                dy[i][j] = pbc(r[i][1] - r[j][1]);
                dz[i][j] = pbc(r[i][2] - r[j][2]);
                L[i][j] = sqrt(dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j] + dz[i][j] * dz[i][j]);
            }
        }
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                v[i][j] += 0.5 * dt * a[i][j];
            }
        }
        // initalize the acceleration
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                a[i][j] = 0;
            }
        }
        for(i = 0;i < N;i++){
             for(j = 0;j < N;j++){
                 if( (i != j)  && (L[i][j] <= L2)){
                  a[i][0] += aLJ(dx[i][j] , L[i][j]);
                  a[i][1] += aLJ(dy[i][j] , L[i][j]);
                  a[i][2] += aLJ(dz[i][j] , L[i][j]);
                 }
             }
         }
       for(i = 0;i < N;i++){
               for(j = 0;j < 3;j++){
                   v[i][j] += 0.5 * dt * a[i][j];
               }
           }
    }
    //  original velocity verlet
    for(k = 0;k < 1000001;k++){
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                r[i][j] += dt * v[i][j] + 0.5 * a[i][j] * dt * dt;
                r[i][j] = pbc(r[i][j]);
            }
        }
        for(i = 0; i < N;i++){
            for(j = 0; j < N;j++){
                dx[i][j] = pbc(r[i][0] - r[j][0]);
                dy[i][j] = pbc(r[i][1] - r[j][1]);
                dz[i][j] = pbc(r[i][2] - r[j][2]);
                L[i][j] = sqrt(dx[i][j] * dx[i][j] + dy[i][j] * dy[i][j] + dz[i][j] * dz[i][j]);
            }
        }
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                v[i][j] += 0.5 * dt * a[i][j];
            }
        }
        // initalize the acceleration
        for(i = 0;i < N;i++){
            for(j = 0;j < 3;j++){
                a[i][j] = 0;
            }
        }
        for(i = 0;i < N;i++){
             for(j = 0;j < N;j++){
                 if( (i != j)  && (L[i][j] <= L2)){
                  a[i][0] += aLJ(dx[i][j] , L[i][j]);
                  a[i][1] += aLJ(dy[i][j] , L[i][j]);
                  a[i][2] += aLJ(dz[i][j] , L[i][j]);
                 }
             }
         }
       for(i = 0;i < N;i++){
               for(j = 0;j < 3;j++){
                   v[i][j] += 0.5 * dt * a[i][j];
               }
           }
      // calculate the kinetic energy
      for (i = 0; i < N; i++ ){
              K += 0.5 * M * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
          }
     // calculate the potential energy
          for (i = 0; i < N-1 ; i++ ){
              for (j = i+1; j < N; j++ ){
                  if (L[i][j] < L2 ){
                  V += 4 * epsilon * (pow(sigma/L[i][j],12)-pow(sigma/L[i][j],6));
                  }
              }
          }
   //  calculate the total energy
          H = V + K;
          t = k * dt;
        if( k % 1 == 0){
        printf("%f, %f, %f, %f,\n",t , H , K , V);
        }
 //
          H = 0;
          V = 0;
          K = 0;
        for(i = 0; i < N;i++){
            for(j = 0; j < N;j++){
                L[i][j] = 0;
            }
        }
      }
    return 0;
}
