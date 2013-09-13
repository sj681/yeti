//-------------------------------------------------------------
//
//                   Program name goes here
//
//         (c) Sarah Jenkins and Holly Hathrell 2013
//
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or 
//  without modification, are permitted provided that the 
//  following conditions are met:
//
//   o Redistributions of source code must retain the above 
//     copyright notice, this list of conditions and the 
//     following disclaimer.
//
//   o Redistributions in binary form must reproduce the above 
//     copyright notice, this list of conditions and the 
//     following disclaimer in the documentation and/or other 
//     materials provided with the distribution.
//
//   o The names of its contributors may be used to endorse or 
//     promote products derived from this software without 
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
//  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
//  INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
//  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
//  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
//  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
//  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
//
//----------------------------------------------------------------
//
//   xx is a program to create simulated STEM images of in-
//   situ samples containing objects of different denities. 
//
//   Rutherford scattering of the electrons is assumed. An
//   electron beam moves over the sample and the high angle  
//   scattered electrons are recorded on a plane detector. 
//
//   
//-----------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <math.h>
#include<vector>
#include "mtrand.hpp"
#include "array4d.h"

using namespace std;
// This code is set to a nm = 1 pixel^2 scale
double res = 1;

int xa = 250/res; 					//xa, ya, za represent the maximum x,y,z co-ordinates of the sample size. Assuming the sample starts at (0,0,0).
int ya = 250/res;
int za = 250/res;

double r1 = 1e7; 					//sets the inside radius of the HAADF detector

int j = 25;						//sets the number of electrons going through each pixel


double xangle = M_PI/2;
double zangle = -M_PI/2;

// DENSITY - a function for choosing the density, Atomic weight and Atomic number of each pixel in the array at the moment we can only create geometric shapes.
//at the moment this has 3 circles of Gold and sets everywhere else to be water.

double density(double x, double y,double z, int xa, int ya, int za, double& A, double& Z, int zinitial){
    double p;

	if (((x*x)+((y/res)*(y/res))+((z/res)*((z/res)))) < 800/(res*res)){p = 19.2;A = 196;Z = 79;}
else if (((x*x)+((y/res)*(y/res))+((z+60/res)*((z+60/res)))) < 800/(res*res)){p = 19.2;A = 196;Z = 79;}

		//else if ((z<50/res) || (z> (za-50/res))) {p = 3.3;A = 140;Z= 11;} 
		//else {p = 1;A = 18;Z = 7.42;} 
//}
else  { p = 0.001; A  = 0.0001; Z = 0.0001;}
	  
	 // {p = 10.49;A = 107;Z= 47;}
	//else  
	//{p = 1;A = 18;Z = 7.42;} 
    return (p);
}


// LENGTH - to calculate the mean free path (l) for each material, this is calculated using the rutherford cross section of the material and the materials defect.
//the free path for each material is calculated using a random number times the mean free path - this is called the step.

void length(double Z, double E,double A,double NA,double p,double RND1, double& defect,double cs, double& l, double res){
    defect = 3.4*1e-3 * pow(Z,0.67)/E;
    cs = (Z*Z/(E*E)) * 5.21e-21 * ((4.0*M_PI )/(defect*(defect+1.0))) * pow(((E + 511.0)/(E+1024.0)),2.0);
    l = A/(NA * p * cs)*(1e7/res);
}


// ANGLE -  calculates the theta and phi angles, using random numbers and assigns and calculates the scattering potential

void angle(int& k, int& q, int& s,double RND2, double defect, double RND3, double &theta, double& phi, std::vector<double>& bin, double& probability, double& prob){
    
    double R;
    double R2;
    double theta2;
    double P;
    double p; 						//defines all the variables only used in this function
    
    phi = 2*RND3*M_PI; 				// calculates the angle phi, this is the angle in the x,y direction.
    theta2 = acos(1- ((2*defect* RND2)/(1+defect-RND2))); 				//calculates the angles - randomly
   // std::cout<<theta2*180/M_PI<<endl;
	prob = 1 - RND2;			
    int idx = int((theta2*18000)/M_PI);			//saves the theta value as a interger number in degrees, this is timesed by 100 so we can increase the precision to 0.01 degrees
    bin[idx] = bin[idx] + 1; 					//by adding one to the bin of the the theta value it was scattered by we can create an angluar probability graph.
    P= probability*prob;						// P calculates the culmative probability for each angle, this is to account for the forced high angle scattering.
    probability = P;
    k++;    									//counts the number of scatteirng events per electron
    q++;    									//counts the total number of scattering events
    s++;    									//resets the random scattering by changing the seed number
    theta2 =  theta2 + theta;
    theta = theta2;
	//std::cout<<theta*180/M_PI<<endl;
}


// DEDS - calculates the change in energy after each scattering event. This is dependant on the material properties and the initial energy 

double deds(double p, double Z, double A, double& changeinE, double J, double E){
    J = ((9.76*Z) + (58.5/pow(Z,0.19)))*10e-3;
    
    changeinE = -78500.0*(p*Z/(A*E*1e7))*log((1.166*E/J) +1.0);
    return(changeinE);
}


// POSITION - calculates the unit vectors in x,y,z so we can increase in incremental steps up to the step length seeing if it has changes density

void position(double theta, double phi, double& x, double& y, double& z,double& x3, double& y3, double& z3, double G, double& dP, double& P, double l, double res){
    double unitx = sin(theta)*cos(phi); 				// calculates the unit vecotrs
    double unity = sin(theta)*sin(phi);
    double unitz = cos(theta);
    x3 = x + unitx*G;   						//G increases until the step length is reached, at that point if the
    y3 = y + unity*G;   						//density hasnt changed, x,y,z have increased by the step length
    z3 = z + unitz*G;
    
    dP = (exp(-G/l))/(l*res*res);    						//calc the change in scattering potential over length dl
    P = P-dP;            						//minuses dP from 1, when it equals the minimum P the step is over
}


// FINAL - calulates the final position of the electron hitting the detector 

void final(double xdiff, double ydiff, double& x3, double& y3,double RTran[][3],double RMat[][3],std::vector<double> xyz,std::vector<double> xyzdash,double theta, double phi, double& x, double& y, double& z,double& x2, double& y2, double& q2, double res, double& thetasample){

  theta = theta + thetasample;
  //std::cout<<theta*180/M_PI<<endl;
  //  for (int i = 0; i<3; i++){
//	for (int j = 0; j<3; j++){
	// RTran[i][j] = RMat[j][i];
	//}}

										
	  xyz[0] = sin(theta)*cos(phi);
	  xyz[1] = sin(theta)*sin(phi);
	  xyz[2] = cos(theta);

	  //std::cout<<"initial"<<"\t"<<xyz[0]<<"\t"<<xyz[1]<<"\t"<<xyz[2]<<endl;
  for (int i = 0; i<3; i++){
	  double total2 = 0;
	  for (int j = 0; j<3; j++){

	  total2 = total2 + xyz[j]*RMat[i][j];
	  xyzdash[i] = total2;}}
  
	double unitx = xyz[0];
	double unity = xyz[1];
	double unitz = xyz[2];
	//std::cout<<"dash"<<"\t"<<xyzdash[0]<<"\t"<<xyzdash[1]<<"\t"<<xyzdash[2]<<endl;
    
    int z2 = 1e8;               								//the screen is 30cm away
    double r = ((z*res) - z2)/unitz;     				//calculates the distance from the electron to the screen
    x2 = x*res + unitx*r;       						//calculates the final position on the screen
    y2 = y*res + unity*r;
    double q = (x2*x2) + (y2*y2);
    q2 = sqrt(q); 
	//std::cout<<q2<<endl;
}
  // RESET - resets all the co-ordinates for a new electron

void reset(double zinitial2, double yinitial2, double xinitial2,double thetasample, double& theta, int& k, double& z, double& E, double& probability, double&x, double& y, double xinitial, double yinitial){
    E = 100.0;
    probability = 1.0;
    k = 0;
	z = zinitial2;
	theta = -thetasample;
	x = xinitial2;
	y = yinitial2;
}
void coord(double& x, double& y, double& z, double x3, double y3, double z3){
  x = x3;
  y = y3;
  z =z3;
}

void RMatrix(double& phisample, double& unitx, double& unity, double& unitz, double Rz[][3],double Rx[][3],double RMat[][3],std::vector<double> xyz,std::vector<double> xyzdash,double RND4, double RND5,double& theta, double& thetasample){
for(int i=0;i<3;i++){ 								//a for loop to calculate the rotation matrix by times the rotation matrixs in x and z.
   for(int j=0;j<3;j++){
   double total = 0;
    for(int k=0;k<3;k++){
     total = total + Rx[i][k] * Rz[k][j]; 
	  RMat[i][j]=total;
   }}}
	xyz[0] = 0;
	xyz[1] = 0;
	xyz[2] = 1;
	
	for (int i = 0; i<3; i++){ 						//This loop takes the incoming direction cosines and calculates the new direstion cosines using the rotation matrix
	  double total2 = 0;
	  for (int j = 0; j<3; j++){

		total2 = total2 + xyz[j]*RMat[i][j];
		xyzdash[i] = total2;
	  }}

  thetasample = acos(xyzdash[2]);					//for the direction cosines theta and phi are calculated
  phisample = 1.57079633;
  unitx = sin(thetasample)*cos(phisample);	//the unit vectors are worked out.
  unity = sin(thetasample)*sin(phisample);
  unitz = cos(thetasample);
}
void rotate(double& xinitial2, double& yinitial2, double& zinitial2, double& z, double RMat[][3],double unitx, double unity, double unitz,std::vector<double> xyz,std::vector<double> xyzdash, double& x, double& y, double phisample, double thetasample, double xinitial, double yinitial){

	  xyz[0] = xinitial;
	  xyz[1] = yinitial;
	  xyz[2] = 1.5*za;

for (int i = 0; i<3; i++){					 		//this loop takes the x,y,z positions and calculates the rotated x,y,z positions
	  double total2 = 0;
	  for (int j = 0; j<3; j++){
	
	  total2 = total2 + xyz[j]*RMat[i][j];
	  xyzdash[i] = total2;}}
	  //std::cout<<xyzdash[0]<<"\t"<<xyzdash[1]<<"\t"<<xyzdash[2]<<"\t"<<endl;
	  for (int r = 1; r<4000; r++){
	  if (xyzdash[1] >0) {
	  
	  y = int(xyzdash[1] - (unity*r));
	  x = int(xyzdash[0] - (unitx*r));
	  z = int(xyzdash[2] - (unitz*r));}
	  	  if (xyzdash[1] <0) {
	  
	  y = int(xyzdash[1] + (unity*r));
	  x = int(xyzdash[0] + (unitx*r));
	  z = int(xyzdash[2] + (unitz*r));}
	  
	  if (( y == ya/2) || ( y == - ya/2) || ( z == za/2) || ( z == -za/2) || ( x == xa/2) || ( x == -xa/2)) {
			xinitial2 = x;
			yinitial2 = y;
			zinitial2 = z;
			break;
	  }}}
	  
//----------------------------------------------------------- main code-------------------------------------------------------








int main(){
  
    ofstream myfile;
    myfile.open ("eloss.txt");
    ofstream myfile2;
    myfile2.open ("theta.txt"); 		//saves the x,y,z co-ordinates at each scattering event
    ofstream myfile3;
    myfile3.open ("screen.txt"); 		// saves the x,y,z cordinates when the electron hits the screen 10cm away from the block
    ofstream myfile7;
    myfile7.open ("scan.txt"); 		//saves the final images
    ofstream myfile1;
    myfile1.open ("silver.txt"); 		//creates a angular distribution graph.
	ofstream myfile8;
    myfile8.open ("initial.txt"); 		//creates a angular distribution graph.
	ofstream myfile9;
    myfile9.open ("rotate.txt"); 		//creates a angular distribution graph.
    ofstream myfile10;
    myfile10.open ("rotated.txt"); 		//creates a angular distribution graph.
	
    const double NA = 6.0221413e23; //Avagadros constant
    double RND1, RND2, RND3, RND4, RND5; 		//random numbers
    double defect,step,cs,l;
    double p,A,Z;
    double x,y,z;
    double x3,y3,z3, q2;
    double E, changeinE,J;
    double unitx, unity, unitz;
    double r,x2,y2,z2;
    double theta, phi;
    int s=0.0; 					//counts the number of scattering event
    int q = 0;
    int k = 0.0;
    double total;
    double probability, P , prob, dP;
	double thetasample,phisample;
	double xdash, ydash;
	double zinitial, xinitial2, yinitial2, zinitial2;
	double xdiff, ydiff,zdiff;
    
    //Array4D<double> scan(xa,ya,za,3);				//a 10 by 10 array each containing the density, atomic number and atomic weight
    
    std::vector<double> bin(36001,0.0);		 	//starts a vector with 361 points, one for each theta. each with a starting co-ordinate 0,0
    
    std::vector<std::vector<double> >scan2;   	//sets up an array to count how many electrons hit the detector from each pixel.
    scan2.resize(xa+1);
    for(int i=0; i<scan2.size(); i++) scan2[i].resize(ya+1);
	
    double Rz[3][3] = {{cos(zangle),-sin(zangle),0},{sin(zangle),cos(zangle),0},{0,0,1}};
	double Rx[3][3] = {{1,0,0},{0,cos(xangle),-sin(xangle)},{0,sin(xangle),cos(xangle)}};

	double RMat[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double RTran[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	
	std::vector<double> xyz(3,0.0);	
	std::vector<double> xyzdash(3,0.0);	
	
	RMatrix(phisample,unitx,unity,unitz, Rz,Rx,RMat, xyz,xyzdash, RND4,RND5,theta,thetasample);
	
    for (int yinitial = -ya/2; yinitial<ya/2; yinitial ++){      						//loops over y coordinates
        std::cout<<yinitial<<endl;                         						//prints out yintitial value to guage progression
        for (int xinitial = - xa/2; xinitial<xa/2; xinitial ++){   						//loops over x coordinates
 
		  rotate(xinitial2, yinitial2, zinitial2,z,RMat, unitx, unity,unitz,xyz, xyzdash, x,y,phisample,thetasample,xinitial, yinitial); 
		 //std::cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<xinitial<<"\t"<<yinitial<<endl;
		  for (int n=0; n<j; n++){    
				MTRand grnd;            										//declare a class of type MTRand called grnd
                int seed=(12345+s);    											//declare and seed RNG. Rnd# are not random, seed determines order
                grnd.seed(seed);
                reset(zinitial2, yinitial2,xinitial2,thetasample,theta,k,z,E,probability,x,y,xinitial, yinitial); 								//calls the reset function
                    double RND1 = grnd(); 										//selects random numbers
                    double RND2 = grnd();
                    double RND3 = grnd();
					double RND4 = grnd();
					double RND5 = grnd();

                    p = density(x,y,z,xa,ya,za,A,Z, zinitial);    				//outputs p,A,Z for the x,y,z positions
                    length(Z, E,A, NA, p, RND1, defect,cs, l,res); 				//calculates the mean free path,
                    angle(k,q,s,RND2, defect, RND3, theta, phi, bin, probability, prob); //randomly outputs theta and phi, calculated using p,Z,A
					//std::cout<<xinitial<<"\t"<<yinitial<<"\t"<<x<<"\t"<<y<<"\t"<<theta*180/M_PI<<endl;
					for (double G=1 ; G<100000000; G = G + 1/res){ 	 //a function increasing x,y,z in nm steps to check it doesnt change density before the end of the step length
					        //std::cout<<x3<<"\t"<<y3<<"\t"<<z3<<"\t"<<theta<<endl;   
							position(theta, phi, x, y, z, x3, y3,z3,G, dP,P, l,res);     //functions to calc. how x,y,z progress over time
							if ((z3>za/2) || ( y3>((ya/2))) || (x3>((xa/2))) ||(x3<((-xa/2)))||( z3<(-za/2)) || (y3<-(xa/2))) { 
							 //std::cout<<"bottom"<<"\t"<<x3<<"\t"<<y3<<"\t"<<z3<<endl;
							  	coord(x, y,z,x3, y3, z3);
								P= 1/res;											//if it has left and is traveling in the +ve z direction it will hit the detector
                                final(xdiff, ydiff, x3,y3,RTran,RMat,xyz,xyzdash,theta,phi,x,y,z,x2,y2,q2,res, thetasample);
								if (q2>r1){     									//if this is greater than r1 (the inside detector radius) it is saved to the array scan2 and to a file
                                    myfile3 <<x2<< "\t" <<y2<<"\t" <<z2<<"\t" <<r<< endl;
									//std::cout<<xinitial<<"\t"<<yinitial<<"\t"<<xinitial2<<"\t"<<yinitial2<<"\t"<<zinitial2<<endl;
									scan2[xinitial+xa/2][yinitial+ya/2] = scan2[xinitial+xa/2][yinitial+ya/2] + pow(probability,1/k)*10000;
									myfile <<x2<<"\t" << y2<< "\t" <<E<<endl;
								}
							break;
							}
							if (P<(RND1/res)){   											//if the scattering potential reaches the min value is RND1 because it is (e^-G/l)/l = RND1
							//std::cout<<"end"<<"\t"<<x3<<"\t"<<y3<<"\t"<<z3<<endl;
							  coord(x, y,z,x3, y3, z3);
							  myfile2 << n << "\t"<< x<< "\t" << y<<"\t" << z<<"\t"<<theta*180/M_PI<< endl;
							  s++;        										//changes the seed and random number order
							  P =1/res;       										//resets P
							  E = E + deds( p, Z, A, changeinE, J,  E); 		//calculates the change in energy
							  double RND1 = grnd(); 										//selects random numbers
							  double RND2 = grnd();
							  double RND3 = grnd();

							  p = density(x,y,z,xa,ya,za,A,Z,zinitial);    					//outputs p,A,Z for the x,y,z positions
							  length(Z, E,A, NA, p, RND1, defect,cs, l,res); 		//calculates the mean free path,
							  angle(k,q,s,RND2, defect, RND3, theta, phi, bin, probability, prob); //randomly outputs theta and phi, calculated using p,Z,A
  
                        }
             
/*
							else if (x3>xa-15){
							 //std::cout<<"pac"<<"\t"<<x3<<"\t"<<y3<<"\t"<<z3<<endl;
							  x = x - xa;}
							else if (y3>ya-15){
							   //std::cout<<"pac"<<"\t"<<x3<<"\t"<<y3<<"\t"<<z3<<endl;
							  y = y - ya;}
							else if (x3<15){
							   //std::cout<<"pac"<<"\t"<<x3<<"\t"<<y3<<"\t"<<z3<<endl;
							  x = x+xa;}
							else if (y3<15){
							   //std::cout<<"pac"<<"\t"<<x3<<"\t"<<y3<<"\t"<<z3<<endl;
							  y = y +ya;}*/
						   
                        
							else if (density(x3,y3,z3,xa,ya,za,A,Z,zinitial) != density(x,y,z,xa,ya,za,A,Z,zinitial)){ //if the denstiy changes evertyhing is recalculated for the new material
								//std::cout<<"change"<<"\t"<<theta*180/M_PI<<"\t"<<y3<<"\t"<<z3<<endl;
								coord(x, y,z,x3, y3, z3);
								G = 0;
								p = density(x,y,z,xa,ya,za,A,Z,zinitial);    			//recalculates the p,A,Z for the new co-ordinates
								length(Z,  E, A,NA, p,RND1,defect,cs,l,res);   	//calculates how far it can travel, P is not re-set so it wont travel the full step length
								
								//Random numbers are not reset so it carries on going in the same direction and the min value in the same
							}
					}
                
                //========================= end of scattering in the material ===========================
            }
        }
    }
    
    
    
    for(int i = 0; i<scan2.size(); i++){    //loops over the 4D array size (x and y coords) and saves how many hit the detector for each pixel (nm^2)
        for(int k=0; k<scan2.size(); k++){
            
            if((i<xa) & (k<ya) & (k>0) & (i>0)){    //if the electron is within the x and y coordinates
                myfile7 << i*res-xa/2 << "\t"<< k*res-ya/2 << "\t"<<scan2[i][k] << std::endl;     //writes these to a file
			}
        }
        myfile7 <<std::endl;
    }
    
    
  //loops over the 4D array size (x and y coords) and saves how many hit the detector for each pixel (nm^2)
    
    for ( double idx = 0; idx<36001; idx++){    //loops over each (decimal) theta and saves how many are in the specific bin to a file
        myfile1 << idx/100.0 << "\t" <<((bin[idx]*10000)/q)<<endl; // goes through all the bin numbers and stores how many are in each bin.q
    }
    
    return 0;
    
}   