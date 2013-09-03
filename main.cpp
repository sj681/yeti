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
//   yeti is a program to create simulated STEM images of in-
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

int xa = 250; 					//xa, ya, za represent the maximum x,y,z co-ordinates of the sample size. Assuming the sample starts at (0,0,0).
int ya = 250;
int za = 250;


// DENSITY - a function for choosing the density, Atomic weight and Atomic number of each pixel in the array at the moment we can only create geometric shapes.
//at the moment this has 3 circles of Gold and sets everywhere else to be water.

double density(double x, double y,double z, int xa, int ya, int za, double& A, double& Z){
    double p;
    
    if ((z<50) || (z> (za-50))) {p = 3.3;A = 140;Z= 11;} 	// Silicon Nitride membrane at the top and bottom of the material
    else if ((((x-(xa/2))*(x-(xa/2)))+(((y-(ya/2)))*((y-(ya/2))))+(((z-200))*((z-200)))) < 200){p = 19.2;A = 196.0;Z = 79;} 	//circle of radius 30, centered at (125,125,50)
    //else if ((((x-50)*(x-50))+(((y-50))*((y-50)))+(((z-100))*((z-100)))) < 200){p = 19.2;A = 196.0;Z = 79;} 	//circle of radius 26 centered at (50,50,50)
    //else if ((((x-(200))*(x-200))+(((y-200))*((y-900)))+(((z-900))*((z-180)))) < 200){p = 19.2;A = 196.0;Z = 79;} // circle of radius 20 centered at (200,200,50)
	else
    {p = 1;A = 18;Z = 7.42;} 	//everything else is water
   //{p = 19.2;A = 196.0;Z = 79;}
    return (p);
}


// LENGTH - to calculate the mean free path (l) for each material, this is calculated using the rutherford cross section of the material and the materials defect.
//the free path for each material is calculated using a random number times the mean free path - this is called the step.

void length(double Z, double E,double A,double NA,double p,double RND1, double& defect,double cs, double &l, double& step){
    defect =3.4*1e-3 * pow(Z,0.67)/E;
    cs = (Z*Z/(E*E)) * 5.21e-21 * ((4.0*M_PI )/(defect*(defect+1.0))) * pow(((E + 511.0)/(E+1024.0)),2.0);
    l = A/(NA * p * cs)*1e7;
    step = - l*log(RND1);
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
    
    //R = (1-((cos(theta2)+ defect*cos(theta2) -1 - defect)/(cos(theta2)-1-2*defect)));		//this is a re-arranged version of the previous theta equation
    //R2 = (1-((cos(theta2+0.01)+ defect*cos(theta2+0.01) -1 - defect)/(cos(theta2+0.01)-1-2*defect)));
    //prob = sqrt((R-R2)*(R-R2)); 				//the probability of each angle is equal to the intergral of the culmative R, theta graph for that point.
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
}


// DEDS - calculates the change in energy after each scattering event. This is dependant on the material properties and the initial energy 

double deds(double p, double Z, double A, double& changeinE, double J, double E){
    J = ((9.76*Z) + (58.5/pow(Z,0.19)))*10e-3;
    
    changeinE = -78500.0*(p*Z/(A*E*1e7))*log((1.166*E/J) +1.0);
    return(changeinE);
}


// POSITION - calculates the unit vectors in x,y,z so we can increase in incremental steps up to the step length seeing if it has changes density

void position(double unitx, double unity, double unitz, double theta, double phi, double& x, double& y, double& z,double& x3, double& y3, double& z3, double G, double& dP, double& P, double l){
    unitx = sin(theta)*cos(phi); 				// calculates the unit vecotrs
    unity = sin(theta)*sin(phi);
    unitz = cos(theta);
    x3 = x + unitx*G;   						//G increases until the step length is reached, at that point if the
    y3 = y + unity*G;   						//density hasnt changed, x,y,z have increased by the step length
    z3 = z + unitz*G;
    
    dP = exp(-G/l)/l;    						//calc the change in scattering potential over length dl
    P = P-dP;            						//minuses dP from 1, when it equals the minimum P the step is over
}


// FINAL - calulates the final position of the electron hitting the detector 

void final(double theta, double phi, double& x, double& y, double& z,double& x2, double& y2, double& q2){
    double unitx = sin(theta)*cos(phi);			// calculates the unit vecotrs
    double unity = sin(theta)*sin(phi);
    double unitz = cos(theta);
    int z2 = 1e8;               				//the screen is 10cm away
    double r = (z + z2)/unitz;     				//calculates the distance from the electron to the screen
    x2 = x + unitx*r;       					//calculates the final position on the screen
    y2 = y + unity*r;
    double q = (x2*x2) + (y2*y2);
    q2 = (pow(q,0.5));  						//calculates the radial length of the x,y value from the centre
}


// RESET - resets all the co-ordinates for a new electron

void reset(int& k, double& x, double& y, double& z,double& x3, double& y3, double& z3,int& xinitial, int& yinitial, double& E, double& theta, double& probability, double& prob, double& P, double RND4){
    x = xinitial;
    y = yinitial;
    z = 0;
    E = 100.0;
    theta = 24e-3*((RND4*2)-1);
    probability = 1.0;
    P = 0;
    k = 0;
}
void coord(double& x, double& y, double& z, double x3, double y3, double z3){
  x = x3;
  y = y3;
  z =z3;
}







//========================================================= main code ===================================================================








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
    myfile1.open ("bin.txt"); 		//creates a angular distribution graph.
    
    const double NA = 6.0221413e23; //Avagadros constant
    double RND1, RND2, RND3; 		//random numbers
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

    double r1 = 2e7; 				//sets the inside radius of the HAADF detector
    
    int j = 2000;					//sets the number of electrons going through each pixel
    
    Array4D<double> scan(xa,ya,za,3);				//a 10 by 10 array each containing the density, atomic number and atomic weight
    
    std::vector<double> bin(36001,0.0);		 	//starts a vector with 361 points, one for each theta. each with a starting co-ordinate 0,0
    
    std::vector<std::vector<double> >scan2;   	//sets up an array to count how many electrons hit the detector from each pixel.
    scan2.resize(xa+1);
    for(int i=0; i<scan2.size(); i++) scan2[i].resize(ya+1);
    
    
    for (int yinitial = 124; yinitial<125; yinitial ++){      						//loops over y coordinates
        std::cout<<yinitial<<endl;                         						//prints out yintitial value to guage progression
        for (int xinitial =124; xinitial<125; xinitial ++){   					//loops over x coordinates
            
            for (int n=0; n<j; n++){    										//for loop for the number of electrons up to j
                
                MTRand grnd;            											//declare a class of type MTRand called grnd
                int seed=(12345+s);    											//declare and seed RNG. Rnd# are not random, seed determines order
                grnd.seed(seed);
				double RND4 = grnd();
                reset( k, x, y, z, x3, y3, z3,xinitial, yinitial,  E,  theta,  probability, prob,  P, RND4); //calls the reset function

                while((E>0.1) & (x<xa) & (y<ya) & (z<za) & (x>-0.01) & (y>-0.01) & (z>-0.01)){ //while the electron is inside the box membrane and it has sufficient energy
                    
                    double RND1 = grnd(); 										//selects random numbers
                    double RND2 = grnd();
                    double RND3 = grnd();

                    p = density(x,y,z,xa,ya,za,A,Z);    						//outputs p,A,Z for the x,y,z positions
                    length(Z,  E, A, NA, p, RND1, defect, cs, l, step); 		//calculates the mean free path,
                    angle(k,q,s,RND2, defect, RND3, theta, phi, bin, probability, prob); //randomly outputs theta and phi, calculated using p,Z,A
                    for (int G=1 ; G<100000000; G++){  							//a function increasing x,y,z in nm steps to check it doesnt change density before the end of the step length
                         position( unitx,unity, unitz,theta, phi, x, y,z, x3, y3,z3,G, dP,P, l);     //functions to calc. how x,y,z progress over time
                        if (E<0.1){break;}
                        if ((x3>xa) || (y3>ya) || (z3>za) || (x3<-0.01) || (y3<-0.01) || (z3<-0.01)){           //if the particle has left the box
							coord(x, y,z,x3, y3, z3);
                            P= 1;       										//resets the P for each electron
                            if (theta<(3*M_PI/2) && theta>(M_PI/2)) { 			//if the electron is travelling in the -ve z direction it will hit the back scatter detector
                                scan2[xinitial][yinitial] = scan2[xinitial][yinitial] + pow(probability,1/k)*10000000;  //adds the number of electrons hitting a certain point on the back scatter detector weighted by the kth root of the probabilty
							
							}
                            
							else{ 												//if it has left and is traveling in the +ve z direction it will hit the detector
                                final(theta,phi, x, y, z,x2, y2, q2);   		//calc the unit x,y,z values
                                if (q2>r1){     								//if this is greater than r1 (the inside detector radius) it is saved to the array scan2 and to a file
                                    myfile3 <<x2<< "\t" <<y2<<"\t" <<z2<<"\t" <<r<< endl;
                                    scan2[xinitial][yinitial] = scan2[xinitial][yinitial] + pow(probability,1/l)*10000000;
									
									myfile <<x2<<"\t" << y2<< "\t" <<E<<endl;
								}
                            }
                            break;  											//breaks the 'if the particle has left the box loop'
                        }
                        
                        else if (density(x3,y3,z3,xa,ya,za,A,Z) != density(x,y,z,xa,ya,za,A,Z)){ //if the denstiy changes evertyhing is recalculated for the new material
                            
							coord(x, y,z,x3, y3, z3);
                            G = 0;
                            p = density(x,y,z,xa,ya,za,A,Z);    			//recalculates the p,A,Z for the new co-ordinates
                            length(Z,  E, A,NA, p,RND1,defect,cs,l,step);   //calculates how far it can travel, P is not re-set so it wont travel the full step length
                            
                            //Random numbers are not reset so it carries on going in the same direction and the min value in the same
                        }
                        
                        else if (P<RND1){   	//if the scattering potential reaches the min value is RND1 because it is (e^-G/l)/l = RND1
                            
							coord(x, y,z,x3, y3, z3);
                            myfile2 << n << "\t"<< x<< "\t" << y<<"\t" << z<<"\t"<<theta*180/M_PI<< endl;
                            s++;        									//changes the seed and random number order
                            P =1;       									//resets P
                            E = E + deds( p, Z, A, changeinE, J,  E); 		//calculates the change in energy
                            break;
                            
                        }
                    }
                }
                //========================= end of scattering in the material ===========================
            }
        }
    }
    
    
    
    for(int i = 0; i<scan2.size(); i++){    //loops over the 4D array size (x and y coords) and saves how many hit the detector for each pixel (nm^2)
        for(int k=0; k<scan2.size(); k++){
            
            if((i<xa) & (k<ya) & (k>0) & (i>0)){    //if the electron is within the x and y coordinates
                myfile7 << i << "\t"<< k << "\t"<<scan2[i][k] << std::endl;     //writes these to a file
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