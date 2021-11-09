assumes a vector tempgrid with elements 1...nGridSize, which holds the current grid points.
Because of the moving frame of reference, the matrices are only time-dependent in the 2x2 corner 
submatrices.  

double A1matrix(int k, int j){
	if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
	if ((k==1) && (j ==1)) return (tempgrid[1] + tempgrid[2])/(2*(-1*tempgrid[1] + tempgrid[2]));
   else if ((k ==1) && (j == 2)) return (tempgrid[1] + tempgrid[2])/(2*(tempgrid[1] - tempgrid[2]));
   else if ((k==nGridSize) && (j == nGridSize-1)) return (tempgrid[nGridSize] + tempgrid[nGridSize-1])/(2*(-1*tempgrid[nGridSize] + tempgrid[nGridSize-1]));
   else if ((k==nGridSize) && (j == nGridSize)) return (tempgrid[nGridSize] + tempgrid[nGridSize-1])/(2*(tempgrid[nGridSize] - tempgrid[nGridSize-1]));
   else if (j == k-1) return (tempgrid[k] + tempgrid[k-1])/(2*(-1*tempgrid[k] + tempgrid[k-1]));
   else if (j == k) return tempgrid[k]*(-1*tempgrid[k-1] + tempgrid[k+1])/((tempgrid[k] - tempgrid[k-1])*(-1*tempgrid[k] + tempgrid[k+1]));
   else if (j == k+1) return (tempgrid[k] + tempgrid[k+1])/(2*(tempgrid[k] - tempgrid[k+1]));
   else return 0;
};

double A2matrix(int k, int j){
	if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
  	if ((k==1) && (j ==1)) return -1*tempgrid[1]*tempgrid[1]/4 - tempgrid[1]*tempgrid[2]/6 - tempgrid[2]*tempgrid[2]/12;
   else if ((k ==1) && (j == 2)) return -1*tempgrid[1]*tempgrid[1]/12 - tempgrid[1]*tempgrid[2]/6 - tempgrid[2]*tempgrid[2]/4;
   else if ((k ==2) && (j == 1)) return (3*tempgrid[1]*tempgrid[1] + 2*tempgrid[1]*tempgrid[2] + tempgrid[2]*tempgrid[2])/12;
   else if ((k ==2) && (j == 2)) return (tempgrid[1] - tempgrid[3])*(tempgrid[1] + 2*tempgrid[2] + tempgrid[3])/12;
   else if ((k == nGridSize-1) && (j == nGridSize-1)) return (-1*tempgrid[nGridSize] + tempgrid[nGridSize-2])*(tempgrid[nGridSize] + 2*tempgrid[nGridSize-1] + tempgrid[nGridSize-2])/12;
   else if ((k==nGridSize-1)&&(j == nGridSize)) return -1*tempgrid[nGridSize]*tempgrid[nGridSize]/4 - tempgrid[nGridSize]*tempgrid[nGridSize-1]/6 - tempgrid[nGridSize-1]*tempgrid[nGridSize-1]/12;
   else if ((k==nGridSize) && (j == nGridSize-1)) return (tempgrid[nGridSize]*tempgrid[nGridSize] + 2*tempgrid[nGridSize]*tempgrid[nGridSize-1] + 3*tempgrid[nGridSize-1]*tempgrid[nGridSize-1])/12;
   else if ((k==nGridSize) && (j == nGridSize)) return (3*tempgrid[nGridSize]*tempgrid[nGridSize] + 2*tempgrid[nGridSize]*tempgrid[nGridSize-1] + tempgrid[nGridSize-1]*tempgrid[nGridSize-1])/12;
   else if (j == k-1) return (tempgrid[k]*tempgrid[k] + 2*tempgrid[k]*tempgrid[k-1] + 3*tempgrid[k-1]*tempgrid[k-1])/12;
   else if (j == k) return (tempgrid[k-1] - tempgrid[k+1])*(2*tempgrid[k] + tempgrid[k-1] + tempgrid[k+1])/12;
   else if (j == k+1) return -1*tempgrid[k]*tempgrid[k]/12 - tempgrid[k]*tempgrid[k+1]/6 - tempgrid[k+1]*tempgrid[k+1]/4;
   else return 0;
};

double A3matrix(int k, int j){
	if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
	if ((k==1) && (j ==1)) return tempgrid[2]*(tempgrid[1] + tempgrid[2])/12;
   else if ((k ==1) && (j == 2)) return -1*(tempgrid[2]*(tempgrid[1] + tempgrid[2]))/12;
   else if ((k ==2) && (j == 1)) return tempgrid[2]*(tempgrid[1] + 3*tempgrid[2])/12;
   else if ((k ==2) && (j == 2)) return (-1*(tempgrid[1]*tempgrid[2]) + 2*tempgrid[2]*tempgrid[3] + tempgrid[3]*tempgrid[3])/12;
   else if ((k == nGridSize-1) && (j == nGridSize-1)) return (tempgrid[nGridSize]*tempgrid[nGridSize-1] - 2*tempgrid[nGridSize-1]*tempgrid[nGridSize-2] - tempgrid[nGridSize-2]*tempgrid[nGridSize-2])/12;
   else if ((k==nGridSize-1)&&(j == nGridSize)) return -1*(tempgrid[nGridSize-1]*(tempgrid[nGridSize] + 3*tempgrid[nGridSize-1]))/12;
   else if ((k==nGridSize) && (j == nGridSize-1)) return tempgrid[nGridSize-1]*(tempgrid[nGridSize] + tempgrid[nGridSize-1])/12;
   else if ((k==nGridSize) && (j == nGridSize)) return -1*(tempgrid[nGridSize-1]*(tempgrid[nGridSize] + tempgrid[nGridSize-1]))/12;
   else if (j == k-1) return (3*tempgrid[k]*tempgrid[k] + 2*tempgrid[k]*tempgrid[k-1] + tempgrid[k-1]*tempgrid[k-1])/12;
   else if (j == k) return (-1*tempgrid[k-1] + tempgrid[k+1])*(2*tempgrid[k] + tempgrid[k-1] + tempgrid[k+1])/12;
   else if (j == k+1) return -1*tempgrid[k]*tempgrid[k]/4 - tempgrid[k]*tempgrid[k+1]/6 - tempgrid[k+1]*tempgrid[k+1]/12;
   else return 0;
};

double Bmatrix(int k,int j){
	if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
	if ((k==1) && (j ==1)) return -1*tempgrid[1]*tempgrid[1]/4 + tempgrid[1]*tempgrid[2]/6 + tempgrid[2]*tempgrid[2]/12;
   else if ((k ==1) && (j == 2)) return (-1*tempgrid[1]*tempgrid[1] + tempgrid[2]*tempgrid[2])/12;
   else if ((k ==2) && (j == 1)) return (-1*tempgrid[1]*tempgrid[1] + tempgrid[2]*tempgrid[2])/12;
   else if ((k ==2) && (j == 2)) return (-1*tempgrid[1] + tempgrid[3])*(tempgrid[1] + 2*tempgrid[2] + tempgrid[3])/12;
   else if ((k == nGridSize-1) && (j == nGridSize-1)) return (tempgrid[nGridSize] - tempgrid[nGridSize-2])*(tempgrid[nGridSize] + 2*tempgrid[nGridSize-1] + tempgrid[nGridSize-2])/12;
   else if ((k==nGridSize-1)&&(j == nGridSize)) return (tempgrid[nGridSize]*tempgrid[nGridSize] - tempgrid[nGridSize-1]*tempgrid[nGridSize-1])/12;
   else if ((k==nGridSize) && (j == nGridSize-1)) return (tempgrid[nGridSize]*tempgrid[nGridSize] - tempgrid[nGridSize-1]*tempgrid[nGridSize-1])/12;
   else if ((k==nGridSize) && (j == nGridSize)) return tempgrid[nGridSize]*tempgrid[nGridSize]/4 - tempgrid[nGridSize]*tempgrid[nGridSize-1]/6 - tempgrid[nGridSize-1]*tempgrid[nGridSize-1]/12;
   else if (j == k-1) return (tempgrid[k]*tempgrid[k] - tempgrid[k-1]*tempgrid[k-1])/12;
   else if (j == k) return (-1*tempgrid[k-1] + tempgrid[k+1])*(2*tempgrid[k] + tempgrid[k-1] + tempgrid[k+1])/12;
   else if (j == k+1) return (-1*tempgrid[k]*tempgrid[k] + tempgrid[k+1]*tempgrid[k+1])/12;
   else return 0;
};
