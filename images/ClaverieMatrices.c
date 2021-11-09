This code assumes the radial grid defined in a vector 'grid'
of nGridSize elements, counting from 1 to nGridsize,


double Bmat(int k, int j){
if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
if ((k==1) && (j ==1)) return -1*grid[1]*grid[1]/4 + grid[1]*grid[2]/6 + grid[2]*grid[2]/12;
else if ((k ==1) && (j == 2)) return (-1*grid[1]*grid[1] + grid[2]*grid[2])/12;
else if ((k ==2) && (j == 1)) return (-1*grid[1]*grid[1] + grid[2]*grid[2])/12;
else if ((k ==2) && (j == 2)) return (-1*grid[1] + grid[3])*(grid[1] + 2*grid[2] + grid[3])/12;
else if ((k == nGridSize-1) && (j == nGridSize-1)) return (grid[nGridSize] - grid[nGridSize-2])*(grid[nGridSize] + 2*grid[nGridSize-1] + grid[nGridSize-2])/12;
else if ((k==nGridSize-1)&&(j == nGridSize)) return (grid[nGridSize]*grid[nGridSize] - grid[nGridSize-1]*grid[nGridSize-1])/12;
else if ((k==nGridSize) && (j == nGridSize-1)) return (grid[nGridSize]*grid[nGridSize] - grid[nGridSize-1]*grid[nGridSize-1])/12;
else if ((k==nGridSize) && (j == nGridSize)) return grid[nGridSize]*grid[nGridSize]/4 - grid[nGridSize]*grid[nGridSize-1]/6 - grid[nGridSize-1]*grid[nGridSize-1]/12;
else if (j == k-1) return (grid[k]*grid[k] - grid[k-1]*grid[k-1])/12;
else if (j == k) return (-1*grid[k-1] + grid[k+1])*(2*grid[k] + grid[k-1] + grid[k+1])/12;
else if (j == k+1) return (-1*grid[k]*grid[k] + grid[k+1]*grid[k+1])/12;
else return 0;
}

double A1(int k, int j){
if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
if ((k==1) && (j ==1)) return (grid[1] + grid[2])/(2*(-1*grid[1] + grid[2]));
else if ((k ==1) && (j == 2)) return (grid[1] + grid[2])/(2*(grid[1] - grid[2]));
else if ((k==nGridSize) && (j == nGridSize-1)) return (grid[nGridSize] + grid[nGridSize-1])/(2*(-1*grid[nGridSize] + grid[nGridSize-1]));
else if ((k==nGridSize) && (j == nGridSize)) return (grid[nGridSize] + grid[nGridSize-1])/(2*(grid[nGridSize] - grid[nGridSize-1]));
else if (j == k-1) return (grid[k] + grid[k-1])/(2*(-1*grid[k] + grid[k-1]));
else if (j == k) return grid[k]*(-1*grid[k-1] + grid[k+1])/((grid[k] - grid[k-1])*(-1*grid[k] + grid[k+1]));
else if (j == k+1) return (grid[k] + grid[k+1])/(2*(grid[k] - grid[k+1]));
else return 0;
};

double A2(int k, int j){
if ((k > nGridSize) || (j > nGridSize) || (k < 1) || (j < 1)) return 0;
if ((k==1) && (j ==1)) return -1*grid[1]*grid[1]/4 - grid[1]*grid[2]/6 - grid[2]*grid[2]/12;
else if ((k ==1) && (j == 2)) return -1*grid[1]*grid[1]/12 - grid[1]*grid[2]/6 - grid[2]*grid[2]/4;
else if ((k ==2) && (j == 1)) return (3*grid[1]*grid[1] + 2*grid[1]*grid[2] + grid[2]*grid[2])/12;
else if ((k ==2) && (j == 2)) return (grid[1] - grid[3])*(grid[1] + 2*grid[2] + grid[3])/12;
else if ((k == nGridSize-1) && (j == nGridSize-1)) return (-1*grid[nGridSize] + grid[nGridSize-2])*(grid[nGridSize] + 2*grid[nGridSize-1] + grid[nGridSize-2])/12;
else if ((k==nGridSize-1)&&(j == nGridSize)) return -1*grid[nGridSize]*grid[nGridSize]/4 - grid[nGridSize]*grid[nGridSize-1]/6 - grid[nGridSize-1]*grid[nGridSize-1]/12;
else if ((k==nGridSize) && (j == nGridSize-1)) return (grid[nGridSize]*grid[nGridSize] + 2*grid[nGridSize]*grid[nGridSize-1] + 3*grid[nGridSize-1]*grid[nGridSize-1])/12;
else if ((k==nGridSize) && (j == nGridSize)) return (3*grid[nGridSize]*grid[nGridSize] + 2*grid[nGridSize]*grid[nGridSize-1] + grid[nGridSize-1]*grid[nGridSize-1])/12;
else if (j == k-1) return (grid[k]*grid[k] + 2*grid[k]*grid[k-1] + 3*grid[k-1]*grid[k-1])/12;
else if (j == k) return (grid[k-1] - grid[k+1])*(2*grid[k] + grid[k-1] + grid[k+1])/12;
else if (j == k+1) return -1*grid[k]*grid[k]/12 - grid[k]*grid[k+1]/6 - grid[k+1]*grid[k+1]/4;
else return 0;
};
