#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define xmin 0.0
#define ymin 0.0
#define xmax 2.0
#define ymax 2.0
#define nx 100
#define ny 100
#define tm 2.0
#define pi M_PI
#define alpha 0.9  //under relaxation
#define tol pow(10,-8)

double analytical(double x, double y)
{
    double t;
    t=tm*sin(pi*x/xmax)*sinh(pi*y/xmax) / sinh(pi*ymax/xmax);
    return t;
}

int main(){
    double dx=(xmax-xmin)/(nx-1.0);
    double dy=(ymax-ymin)/(ny-1.0);
    double p=dx/dy;
    int n=nx*ny;
    double x[n], y[n], b[n];
    double ap[n], ae[n], aw[n], an[n], as[n];
    
    FILE *fp; 
    fp=fopen("Numerical.txt","w");
    FILE *fp1;
    fp1=fopen("Analytical.txt","w");
    FILE *fp2;
    fp2=fopen("Mid_vertical_numerical.txt","w");
    FILE *fp3;
    fp3=fopen("Mid_vertical_analytical.txt","w");
    FILE *fp4;
    fp4=fopen("Mid_horizontal_numerical.txt","w");
    FILE *fp5;
    fp5=fopen("Mid_horizontal_analytical.txt","w");

        
    //initializing all vectors at index 1
    as[0]=0.0;
    an[0]=p*p;
    ap[0]=-2.0*(1.0+p*p);
    ae[0]=1.0;
    aw[0]=0.0;
    b[0]=0.0;
    for(int i=1;i<ny;i++){
        aw[i]=0.0;
        ae[i]=1.0;
        as[i]=p*p;
        an[i]=p*p;
        ap[i]=-2.0*(1.0+p*p);
        b[i]=0.0;
    }

    //initializing all vectors till index n-ny
    for(int i=ny;i<=(n-ny-1);i++){
        aw[i]=1.0;
        ae[i]=1.0;
        as[i]=p*p;
        an[i]=p*p;
        ap[i]=-2.0*(1.0+p*p);
        b[i]=0.0;
    }

    //initializing all vectors till index n-1
    for(int i=n-ny;i<(n-1);i++){
        aw[i]=1.0;
        ae[i]=0.0;
        as[i]=p*p;
        an[i]=p*p;
        ap[i]=-2.0*(1.0+p*p);
        b[i]=0.0;
    }

    //initializing all vectors at index n
    aw[n-1]=1.0;
    ae[n-1]=0.0;
    as[n-1]=p*p;
    an[n-1]=0.0;
    ap[n-1]=-2.0*(1.0+p*p);
    b[n-1]=0.0;

    //updating vectors at boundary points
    for(int i=0;i<nx;i++)
    {
        double xpos=xmin+i*dx;
        for(int j=0;j<ny;j++)
        {
            int l=i*ny+j; 
           if(i==0 || i==(nx-1) || j==0 || j==(ny-1))
           {
               ap[l]=1.0;
               ae[l]=aw[l]=an[l]=as[l]=0.0;

               //Applying boundary condition at upper boundary
               if(j==(ny-1)){
               b[l]=tm*sin(pi*xpos/xmax);
               }
           }
        }
    }

    double lw[n], ls[n], lp[n], un[n], ue[n], mnw[n], mse[n];
    double temp, err, sum, RHS[n];
    double np[n], nn[n], ns[n], ne[n], nw[n];

    //initializing all the vector till index ny
    lw[0]=aw[0];    ls[0]=as[0];    lp[0]=ap[0];
    un[0]=an[0]/lp[0];  ue[0]=ae[0]/lp[0];

    for(int i=1;i<ny;i++){
        lw[i]=aw[i];    ls[i]=as[i] / (1+alpha*ue[i-1]);
        lp[i]=ap[i] + alpha*ls[i]*ue[i-1] - ls[i]*un[i-1];
        un[i]=an[i]/lp[i];  ue[i]=(ae[i] - alpha*ls[i]*ue[i-1]) / lp[i];
    }

    //initializing all the vectors till index n
    for(int i=ny;i<n;i++){  
        lw[i]=aw[i] / (1+alpha*un[i-ny]);   ls[i]=as[i] / (1+alpha*ue[i-1]);
        lp[i]=ap[i] + alpha*(lw[i]*un[i-ny] + ls[i]*ue[i-1]) - lw[i]*ue[i-ny] - ls[i]*un[i-1];
        un[i]=(an[i] - alpha*lw[i]*un[i-ny]) / lp[i];   ue[i]=(ae[i] - alpha*ls[i]*ue[i-1]) / lp[i];
    }

    //initializing mnw and mse till index ny
    mnw[0]=mse[0]=0.0;
    for(int i=1;i<ny;i++){
        mnw[i]=0.0;
        mse[i]=ls[i]*ue[i-1];
    }

    //initializing mnw and mse till index n
    for(int i=ny;i<n;i++){
        mnw[i]=lw[i]*un[i-ny];
        mse[i]=ls[i]*ue[i-1];
    }

    //Calculating remaining 5 diagonals of N vector
    for(int i=0;i<n;i++){
        np[i]=alpha*(mnw[i]+mse[i]);    nn[i]=-alpha*mnw[i];
        ns[i]=-alpha*mse[i];    ne[i]=-alpha*mse[i];
        nw[i]=-alpha*mnw[i];
    }

    //taking inital guess as 0 at all the grid points
    for(int i=0;i<n;i++){
        x[i]=0.0;
    }
    
    int it=0;
    do
    {
        err=0.0;    sum=0; 
        double Nx[n];

        //initializing Nx till index ny-1
        Nx[0]=np[0]*x[0] + nn[0]*x[1] + ne[0]*x[ny] + mse[0]*x[ny-1];
        for(int i=1;i<ny-1;i++){
            Nx[0]=ns[i]*x[i-1] + np[i]*x[i] + nn[i]*x[i+1] + ne[i]*x[i+ny] + mse[i]*x[i+ny-1];
        }

        //initializing Nx at index ny
        Nx[ny-1]=ns[ny-1]*x[ny-2] + np[ny-1]*x[ny-1] + nn[ny-1]*x[ny] + ne[ny-1]*x[ny-1+ny] + mse[ny-1]*x[ny-1+ny-1] + mnw[ny-1]*x[0];

        //initializing Nx till index n-ny
        for(int i=ny;i<n-ny;i++){
            Nx[i]=nw[i]*x[i-ny] + ns[i]*x[i-1] + np[i]*x[i] + nn[i]*x[i+1] + ne[i]*x[i+ny] + mse[i]*x[i+ny-1] + mnw[i]*x[i-ny+1];
        }        

        //initializing Nx at index n-ny+1
        Nx[n-ny]=nw[n-ny]*x[n-ny-ny] + ns[n-ny]*x[n-ny-1] + np[n-ny]*x[n-ny] + nn[n-ny]*x[n-ny+1] + mse[n-ny]*x[n-1] + mnw[n-ny]*x[n-ny-ny+1];

        //initializing Nx till index n-1
        for(int i=n-ny;i<n-1;i++){
            Nx[i]=nw[i]*x[i-ny] + ns[i]*x[i-1] + np[i]*x[i] + nn[i]*x[i+1] + mnw[i]*x[i-ny+1];
        }

        //initializing Nx at index n
        Nx[n-1]=nw[n-1]*x[n-1-ny] + ns[n-1]*x[n-2] + np[n-1]*x[n-1] + mnw[n-1]*x[n-1-ny+1];

        //Initializing RHS vector
        for(int i=0;i<n;i++){
            RHS[i]=b[i]+Nx[i];
        }

        //FORWARD SUBSTITUTION
        //finding y (=Ux) till index ny
        y[0]=RHS[0]/lp[0];
        for(int i=1;i<ny;i++){
            y[i]=(RHS[i] - y[i-1]*ls[i])/lp[i];
        }

        //finding y till index n
        for(int i=ny;i<n;i++){
            y[i]=(RHS[i] - y[i-1]*ls[i] - y[i-ny]*lw[i])/lp[i];
        }

        //BACKWARD SUBSTITUTION
        //finding x till index n-1-ny
        temp=y[n-1];
        err+=fabs(temp-x[n-1]);
        x[n-1]=temp;
        for(int i=n-2;i>(n-1-ny);i--){
            temp=y[i] - x[i+1]*un[i];
            err+=fabs(temp-x[i]);
            x[i]=temp;
        }

        //finding x till index 1
        for(int i=n-1-ny;i>=0;i--){
            temp=y[i] - x[i+1]*un[i] - x[i+ny]*ue[i];
            err+=fabs(temp-x[i]);
            x[i]=temp;
        }
        it++;
        printf("%d: %.11f\n",it,err);
    } while (err>tol);

    //printing result in output file
        fprintf(fp,"x\ty\tT\n");
        fprintf(fp1,"x\ty\tT\n");
        for(int i=0;i<nx;i++){
            double xpos=xmin+i*dx;
            for(int j=0;j<ny;j++)
            {
                int l=i*ny+j;
                double ypos=ymin+j*dy;

                //Print Numerical solution
                fprintf(fp,"%f\t%f\t%0.8f\n",xpos,ypos,x[l]);

                //Print Analytical solution
                double t_ana=analytical(xpos,ypos);
                fprintf(fp1,"%f\t%f\t%0.8f\n",xpos,ypos,t_ana);
            }
        }

        int i1, i2, imid, j1, j2, jmid;
        double xpos, ypos;

        //Print temperature along mid-verical line
        fprintf(fp2,"T\ty\n");
        fprintf(fp3,"T\ty\n");
        if(nx%2==0){
            i1=nx/2;
            i2=i1-1;
            xpos = ((xmin+i1*dx)+(xmin+i2*dx))*0.5;
            for(int j=0;j<ny;j++)
            {
                int l1=i1*ny+j;  int l2=i2*ny+j;
                ypos=ymin+j*dy;
                double t = (x[l1]+x[l2])*0.5;
                double t_ana=analytical(xpos,ypos);
                //Print Numerical solution along mid-vertical line
                fprintf(fp2,"%f\t%0.8f\n",t,ypos);
                //Print Analytical solution along mid-vertical line
                fprintf(fp3,"%f\t%0.8f\n",t_ana,ypos);
            }
        }
        else{
            imid=(nx-1)/2;
            xpos = xmin+imid*dx;
            for(int j=0;j<ny;j++){
                ypos=ymin+j*dy;
                int l=imid*ny+j;
                double t_ana=analytical(xpos,ypos);
                //Print Numerical solution along mid-vertical line
                fprintf(fp2,"%f\t%0.8f\n",x[l],ypos);
                //Print Analytical solution along mid-vertical line
                fprintf(fp3,"%f\t%0.8f\n",t_ana,ypos);
            }
        }

        //Print temperture along mid-horizontal line
        fprintf(fp4,"x\tT\n");
        fprintf(fp5,"x\tT\n");
        if(ny%2==0){
            j1=ny/2;
            j2=j1-1;
            ypos = ((ymin+j1*dy)+(ymin+j2*dy))*0.5;
            for(int i=0;i<nx;i++)
            {
                int l1=i*ny+j1;  int l2=i*ny+j2;
                xpos=xmin+i*dx;
                double t = (x[l1]+x[l2])*0.5;
                double t_ana=analytical(xpos,ypos);
                //Print Numerical solution along mid-horizontal line
                fprintf(fp4,"%f\t%0.8f\n",xpos,t);
                //Print Analytical solution along mid-horizontal line
                fprintf(fp5,"%f\t%0.8f\n",xpos,t_ana);
            }
        }
        else{
            jmid=(ny-1)/2;
            ypos = ymin+jmid*dy;
            for(int i=0;i<nx;i++){
                xpos=xmin+i*dx;
                int l=i*ny+jmid;
                double t_ana=analytical(xpos,ypos);
                //Print Numerical solution along mid-horizontal line
                fprintf(fp4,"%f\t%0.8f\n",xpos,x[l]);
                //Print Analytical solution along mid-horizontal line
                fprintf(fp5,"%f\t%0.8f\n",xpos,t_ana);
            }
        }
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);

    return 0;
}