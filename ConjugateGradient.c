#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define xmin 0.0
#define ymin 0.0
#define xmax 1.0
#define ymax 1.0
#define nx 100
#define ny 100
#define pi M_PI
#define tol pow(10,-8)

//Calculate Analyticat
double analytical(double x, double y)
{
    double T=sin(2.0*pi*x)*sin(2.0*pi*y);
    return T;
}

//Caluclate RHS of Governing equation
double f(double x, double y)
{
    double val=-8.0*pi*pi*(sin(2.0*pi*x)*sin(2.0*pi*y));
    return val;
}

int main()
{
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

        
    /*Initializing all vectors at index 1, b vector is 
    initialized to 0 at all points due to the given boundary condition*/

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

    //updating ap, ae, aw, as, an at boundary points and b at interior points
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            int l=i*ny+j; 
           if(i==0 || i==(nx-1) || j==0 || j==(ny-1))
           {
               ap[l]=1.0;
               ae[l]=aw[l]=an[l]=as[l]=0.0;

           }
           else
           {
               double xpos=xmin+i*dx;   double ypos=ymin+j*dy;
               b[l]=f(xpos,ypos)*dx*dx;
           }
        }
    }

    //taking initial guess for Conjugate gradient method as 0 at all grid points
    for(int i=0;i<n;i++)
    {
        x[i]=0.0;
    }
    double d[n], r[n], Ax[n], alpha, beta;
    double err;
    int it=0;

    //initializing Ax till index ny
    Ax[0]=ap[0]*x[0] + an[0]*x[1] + ae[0]*x[ny];
    for (int i=1; i<ny; i++)
    {
        Ax[i]=as[i]*x[i-1] + ap[i]*x[i] + an[i]*x[i+1] + ae[i]*x[i+ny];
    }

    //initializing Ax till index n-ny
    for (int i=ny; i<(n-ny); i++)
    {
        Ax[i]=aw[i]*x[i-ny] + as[i]*x[i-1] + ap[i]*x[i] + an[i]*x[i+1] + ae[i]*x[i+ny];
    }

    //initalizing Ax till index n-1
    for(int i=n-ny; i<(n-1); i++){
        Ax[i]=aw[i]*x[i-ny] + as[i]*x[i-1] + ap[i]*x[i] + an[i]*x[i+1];
    }

    //initializing Ax at index n
    Ax[n-1]=aw[n-1]*x[n-1-ny] + as[n-1]*x[n-1-1] + ap[n-1]*x[n-1];

    for(int i=0;i<n;i++){
        d[i]=r[i]=b[i]-Ax[i];
    }

    do
    {
        double num=0, den=0, Ad[n];
        err=0;
        
        //initializing Ad till index ny
        Ad[0]=ap[0]*d[0] + an[0]*d[1] + ae[0]*d[ny];
        num+=r[0]*r[0];     den+=d[0]*Ad[0];
        for (int i=1; i<ny; i++)
        {
            Ad[i]=as[i]*d[i-1] + ap[i]*d[i] + an[i]*d[i+1] + ae[i]*d[i+ny];
            num+=r[i]*r[i];     den+=d[i]*Ad[i];
        }

        //initializing Ad till index n-ny
        for (int i=ny; i<(n-ny); i++)
        {
            Ad[i]=aw[i]*d[i-ny] + as[i]*d[i-1] + ap[i]*d[i] + an[i]*d[i+1] + ae[i]*d[i+ny];
            num+=r[i]*r[i];     den+=d[i]*Ad[i];
        }

        //initalizing Ad till index n-1
        for(int i=n-ny; i<(n-1); i++){
            Ad[i]=aw[i]*d[i-ny] + as[i]*d[i-1] + ap[i]*d[i] + an[i]*d[i+1];
            num+=r[i]*r[i];     den+=d[i]*Ad[i];
        }

        //initializing Ad at index n
        Ad[n-1]=aw[n-1]*d[n-1-ny] + as[n-1]*d[n-1-1] + ap[n-1]*d[n-1];
        num+=r[n-1]*r[n-1];     den+=d[n-1]*Ad[n-1];
        
        alpha=num/den;
        num=0; den=0;
        for (int i=0; i<n; i++)
        {
            x[i]=x[i] + alpha*d[i];
            den+=d[i]*r[i];
            r[i]=r[i] - alpha*Ad[i];
            num+=r[i]*r[i];
            err+=r[i]*r[i];
        }
        beta=num/den;
        for (int i=0; i<n; i++)
        {
            d[i] = r[i]+beta*d[i];
        }
        it++;
        err=sqrt(err);
        printf("%d: %.8f\n",it,err);
    } while (err>tol);  //checking norm of residue for convergence


    //Printing result in output file
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

    //Print temperature along mid-vertical line
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
            fprintf(fp2,"%f\t%f\n",t,ypos);
            //Print Analytical solution along mid-vertical line
            fprintf(fp3,"%f\t%f\n",t_ana,ypos);
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
            fprintf(fp2,"%f\t%f\n",x[l],ypos);
            //Print Analytical solution along mid-vertical line
            fprintf(fp3,"%f\t%f\n",t_ana,ypos);
        }
    }

    //Print temperature along mid-horizontal line
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
            fprintf(fp4,"%f\t%f\n",xpos,t);
            //Print Analytical solution along mid-horizontal line
            fprintf(fp5,"%f\t%f\n",xpos,t_ana);
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
            fprintf(fp4,"%f\t%f\n",xpos,x[l]);
            //Print Analytical solution along mid-horizontal line
            fprintf(fp5,"%f\t%f\n",xpos,t_ana);
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