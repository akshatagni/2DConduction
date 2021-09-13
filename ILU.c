#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define xmin 0.0
#define ymin 0.0
#define xmax 1.0
#define ymax 1.0
#define nx 120
#define ny 120
#define pi M_PI
#define tol pow(10,-8)

double analytical(double x, double y)
{
    double t=0, dt;
    double sum=0.0;
    int n=1;
    do
    {
        sum+=(2.0/pi) * ( (pow(-1,n+1)+1) / n ) * sin(n*pi*x) * sinh(n*pi*(y)) / sinh(n*pi);
        dt=fabs(t-sum);
        t=sum;
        n++;
    } while (dt>tol);
    return t;
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
        for(int j=0;j<ny;j++)
        {
            int l=i*ny+j;
           if(i==0 || i==(nx-1) || j==0 || j==(ny-1))
           {
               ap[l]=1.0;
               ae[l]=aw[l]=an[l]=as[l]=0.0;

               //T=1 at upper boundary
               if(j==(ny-1)){
               b[l]=1.0;
               }
           }
        }
    }

    //lw=aw and ls=as
    double lp[n], un[n], ue[n], mnw[n], mse[n], RHS[n], temp, err, sum;


    //initializing lp and un till index ny
    lp[0]=ap[0];
    un[0]=an[0]/lp[0];
    for(int i=1;i<ny;i++){
        lp[i]=ap[i]-as[i]*un[i-1];
        un[i]=an[i]/lp[i];
    }

    //initializing ue till index ny
    for(int i=0;i<ny;i++){
        ue[i]=ae[i]/lp[i];
    }

    //initializing lp, un and ue til index n
    for(int i=ny;i<n;i++){
        lp[i]=ap[i] - aw[i]*ue[i-ny] - as[i]*un[i-1];
        un[i]=an[i]/lp[i];
        ue[i]=ae[i]/lp[i];
    }

    //initializing mnw and mse till index ny
    mnw[0]=mse[0]=0.0;
    for(int i=1;i<ny;i++){
        mnw[i]=0.0;
        mse[i]=as[i]*ue[i-1];
    }

    //initializing mnw and mse till index n
    for(int i=ny;i<n;i++){
        mnw[i]=aw[i]*un[i-ny];
        mse[i]=as[i]*ue[i-1];
    }

    //taking inital guess as 0 at all the grid points
    for(int i=0;i<n;i++){
        x[i]=0.0;
    }
    int it=0;
    do
    {
        err=0.0;    sum=0;
        //initializing RHS till index ny-1
        for(int i=0;i<(ny-1);i++){
        RHS[i]=b[i]+mse[i]*x[i+ny-1];
        }

        //initializing RHS till index n-ny+1
        for(int i=ny-1;i<=(n-ny);i++){
            RHS[i]=b[i] + mse[i]*x[i+ny-1] + mnw[i]*x[i-ny+1];
            }

        //initializing RHS till index n
        for(int i=n-ny+1;i<n;i++){
            RHS[i]=b[i] + mnw[i]*x[i-ny+1];
            }

        //FORWARD SUBSTITUTION
        //finding y (=Ux) till index ny
        y[0]=RHS[0]/lp[0];
        for(int i=1;i<ny;i++){
            y[i]=(RHS[i] - y[i-1]*as[i])/lp[i];
        }

        //finding y till index n
        for(int i=ny;i<n;i++){
            y[i]=(RHS[i] - y[i-1]*as[i] - y[i-ny]*aw[i])/lp[i];
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
