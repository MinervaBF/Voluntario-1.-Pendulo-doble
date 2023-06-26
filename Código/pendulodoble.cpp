#include <cmath>
#include <math.h>
#include <fstream>

#define g 9.8

double thetapuntopunto1(double theta1, double theta2, double thetapunto1, double thetapunto2);
double thetapuntopunto2(double theta1, double theta2, double thetapunto1, double thetapunto2);
double thetapunto(double theta1, double theta2, double hamiltoniano);

int main()
{
    double theta1, theta2, vtheta1, vtheta2; 
    double pospendulo1[2], pospendulo2[2];
    double k1[2], k2[2], k3[2], k4[2], l1[2], l2[2], l3[2], l4[2], h, t, E;
    int pasos;
    FILE *f1, *f2, *f3, *f4;


    //Abrimos los archivos
    f1=fopen("pendulodoble.txt","w");
    f2=fopen("poincaretheta1theta2.txt","w");
    f3=fopen("poincaretheta2theta2punto.txt","w");
    f4=fopen("poincaretheta2theta1punto.txt","w");

    //Escribimos el tiempo durante el que se representará la animación:
    pasos=1*3600; //Pasamos horas a segundos.

    //Introducimos las condiciones iniciales:
    theta1=45*M_PI/180;
	theta2=45*M_PI/180;
	E=1.;
	vtheta2=0;
	vtheta1=thetapunto(theta1,theta2,E);

    //Calculamos las posiciones iniciales de los péndulos:
    pospendulo1[0]=sin(theta1);
    pospendulo1[1]=-cos(theta1);
    pospendulo2[0]=sin(theta1)+sin(theta2);
    pospendulo2[1]=-cos(theta1)-cos(theta2);

    fprintf(f1, "%i%c\t%i\n", 0, 44, 0);
    fprintf(f1, "%e%c\t%e\n", pospendulo1[0], 44, pospendulo1[1]);
    fprintf(f1, "%e%c\t%e\n", pospendulo2[0], 44, pospendulo2[1]);
    fprintf(f1, "\n");

    h=0.01;

    //Comenzamos el ciclo:

    for(t=0.01; t<=pasos; t=t+h)
    {
        //Calculamos k_1 y l1:
        k1[0]=h*(vtheta1);
        l1[0]=h*thetapuntopunto1(theta1, theta2, vtheta1, vtheta2);
        k1[1]=h*vtheta2;
        l1[1]=h*thetapuntopunto2(theta1, theta2, vtheta1, vtheta2);

        //Calculamos k_2 y l2:
        k2[0]=h*(vtheta1+l1[0]/2.);
        l2[0]=h*thetapuntopunto1(theta1+k1[0]/2., theta2+k1[1]/2., vtheta1+l1[0]/2., vtheta2+l1[1]/2.);
        k2[1]=h*(vtheta2+l1[1]/2.);
        l2[1]=h*thetapuntopunto2(theta1+k1[0]/2., theta2+k1[1]/2., vtheta1+l1[0]/2., vtheta2+l1[1]/2.);

        //Calculamos k_3 y l3:
        k3[0]=h*(vtheta1+l2[0]/2.);
        l3[0]=h*thetapuntopunto1(theta1+k2[0]/2., theta2+k2[1]/2., vtheta1+l2[0]/2., vtheta2+l2[1]/2.);
        k3[1]=h*(vtheta2+l2[1]/2.);
        l3[1]=h*thetapuntopunto2(theta1+k2[0]/2., theta2+k2[1]/2., vtheta1+l2[0]/2., vtheta2+l2[1]/2.);

        //Calculamos k_4 y l4:
        k4[0]=h*(vtheta1+l3[0]);
        l4[0]=h*thetapuntopunto1(theta1+k3[0], theta2+k3[1], vtheta1+l3[0], vtheta2+l3[1]);
        k4[1]=h*(vtheta2+l3[1]);
        l4[1]=h*thetapuntopunto2(theta1+k3[0], theta2+k3[1], vtheta1+l3[0], vtheta2+l3[1]);

        //Calculamos los nuevos valores de theta1, theta2, vtheta1 y vtheta2:
        theta1=theta1+1./6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        theta2=theta2+1./6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        vtheta1=vtheta1+1./6*(l1[0]+2*l2[0]+2*l3[0]+l4[0]);
        vtheta2=vtheta2+1./6*(l1[1]+2*l2[1]+2*l3[1]+l4[1]);

        //Calculamos la posición de los péndulos.
        //Tenemos que tener en cuenta la posición de los ejes. La y apunta a arriba y el angulo es con respecto a la parte positiva de y.
        pospendulo1[0]=sin(theta1);
        pospendulo1[1]=-cos(theta1);
        pospendulo2[0]=sin(theta1)+sin(theta2);
        pospendulo2[1]=-cos(theta1)-cos(theta2);

        //Escribimos los datos en un fichero, separados por una coma.
        fprintf(f1, "%i%c\t%i\n", 0, 44, 0);
        fprintf(f1, "%e%c\t%e\n", pospendulo1[0], 44, pospendulo1[1]);
        fprintf(f1, "%e%c\t%e\n", pospendulo2[0], 44, pospendulo2[1]);
        fprintf(f1, "\n"); //Los datos tienen que estar separados en bloques.

        //POINCARÉ THETA 1 Y 2

        //Fijamos thetapunto1 y thetapunto2:
        if (vtheta1>0.8 && vtheta1<0.9){
            if (vtheta2>0.5 && vtheta2<0.6)
                fprintf(f2,"%e%c\t%e\n", theta1, 44, theta2);
        }


        //POINCARÉ THETA 2 Y THETAPUNTO2
        if (vtheta1>0.8 && vtheta1<0.9){
            if (theta1>0.5 && theta1<0.6)
                fprintf(f3,"%e%c\t%e\n", theta2, 44, vtheta2);
        }

        //POINCARÉ THETA2 Y THETAPUNTO1
        if (vtheta2>0.8 && vtheta2<0.9){
            if (theta1>0.5 && theta1<0.6)
                fprintf(f2,"%e%c\t%e\n", theta2, 44, vtheta1);
        }


    }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    return 0;
}

//FUNCIONES PARA CADA VARIABLE

double thetapuntopunto1(double theta1, double theta2, double thetapunto1, double thetapunto2)
{
	double s;

	s=(-3*g*sin(theta1)-g*sin(theta1-2*theta2)-2*sin(theta1-theta2)*(thetapunto2*thetapunto2+thetapunto1*thetapunto1*cos(theta1-theta2)))/(3-cos(2*theta1-2*theta2));


	return s;
}


double thetapuntopunto2(double theta1, double theta2, double thetapunto1, double thetapunto2)
{
	double s;

	s=(2*sin(theta1-theta2)*(2*g*cos(theta1)+2*thetapunto1*thetapunto1+thetapunto2*thetapunto2*cos(theta1-theta2)))/(3-cos(2*theta1-2*theta2));


	return s;
}

double thetapunto(double theta1, double theta2, double E)
{

double s;

s=-sqrt(E+2*g*(cos(theta1))+g*(cos(theta2)));
return s;

}
