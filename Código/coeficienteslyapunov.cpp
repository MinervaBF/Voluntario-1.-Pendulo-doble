#include <cmath>
#include <math.h>
#include <fstream>

#define g 9.8

double thetapuntopunto1(double theta1, double theta2, double thetapunto1, double thetapunto2);
double thetapuntopunto2(double theta1, double theta2, double thetapunto1, double thetapunto2);
double thetapunto(double theta1, double theta2, double hamiltoniano);

int main()
{
    double theta11, theta21, vtheta11, vtheta21, theta12, theta22, vtheta12, vtheta22;
    double exponente, dist_0, distancia, expmedio, perturbacion;
    double k1[2], k2[2], k3[2], k4[2], l1[2], l2[2], l3[2], l4[2], h, t, E1;
    double c1[2], c2[2], c3[2], c4[2], m1[2], m2[2], m3[2], m4[2], E2;

    int pasos;
    FILE *f1, *f2, *f3, *f4, *f5;


    //Abrimos los archivos
    f1=fopen("exponentelyapunov.txt","w");
    f2=fopen("pendulo1theta1.txt","w");
    f3=fopen("pendulo2theta1.txt","w");
    f4=fopen("pendulo1theta2.txt","w");
    f5=fopen("pendulo2theta2.txt","w");


    //Para calcular los coeficientes de Lyapunov emplearemos la fórmula:
    //lambda_i=1/h*ln(d(t_i)/d_0)

    //Tomaremos dos péndulos dobles con condiciones iniciales próximas e iremos calculando las distancias de uno a otro

    //Escribimos el tiempo durante el que se representará la animación:
    pasos=0.25*3600; //Pasamos horas a segundos.

    //Introducimos las condiciones iniciales del primer péndulo doble:
    theta11=2.*M_PI/180;
	theta21=0.*M_PI/180;
	E1=5.;
	vtheta21=0;
	vtheta11=thetapunto(theta11,theta21,E1);

	perturbacion=10e-8;

	//Introducimos las condiciones iniciales del segundo péndulo doble:
	theta12=theta11+perturbacion;
	theta22=theta12+perturbacion;
	E2=5.;
	vtheta22=0;
	vtheta12=thetapunto(theta12,theta22,E2);

	h=0.01;

	//Calculamos la distancia inicial
	dist_0=sqrt((theta11-theta12)*(theta11-theta12)+(theta21-theta22)*(theta21-theta22)+h*h*(vtheta11-theta12)+h*h*(vtheta21-theta22));

    fprintf(f1, "%e\t%e\n", 0., 0.);
    fprintf(f2, "%e\t%e\n", 0., theta11);
    fprintf(f3, "%e\t%e\n", 0., theta12);
    fprintf(f4, "%e\t%e\n", 0., theta21);
    fprintf(f5, "%e\t%e\n", 0., theta22);

    expmedio=0.;

    //Comenzamos el ciclo:

    for(t=0.01; t<=pasos; t=t+h)
    {
        //APLICAMOS EL MÉTODO DE RUNGE-KUTTA 4.
        //Calculamos k_1, l1, c1, m1:
        k1[0]=h*(vtheta11);
        l1[0]=h*thetapuntopunto1(theta11, theta21, vtheta11, vtheta21);
        k1[1]=h*vtheta21;
        l1[1]=h*thetapuntopunto2(theta11, theta21, vtheta11, vtheta21);

        c1[0]=h*(vtheta12);
        m1[0]=h*thetapuntopunto1(theta12, theta22, vtheta12, vtheta22);
        c1[1]=h*vtheta22;
        m1[1]=h*thetapuntopunto2(theta12, theta22, vtheta12, vtheta22);

        //Calculamos k_2, l2, c2, m2:
        k2[0]=h*(vtheta11+l1[0]/2.);
        l2[0]=h*thetapuntopunto1(theta11+k1[0]/2., theta21+k1[1]/2., vtheta11+l1[0]/2., vtheta21+l1[1]/2.);
        k2[1]=h*(vtheta21+l1[1]/2.);
        l2[1]=h*thetapuntopunto2(theta11+k1[0]/2., theta21+k1[1]/2., vtheta11+l1[0]/2., vtheta21+l1[1]/2.);

        c2[0]=h*(vtheta12+m1[0]/2.);
        m2[0]=h*thetapuntopunto1(theta12+c1[0]/2., theta22+c1[1]/2., vtheta12+m1[0]/2., vtheta22+m1[1]/2.);
        c2[1]=h*(vtheta22+m1[1]/2.);
        m2[1]=h*thetapuntopunto2(theta12+c1[0]/2., theta22+c1[1]/2., vtheta12+m1[0]/2., vtheta22+m1[1]/2.);

        //Calculamos k_3, l3, c3, m3:
        k3[0]=h*(vtheta11+l2[0]/2.);
        l3[0]=h*thetapuntopunto1(theta11+k2[0]/2., theta21+k2[1]/2., vtheta11+l2[0]/2., vtheta21+l2[1]/2.);
        k3[1]=h*(vtheta21+l2[1]/2.);
        l3[1]=h*thetapuntopunto2(theta11+k2[0]/2., theta21+k2[1]/2., vtheta11+l2[0]/2., vtheta21+l2[1]/2.);

        c3[0]=h*(vtheta12+m2[0]/2.);
        m3[0]=h*thetapuntopunto1(theta12+c2[0]/2., theta22+c2[1]/2., vtheta12+m2[0]/2., vtheta22+m2[1]/2.);
        c3[1]=h*(vtheta22+m2[1]/2.);
        m3[1]=h*thetapuntopunto2(theta12+c2[0]/2., theta22+c2[1]/2., vtheta12+m2[0]/2., vtheta22+m2[1]/2.);

        //Calculamos k_4, l4, c4, m4:
        k4[0]=h*(vtheta11+l3[0]);
        l4[0]=h*thetapuntopunto1(theta11+k3[0], theta21+k3[1], vtheta11+l3[0], vtheta21+l3[1]);
        k4[1]=h*(vtheta21+l3[1]);
        l4[1]=h*thetapuntopunto2(theta11+k3[0], theta21+k3[1], vtheta11+l3[0], vtheta21+l3[1]);

        c4[0]=h*(vtheta12+m3[0]);
        m4[0]=h*thetapuntopunto1(theta12+c3[0], theta22+c3[1], vtheta12+m3[0], vtheta22+m3[1]);
        c4[1]=h*(vtheta22+m3[1]);
        m4[1]=h*thetapuntopunto2(theta12+c3[0], theta22+c3[1], vtheta12+m3[0], vtheta22+m3[1]);

        //Calculamos los nuevos valores de theta1, theta2, vtheta1 y vtheta2:
        theta11=theta11+1./6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        theta21=theta21+1./6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        vtheta11=vtheta11+1./6*(l1[0]+2*l2[0]+2*l3[0]+l4[0]);
        vtheta21=vtheta21+1./6*(l1[1]+2*l2[1]+2*l3[1]+l4[1]);

        theta12=theta12+1./6*(c1[0]+2*c2[0]+2*c3[0]+c4[0]);
        theta22=theta22+1./6*(c1[1]+2*c2[1]+2*c3[1]+c4[1]);
        vtheta12=vtheta12+1./6*(m1[0]+2*m2[0]+2*m3[0]+m4[0]);
        vtheta22=vtheta22+1./6*(m1[1]+2*m2[1]+2*m3[1]+m4[1]);

        //Calculamos la distancia actualizada:
        distancia=sqrt((theta11-theta12)*(theta11-theta12)+(theta21-theta22)*(theta21-theta22)+h*h*(vtheta11-theta12)+h*h*(vtheta21-theta22));

        exponente=1./h*log(distancia/dist_0);

        //Escribimos en un archivo la evolución del exponente de Lyapunov
        fprintf(f1, "%e\t%e\n", t, exponente);

        //Escribimos las posiciones de los pendulos en archivos para compararlas después.
        if (t>pasos-20){
            fprintf(f2, "%e\t%e\n", t, theta11);
            fprintf(f3, "%e\t%e\n", t, theta12);
            fprintf(f4, "%e\t%e\n", t, theta21);
            fprintf(f5, "%e\t%e\n", t, theta22);
        }
        expmedio=expmedio+exponente;
    }

    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    //Calculamos el exponente de Lyapunov medio:
    expmedio=expmedio/pasos;

    return 0;
}

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

s=-sqrt(E-2*g*(1-cos(theta1))-g*(1-cos(theta2)));
return s;

}
