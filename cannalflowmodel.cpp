#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

//TODO: Revisar si se pueden pasar parametros por referencia para acelerar ejecucion del codigo
//TODO: Pensar si los calados deberian de ser variables privadas

class River {

    /*
    La clase River contiene las funciones que caracterizan el flujo de un canal de superficie abierta
    */

    public:

        //Obligatorios para inicializar la clase

            int steps;                      // Cantidad de puntos en el eje x
            int mode;                       // 1: Por encima de y_normal. 2: Entre n y_normal e y_critica. 3: Por debajo de y_critica

        //Predefinidos (Se pueden modificar)

            double discharge;                // Caudal
            double length;                   // Longitud del canal
            double bed_slope;                // Pendiente del canal. Si es positiva va hacia abajo
            double side_slope;               // Pendiente de los laterales del canal
            double bottom_width;             // Anchura del fondo
            double manning;                  // Numero de Manning
            double gravity;                  // Aceleracion de la gravedad
            double newton_rhapson_delta;     // Delta para dar por bueno el resultado en el metodo Newton-Rhapson

        //Constructores de clase

            River(
                int steps,
                int mode
                );
            
            River(
                int steps,
                int mode,
                double discharge,
                double length,
                double bed_slope,
                double side_slope,
                double bottom_width,
                double manning,
                double gravity,
                double newton_rhapson_delta
            );
        
        //Metodos de clase publicos
        void printcalados();
        int escribedatos();
        int escribeplot();

    private:

        //Atributos privados de la clase

            double runge_kutta_step;
            double delta_x;
            double *positions = (double *)malloc(steps * sizeof(double));
            double *depth = (double *)malloc(steps * sizeof(double));
            double *height = (double *)malloc(steps * sizeof(double));
            double *bed = (double *)malloc(steps * sizeof(double));
            double *plot_normal_depth = (double *)malloc(steps * sizeof(double));
            double *plot_critical_depth = (double *)malloc(steps * sizeof(double));  
            double normal_depth;
            double critical_depth;
            const char *profile_channel_type;

        //Metodos privados de la clase

            double area(double height);
            double Pw(double height);
            double T(double height);
            double fn(double yn);
            double dfn(double yn);
            double fc(double yc);
            double dfc(double yc);
            
            double newtonrhapson(bool b);
            void boundaryconditions();
            double diffeq(double y); 
            void rungekutta();
            void prepareplot();

};

//Definicion de los metodos de la clase River

    River::River(int _steps, int _mode){
        /*
        Constructor de clase para dos argumentos de entrada
        */

        steps = _steps;
        mode = _mode;
        discharge = 25.0;
        length = 4000.0;
        bed_slope = 0.001;
        side_slope = 1.0;
        bottom_width = 5.0;
        manning = 0.025;
        gravity = 9.80665;
        newton_rhapson_delta = 0.001;
        
        delta_x = length / (steps - 1);
        
        critical_depth = newtonrhapson(1);
        if (bed_slope != 0){
            normal_depth = newtonrhapson(0);
        }
        else{
            normal_depth = critical_depth + 1;
        }
        if (critical_depth < 0 || normal_depth < 0){
            cout << "El algoritmo Newton-Rhapson no converge" << endl;
            exit(1);
        }
        boundaryconditions();
        rungekutta();
        prepareplot();
    }

    River::River(int _steps, int _mode, double _discharge, double _length, double _bed_slope, double _side_slope, double _bottom_width, double _manning, double _gravity, double _newton_rhapson_delta){
        /*
        Constructor de clase para todos argumentos de entrada
        */

        steps = _steps;
        mode = _mode;
        discharge = _discharge;
        length = _length;
        bed_slope = _bed_slope;
        side_slope = _side_slope;
        bottom_width = _bottom_width;
        manning = _manning;
        gravity = _gravity;
        newton_rhapson_delta = _newton_rhapson_delta;

        if (abs(bed_slope) < 0.000001){
                bed_slope = 0;
        }

        critical_depth = newtonrhapson(1);
        if (bed_slope != 0){
            normal_depth = newtonrhapson(0);
        }
        else{
            normal_depth = critical_depth + 1;
        }
        if (critical_depth < 0 || normal_depth < 0){
            cout << "El algoritmo Newton-Rhapson no converge" << endl;
            exit(1);
        }
        boundaryconditions();
        rungekutta();
        prepareplot();
    }

    void River::printcalados(){
        /*
        Imprime en pantalla el valor calculado de los calados critico y normal
        */
        double normal = normal_depth;
        double critical = critical_depth;
        double delta = delta_x;
        cout << "Calado normal:\t" << normal << endl;
        cout << "Calado critico\t" << critical << endl;
        cout << "Delta x\t" << delta << endl;
    }

    double River::area(double height){
        /*
        Devuelve el area de la seccion transversal del canal para la altura dada.
        A = (b + my) * y
        */
        return (bottom_width + side_slope * height) * height;
    }

    double River::Pw(double height){
        /*
        Devuelve el perimetro mojado del canal para la altura dada.
        Pw = b + 2y*sqrt(1 + m**2)
        */
        return bottom_width + 2 * height * sqrt(1 + side_slope*side_slope);
    }

    double River::T(double height){
        /*
        Devuelve la anchura de la superficie libre del canal para la altura dada.
        T = b + 2my
        */
        return bottom_width + 2 * side_slope * height;
    }

    double River::fn(double yn){
        /*
        Implicit nonlinear function, fn(yn)
        */
        return (bed_slope >= 0 ? sqrt(bed_slope) : 0) * pow(area(yn), 5.0/3.0) / (manning * pow(Pw(yn), 2.0/3.0)) - discharge;
    }

    double River::dfn(double yn){
        /*
        Implicit nonlinear function derivative, dfn(yn)/dyn
        */
        return (bed_slope >= 0 ? sqrt(bed_slope) : 0) / manning * (pow(area(yn), 2.0/3.0)) *(-4.0/3.0 * sqrt(1 + side_slope*side_slope) + 5.0/3.0 * T(yn) / pow(Pw(yn), 2.0/3.0));
    }

    double River::fc(double yc){
        /*
        Implicit nonlinear function, fc(yc)
        */
        return pow(area(yc), 3.0/2.0) / (sqrt(T(yc))) - discharge/sqrt(gravity);
    }

    double River::dfc(double yc){
        /*
        Implicit nonlinear function derivative, dfc(yc)/dyc
        */
        return -bed_slope * (area(yc) / pow(sqrt(T(yc)), 3.0/2.0)) + 3.00/2.00 * sqrt(area(yc) * T(yc));
    }

    double River::newtonrhapson(bool b){
        /*
        Calcula normal_depth (b=0) critical_depth (b=1) mediante el metodo de Newton-Rhapson
        */

        int i = 0;
        double y = 1;

        if (b){
            double d = fc(y) / dfc(y);
                while (abs(d) > newton_rhapson_delta && i < 100000){
                y = abs(y - 2*d/log(i+2));
                i += 1;
                d = fc(y) / dfc(y);
        }
        return y;
        }

        double d = fn(y) / dfn(y);
        while (abs(d) > newton_rhapson_delta && i < 100000){
            y = abs(y - 2*d/log(i+2));
            i += 1;
            d = fn(y) / dfn(y);
        }
        return y;
    }

    void River::boundaryconditions(){
        /*
        Escoge depth[0] y delta_x basándose en yn e yc
        */

        if (abs(critical_depth - normal_depth) < 0.01){

            //Critical slope
            switch(mode){
                case 1:
                    runge_kutta_step = -delta_x;
                    depth[0] = 1.1 * normal_depth;
                    profile_channel_type = "C1";
                break;

                case 2:
                    runge_kutta_step = -delta_x;
                    depth[0] = 0.5 * (normal_depth + critical_depth);
                    profile_channel_type = "C2";
                break;

                case 3:
                    runge_kutta_step = delta_x;
                    depth[0] = 0.01 * critical_depth;
                    profile_channel_type = "C3";
                break;

                default:
                    cout << "Configuracion no soportada en caso C." << endl;
                    exit(2);
            }
        }
        else{  

            if (critical_depth < normal_depth){
                if (bed_slope <= 0){
                    normal_depth = critical_depth;
                }

                //Mild slope
                switch(mode){
                case 1:
                    //M1
                    //Condiciones río abajo
                    runge_kutta_step = -delta_x;
                    depth[0] = 1.01 * normal_depth;
                    profile_channel_type = "M1";
                break;

                case 2:
                    //M2
                    //Condiciones río abajo:
                    runge_kutta_step = -delta_x;
                    depth[0] =(bed_slope == 0 ? critical_depth*1.0001 : 0.5001 * (normal_depth + critical_depth));
                    //depth[0] = 1.01 * critical_depth
                    (bed_slope == 0) ? profile_channel_type = "H2" : profile_channel_type = "M2";
                break;

                case 3:
                    //M3
                    //Condiciones río arriba:
                    if (bed_slope == 0){
                        runge_kutta_step = delta_x;
                        depth[0] = 0.1 * critical_depth;
                        profile_channel_type = "H3";
                    }
                    else{
                        runge_kutta_step = delta_x;
                        depth[0] = 0.1 * critical_depth;
                        profile_channel_type = "M3";
                    }
                break;

                default:
                    cout << "Configuracion no soportada en caso M." << endl;
                    exit(3);
                }
            }
            else if (critical_depth > normal_depth){
                //Steep slope
                switch(mode){
                    case 1:
                        //S1
                        //Condiciones río arriba:
                        runge_kutta_step = -delta_x;
                        depth[0] = 1.5 * critical_depth;
                        profile_channel_type = "S1";
                    break;

                    case 2:
                        //S2
                        //Condiciones río arriba:
                        runge_kutta_step = delta_x ;
                        depth[0] = 0.5 * (normal_depth + critical_depth);
                        profile_channel_type = "S2";
                        if (bed_slope < 0){
                            depth[0] = 1.5 * critical_depth;
                            profile_channel_type = "A2";
                        }
                    break;

                    case 3:
                        //S3
                        //Condiciones río abajo:
                        runge_kutta_step = delta_x;
                        depth[0] = 0.1 * normal_depth;
                        (bed_slope >= 0) ? profile_channel_type = "S3" : profile_channel_type = "A3";
                    break;

                    default:
                        cout << "Configuracion no soportada en caso S." << endl;
                        exit(4);
                }
            }
            else{
                cout << "Esta configuracion no es soportada: no coincide con ninguna clasificacion" << endl;
                exit(5);
            }  
        }
    }

    double River::diffeq(double y){
        /*
        Ecuacion de evolucion de la profundidad con el desplazamiento para el metodo de Runge_Kutta
        */
       return (bed_slope - manning*manning * discharge*discharge * (pow(Pw(y), 4.0/3.0) / pow(area(y), 10.0/3.0))) / (1 - (discharge*discharge * T(y) / (gravity * pow(area(y),3))));
    }

    void River::rungekutta(){
        /*
        Método de Runge-Kutta de segundo orden para obtener la altura del fluido en los puntos
        */
        
       int i;
       for(i=0; i < steps; i++){                                                                                                        //TODO: Revisar Runge-Kutta
           depth[i+1] = depth[i] + diffeq(depth[i] + 0.5 * diffeq(depth[i]) * runge_kutta_step) * runge_kutta_step;
       }
    }
    
    void River::prepareplot(){
        /*
        Prepara los arrays para los plots
        */

        int i;
        int mult;
        
        if(runge_kutta_step < 0){
            mult = 1;
        }
        else{
            mult = -1;
        }
        
        for(i=0; i < steps; i++){
            positions[i] = i * delta_x;
        }
        
        double corte = positions[steps-1] * bed_slope;
        for(i=0; i < steps; i++){
            bed[i] = corte-positions[i] * bed_slope;
            height[i] = bed[i] + depth[mult*i];
            plot_normal_depth[i] = bed[i] + normal_depth;
            plot_critical_depth[i] = bed[i] + critical_depth;
        }
    }

    int River::escribedatos(){
        /*
        Escribe los datos necesarios para hacer un plot en gnuplot
        */
        FILE *f = fopen("datos.txt","wt");

        if (f == NULL){
            cout << "Error al guardar datos" << endl;
            return 0;
        } 

        int i;

        for(i=0; i < steps ; i++){
        fprintf(f,"%lf %lf %lf %lf %lf\n", positions[i], bed[i], plot_critical_depth[i], plot_normal_depth[i], height[i]);
        }
        fclose(f);

        return 1;
    }

    int River::escribeplot(){
        /*
        Escribe un script de gnuplot para representar los datos
        */
        FILE *f = fopen("plot_perfil.plt","wt");

        if (f == NULL){
            cout << "Error al guardar datos" << endl;
            return 0;
        } 

        int i;
        fprintf(f,"set grid\n set key box opaque\nset xlabel 'Distance (m)'\nset ylabel 'Elevation (m)'\n set title 'Water surface profile (%s)'\nset xrange[0:%lf]\n", profile_channel_type, length);
        fprintf(f,"p \"datos.txt\" u 1:2 w l lc 8 lw 2 t 'Bed', \"datos.txt\" u 1:3 w l lc 7 lw 2 t 'Critical depth', \"datos.txt\" u 1:4 w l lc 2 lw 2 t 'Normal depth', \"datos.txt\" u 1:5 w l lc 22 lw 2 t 'Water surface'\n");
        
        fclose(f);

        return 1;
    }    


//Fin de metodos de clase River


int main(int argc, char **argv){
    
    switch(argc){ //Leemos argumentos por terminal
      case 3:
      {
        River pato(atoi(argv[1]), atoi(argv[2]));
        pato.printcalados();
        pato.escribedatos();
        pato.escribeplot();
      }
      break;

      case 11:
      {
        River pato(atoi(argv[1]), atoi(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), atof(argv[10]));
        pato.printcalados();
        pato.escribedatos();
        pato.escribeplot();
      }
      break;
      default:
          printf("Error: argumentos incorrectos. Los argumentos son los siguientes\n %s steps, mode, discharge, length, bed_slope, side_slope, bottom_width, manning, gravity, newton_rhapson_delta\n",argv[0]);
          exit(-1);
      }

    return 0;
}
