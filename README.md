# Simulación del flujo en ríos
En este repositorio se puede encontrar el código utilizado para realizar el trabajo de la asignatura Física de Fluidos, del tercer curso de la Universidad de Zaragoza.

## Autores
>- Unai Lería Fortea
>- Jose Segovia Burillo
>- Víctor Loras Herrero

## Archivos
- [main.py](main.py): Código utilizado y completamente funcional para la realización del trabajo.
- [cannalflowmodel.cpp](cannalflowmodel.cpp): Con ganas de seguir programando un poquito más, intentamos *traducir* al lenguaje de C++ el código utilizado en el trabajo. Puede haber algún error.

# Instrucciones para usar [main.py](main.py)
1. Se puede ejecutar con un entorno de programación o desde ventana de comandos: *python3 main.py*
1. Se pueden cambiar los parámetros de la simulación en el *'main'* del programa: cambiando los argumentos de inicialización del objeto de la clase *River*. A diferencia del programa traducido a C++, aquí se pueden poner cualquier número de argumentos (y en desorden), poniendo siempre el *keywarg* que se le pasa.
**Son obligatorios** los argumentos *steps* y *mode*
1. El programa sacará él solo las representaciones gráficas haciendo uso de la librería *matplotlib* y no se guarda ningún fichero con los datos.

# Instrucciones para usar [cannalflowmodel.cpp](cannalflowmodel.cpp)
1. Compilar el programa con un compilador para C++ (g++ por ejemplo)
1. Ejecutar el programa ya compilado, escribiendo por linea de comandos los argumentos necesarios, que pueden escribirse de dos maneras
    1. ./exe steps, mode   -> Tendrá predefinidos los demás args *(ver código)*
    1. ./exe steps, mode, discharge, length, bed_slope, side_slope, bottom_width, manning, gravity, newton_rhapson_delta 
1. El programa creará dos ficheros de texto:
    1. *'datos.txt'* contiene la información para la representación
    1. *'plot_perfil.plt'* es un script de *gnuplot* que grafica los datos de *'datos.txt'*, obteniendo el pefil simulado
    
**Importante**: *el programa sobrescribirá automáticamente los datos cotenidos anteriormente en 'datos.txt'*
