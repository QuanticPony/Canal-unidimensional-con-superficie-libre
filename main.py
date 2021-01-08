import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


class River():
    """
    La clase River contiene las funciones que caracterizan el flujo de un canal de superficie abierta
    """

    def __init__(self, 
        # Obligatorios:
            steps,                      # Cantidad de puntos en el eje x
            mode,                       # 1: Por enciama de y_normal. 2: Entre n y_normal e y_critica. 3: Por debajo de y_critica
        # Predefinidos (Se pueden modificar)
            *, discharge=25.0,            # Caudal
            length=4000.0,                # Longitud del canal
            bed_slope=0.001,            # Pendiente del canal. Si es positiva va hacia abajo
            side_slope=1.0,               # Pendiente de los laterales del canal
            bottom_width=5.0,             # Anchura del fondo
            manning=0.025,              # Numero de Manning
            gravity=9.80665,            # Aceleracion de la gravedad
            newton_rhapson_delta=0.001, # Delta para dar por bueno el resultado en el metodo Newton-Rhapson
        # Modo depuracion
            debug=False
        ):
        """
        Iniciamos el objeto `River` con las variables necesarias dadas por el usuario. 
        Rellenamos con valores predefinidos los parámetros no especificados.
        """

        self.steps = steps
        self.discharge = discharge
        self.length = length
        self.bed_slope = bed_slope
        self.side_slope = side_slope
        self.bottom_width = bottom_width
        self.manning = manning
        self.gravity = gravity
        self.newton_rhapson_delta = newton_rhapson_delta

        self._debug_mode = debug
        if debug:
            self._debug_buffer = []
        
        if abs(bed_slope) < 0.000001:
            self.bed_slope = 0


        self.__init_computation_variables()
        self._newtonrhapson()
        self._boundaryconditions(mode)
        self._rungekutta()
        
        self._prepare_plot()


    def __init_computation_variables(self):
        """
        Crea las variables extras necesarias para realizar la computación
        """
        self.runge_kutta_step : float
        self.delta_x = self.length / (self.steps - 1)
        self.positions = np.arange(0, self.length + self.delta_x, self.delta_x)
        self.depth = np.zeros(self.steps)
        self.height = np.zeros(self.steps)
        self.bed = np.zeros(self.steps)
        self.plot_normal_depth = np.zeros(self.steps)
        self.plot_critical_depth = np.zeros(self.steps)
        self.normal_depth: float
        self.critical_depth: float
        self.profile_channel_type : str


    def area(self, height) -> float:
        """
        Devuelve el area de la seccion transversal del canal para la altura dada.
        A = (b + my) * y
        """
        return (self.bottom_width + self.side_slope * height) * height

    def Pw(self, height) -> float:
        """
        Devuelve el perimetro mojado del canal para la altura dada.
        Pw = b + 2y*sqrt(1 + m**2)
        """
        return self.bottom_width + 2 * height * np.sqrt(1 + self.side_slope**2)

    def T(self, height) -> float:
        """
        Devuelve la anchura de la superficie libre del canal para la altura dada.
        T = b + 2my
        """
        return self.bottom_width + 2 * self.side_slope * height

    def _fn(self, yn):
        """
        implicit nonlinear function, fn(yn)
        """

        return np.sqrt(self.bed_slope if self.bed_slope >= 0 else 0) * pow(self.area(yn), 5.0/3.0) / (self.manning * pow(self.Pw(yn), 2.0/3.0)) - self.discharge

    def _fc(self, yc):
        """
        implicit nonlinear function, fc(yc)
        """
        return pow(self.area(yc), 3.0/2.0) / \
            (np.sqrt(self.T(yc))) - self.discharge/np.sqrt(self.gravity)

    def _dfn(self, yn):
        """
        implicit nonlinear function derivative, dfn(yn)/dyn
        """

        return np.sqrt(self.bed_slope if self.bed_slope >= 0 else 0) / self.manning * (pow(self.area(yn), 2.0/3.0)) * \
            (-4.0/3.0 * np.sqrt(1 + self.side_slope**2) + 5.0/3.0 * self.T(yn) / pow(self.Pw(yn), 2.0/3.0))

    def _dfc(self, yc):
        """
        implicit nonlinear function derivative, dfc(yc)/dyc
        """
        
        return -self.bed_slope * (self.area(yc) / pow(np.sqrt(self.T(yc)), 3.0/2.0)) + \
            3.00/2.00 * np.sqrt(self.area(yc) * self.T(yc))

    def _newtonrhapson(self, *, y_c_attempt=0.01, y_n_attempt=0.01):
        """
        Calcula self.normal_depth y self.critical_depth mediante el método de Newton-Rhapson
        """

        def cdiff(y): 
            "y_critica"
            return self._fc(y) / self._dfc(y)
        def ndiff(y): 
            "y_normal"
            return self._fn(y) / self._dfn(y)

        def _debug_calculate(diff, y):
            i = 0
            d = diff(y)
            buffer = []
            while abs(d) > self.newton_rhapson_delta and i < 100000:
                y = abs(y - 2*d/np.math.log(i+2))
                i += 1
                buffer.append(d)
                d = diff(y)
            print(f"Newthon-Rhapson: {diff.__doc__} ha convergido en {i} iteraciones")
            self._debug_buffer.append(buffer)
            return y

        def calculate(diff, y):
            i = 0
            d = diff(y)

            while abs(d) > self.newton_rhapson_delta and i < 100000:
                y = abs(y - 2*d/np.math.log(i+2))
                i += 1
                d = diff(y)
            return y

        if self._debug_mode:
            calculate=_debug_calculate


        # Para obtener la pendiente crítica
        #for i in range(1, 10000):
        #    self.bed_slope += 1e-8
        #    self.critical_depth = calculate(cdiff, y_c_attempt)
        #    self.normal_depth = calculate(ndiff, y_n_attempt)
        #    if abs(self.critical_depth - self.normal_depth) < 0.000001:
        #        print(f"{self.bed_slope}, {self.critical_depth}, {self.normal_depth}")
        #        exit()
            

        self.critical_depth = calculate(cdiff, y_c_attempt)
        if self.bed_slope != 0:
            self.normal_depth = calculate(ndiff, y_n_attempt)
        else:
            self.normal_depth = self.critical_depth + 1

        if self._debug_mode:
            _fig, _ax = plt.subplots()
            _ax.plot(np.arange(len(self._debug_buffer[0])), self._debug_buffer[0], color="red", label="y_critica")
            if self.bed_slope != 0:
                _ax.plot(np.arange(len(self._debug_buffer[1])), self._debug_buffer[1], color="green", label="y_normal")
            
            _fig.legend()
            _ax.set_xlabel("#iteraciones")
            _ax.set_ylabel("delta y (m)")
            plt.show()

        if self.critical_depth<0 or self.normal_depth<0:
            print("El algoritmo Newton-Rhapson no converge")
            exit()


    def _boundaryconditions(self, mode):
        """
        Escoge self.depth[0] y self.delta_x basándose en yn e yc
        """
        if abs(self.critical_depth - self.normal_depth) < 0.01:
            if mode == 1:
                self.runge_kutta_step = self.delta_x 
                self.depth[0] = 1.01 * self.normal_depth
                self.profile_channel_type = "C1"
            
            if mode == 2:
                self.runge_kutta_step = -self.delta_x 
                self.depth[0] = 0.5 * (self.normal_depth + self.critical_depth)
                self.profile_channel_type = "C2"

            if mode == 3:
                self.runge_kutta_step = self.delta_x 
                self.depth[0] = 0.01 * self.critical_depth
                self.profile_channel_type = "C3"
            return

        if self.critical_depth < self.normal_depth :
            if self.bed_slope <= 0:
                self.normal_depth = self.critical_depth

            #Mild slope
            if mode==1:
                #M1
                #Condiciones río abajo
                self.runge_kutta_step = -self.delta_x 
                self.depth[0] = 1.01 * self.normal_depth
                self.profile_channel_type = "M1"
            elif mode==2: 
                #M2
                #Condiciones río abajo:
                self.runge_kutta_step = -self.delta_x
                self.depth[0] = 0.5001 * (self.normal_depth + self.critical_depth) if not self.bed_slope == 0 else self.critical_depth*1.0001
                #self.depth[0] = 1.01 * self.critical_depth
                self.profile_channel_type = "M2" if not self.bed_slope == 0 else "H2"

            elif mode==3:
                #M3
                #Condiciones río arriba:
                if self.bed_slope == 0:
                    self.runge_kutta_step = self.delta_x
                    self.depth[0] = 0.1 * self.critical_depth
                    self.profile_channel_type = "H3"
                else:
                    self.runge_kutta_step = self.delta_x
                    self.depth[0] = 0.1 * self.critical_depth
                    self.profile_channel_type = "M3"

        elif self.critical_depth > self.normal_depth:
            #Steep slope
            if mode==1: 
                #S1
                #Condiciones río arriba:
                self.runge_kutta_step = self.delta_x 
                self.depth[0] = 1.01 * self.critical_depth
                self.profile_channel_type = "S1"
            elif mode==2: 
                #S2
                #Condiciones río arriba:
                self.runge_kutta_step = self.delta_x 
                self.depth[0] = 0.5 * (self.normal_depth + self.critical_depth)
                self.profile_channel_type = "S2"
                if self.bed_slope < 0:
                    self.depth[0] = 1.5 * self.critical_depth
                    self.profile_channel_type = "A2"


            elif mode==3: 
                #S3
                #Condiciones río abajo:
                self.runge_kutta_step = self.delta_x
                self.depth[0] = 0.1 * self.normal_depth
                self.profile_channel_type = "S3" if self.bed_slope > 0 else "A3"
        
        else:
            print('Esta configuracion no es soportada')
            exit()        
        

    def _rungekutta(self):
        """
        Método de Runge-Kutta de segundo orden para obtener la altura del fluido en los puntos
        """
        def diffeq(y): 
            return (self.bed_slope - self.manning**2 * self.discharge**2 * (self.Pw(y)**(4.0/3.0) / self.area(y)**(10.0/3.0))) / \
                (1 - (self.discharge**2 * self.T(y) / (self.gravity * self.area(y)**3)))

        for i in range(1, len(self.depth)) :
            self.depth[i] = self.depth[i-1] + diffeq(self.depth[i-1] + 0.5 * diffeq(self.depth[i-1]) * self.runge_kutta_step) * self.runge_kutta_step


    def _prepare_plot(self):
        """
        Prepara los arrays para los plots
        """
        corte=self.positions[-1] * self.bed_slope
        for i in range(self.steps) :
            self.bed[i] = corte-self.positions[i] * self.bed_slope
            self.height[i] = self.bed[i] + self.depth[-i if self.runge_kutta_step<0 else i]
            self.plot_normal_depth[i] = self.bed[i] + self.normal_depth
            self.plot_critical_depth[i] = self.bed[i] + self.critical_depth

        


if __name__ == '__main__':
    pato = River(4000, mode=1, bed_slope=0.00717747, length=100, debug=False)

    print(f'calado normal: {pato.normal_depth:.3f}')
    print(f'calado critico: {pato.critical_depth:.3f}')
    
    plt.plot(pato.positions, pato.bed, color = 'black', label = 'Bed')
    plt.plot(pato.positions, pato.plot_critical_depth, color = 'red', label = 'Critical depth')
    plt.plot(pato.positions, pato.plot_normal_depth, color = 'green', label = "Normal depth")
    plt.plot(pato.positions, pato.height, color = 'blue', label = 'Water surface')
    plt.legend()
    
    plt.xlabel("Distance (m)")
    plt.ylabel("Elevation (m)")
    plt.title(f"Water surface profile ({pato.profile_channel_type})")
    plt.grid(True)
    plt.xlim(0, pato.length)
    plt.show()