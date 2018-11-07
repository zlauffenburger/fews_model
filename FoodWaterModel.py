from __future__ import division
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt


class FoodWaterModel(object):

    def __init__(self, **kwargs):
        self.P = kwargs['precipitation'] * 1.e-6
        self.PET = kwargs['PET'] * 1.e-6
        self.i_eff = kwargs['irrig_effic']

        self.A = kwargs['total_area']  # Area of the domain (Km^2)

        self.c = kwargs['streamflow_factor']/ self.A
        self.d = kwargs['streamflow_exponent']
        self.w = kwargs['budyko_parameter']
        self.S_ref = kwargs['ref_storage']

        self.tau = kwargs['ces_tau']
        self.beta_land = kwargs['ces_beta_land']
        self.beta_water = kwargs['ces_beta_water']
        self.epsilon = kwargs['ces_epsilon']
        self.rho = (kwargs['ces_sigma'] - 1) / kwargs['ces_sigma']

        self.virt_wat = kwargs['food_trade']

        self.r = kwargs['population_growth_rate']
        self.mu = kwargs['energy_content_crop']
        self.M = kwargs['energy_usage_per_capita']

        self.f_l_ag = kwargs['fraction_cropland']
        self.f_ag_irr = kwargs['fraction_irr_ag']
        self.cwr = kwargs['crop_water_requirement'] * 1.e-6

        self.f_ppl_act = kwargs['fraction_people_active']
        self.f_time_ag = kwargs['fraction_time_agric']
        N = kwargs['init_population']
        S = kwargs['init_storage']
        self.sol_array = np.array([[N], [S]])
        self.time = [0]  # Hold the times at which the model has been evaluated

        self.Q = []
        self.production = []

    def _AET(self, P, Wag=0):
        PWag = P + Wag

        return (1 + (self.PET / PWag) - (1 + (self.PET / PWag) ** self.w) ** (1 / self.w)) * PWag

    def _Q(self, S):
        S = S.clip(min=self.S_ref)
        return self.c * (S - self.S_ref) ** self.d

    def _ces_production_function(self, Lag, Wag, Hag):

        return self.tau * (self.beta_land * Lag ** self.rho
                           + self.beta_water * Wag ** self.rho
                           + (1 - self.beta_land -
                              self.beta_water) * Hag ** self.rho) ** (self.epsilon / self.rho) \
               - self.virt_wat

    def _water_balance(self, S, W_ag):
        return (self.P
                + self.f_l_ag * self.f_ag_irr * W_ag
                - self._Q(S)
                - (1 - self.f_l_ag) * self._AET(self.P)
                - self.f_l_ag * self._AET(self.P, W_ag)
                - self.f_l_ag * W_ag / self.i_eff
                ) * self.A

    def _K(self, L, W, H):
        return self.mu * self._ces_production_function(L, W, H) / self.M

    def _few_system(self, t, x):

        self.f_l_ag = 1 / (1 + np.exp(-0.000000005 * x[0])) - 0.5
        Lag = self.f_l_ag * self.A
        Wag = self.cwr * self.f_ag_irr
        Hag = self.f_ppl_act * self.f_time_ag * x[0]

        dNdt = self.r * x[0] * (1 - x[0] / self._K(Lag, Wag*Lag, Hag))
        dSdt = self._water_balance(x[1], Wag)

        return np.array([dNdt, dSdt]), Wag, Lag

    def integrate(self, x0, tend, dt, *args):
        """
        Integrates Daisy World system
        """

        r = ode(self._few_system)
        r.set_initial_value(x0)

        while r.successful() and r.t < tend:
            sol = r.integrate(r.t + dt)
            self.time.append(r.t)
            self.sol_array = np.append(self.sol_array, np.array(sol).reshape((2, 1)), axis=1)
            #self.Q = np.append(self.Q, )

        return self.time

    def plot_time_series(self, *args):
        t = np.arange(self.sol_array.shape[1])
        print self.sol_array
        time_years = self.time
        plt.plot(time_years, self.sol_array[0, :], 'b', label='Population')
        plt.xlabel('Time [years]', fontsize=14)
        plt.ylabel('People', fontsize=14)
        p = plt.twinx()
        p.plot(time_years, self.sol_array[1, :], 'g', label='Water Storage')
        plt.ylabel('Storage', fontsize=14)
        plt.title('Population growth in the CONUS', fontsize=14)
        plt.legend(loc='best')
        p.legend(loc='best')
        plt.grid()
        plt.show()

    def plot_phase_space(self, x0, dt, tend):
        plt.figure(figsize=(12, 9))
        X0 = x0
        TIME = self.integrate(X0, tend, dt)
        plt.plot(self.sol_array[0, 1:], self.sol_array[1, 1:], 'ro')

        x = np.linspace(self.sol_array[0, 1:].min()*0.8, self.sol_array[0, 1:].max()*1.2, 50)
        y = np.linspace(self.sol_array[1, 1:].min()*0.8, self.sol_array[1, 1:].max()*1.2, 50)

        X1, Y1 = np.meshgrid(x, y)
        X = np.array([X1, Y1])

        DX, Wag, fl_ag = self._few_system(0, [X1, Y1])
        M = np.hypot(DX[0], DX[1])
        M[M == 0] = 1
        #DX /= M

        plt.title('Trajectories and Direction Fields', fontsize=14)
        plt.quiver(X1, Y1, DX[0], DX[1], np.asarray(fl_ag), angles='xy', pivot='mid')
        cbar = plt.colorbar()
        cbar.set_label('Water for Agriculture', fontsize=14)
        cbar.set_clim(np.min(Wag), np.max(Wag))
        plt.text(0.62, 0.9, '$Initial \ N_b=$' + str(x0[0]), fontsize=14)
        plt.text(0.62, 0.85, '$Initial \ N_w=$' + str(x0[1]), fontsize=14)
        plt.xlabel('Population', fontsize=14)
        plt.ylabel('Water Storage', fontsize=14)

        plt.grid()
        plt.show()
