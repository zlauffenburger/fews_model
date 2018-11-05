import numpy as np
import FoodWaterModel as fw

params = {

    'precipitation': 767,  # Precipitation in mm/year
    'PET': 670,  # Potential evapotranspiration in mm/year
    'irrig_effic': 0.8,
    'streamflow_factor': 0.00002,
    'streamflow_exponent': 1.9,
    'budyko_parameter': 1.5,
    'ref_storage': 10.0,  # km^3
    'ces_tau': 2500000.0,
    'ces_beta_land': 0.2,
    'ces_beta_water': 0.5,
    'ces_epsilon': 0.8,
    'ces_sigma': 0.5,
    'food_trade': -5000,

    'energy_content_crop': 13000,  # kJ/kg metabolizable energy
    'energy_usage_per_capita': 10000,  # kJ per adult
    'total_area': 7600000,  # km2
    'fraction_cropland': 0.028,
    'fraction_irr_ag': 0.2,
    'crop_water_requirement': 1000,
    'fraction_people_active': 0.3,
    'fraction_time_agric': 0.05,
    'population_growth_rate': 1.3,

    'init_population': 100000000,
    'init_storage': 60000,
}

x0 = np.array([params['init_population'], params['init_storage']])
t0 = 0  # Initial Time [years]
dt = 1.  # Time step [years]
tend = 100  # End time [years]


few = fw.FoodWaterModel(**params)
few.integrate(x0, tend, dt)

few.plot_time_series()
few.plot_phase_space(x0, dt, tend)
