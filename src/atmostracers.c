/*
This source file is part of the atmostracers library, which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/atmostracers
*/

#include "../include/atmostracers.h"
#include <math.h>
#include <stdlib.h>
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_V 0.01801527
#define R (N_A*K_B)
#define R_V (R/M_V)
#define T_0 273.15
#define EPSILON (1e-10)
#define DENSITY_WATER 1024.0
#define GRAVITY 9.8

double ret_sink_velocity(int solid_or_liquid, double radius, double air_density)
{
	double dry_air_kinematic_viscosity = 14.8e-6;
	
	double result = 0;
	// The solid case.
	if (solid_or_liquid == 0)
	{
		result = 2*M_PI*pow(radius, 2)*DENSITY_WATER*GRAVITY/(9*M_PI*air_density*dry_air_kinematic_viscosity);
	}
	
	// The liquid case.
	if (solid_or_liquid == 1)
	{
		result = 2*M_PI*pow(radius, 2)*DENSITY_WATER*GRAVITY/(9*M_PI*air_density*dry_air_kinematic_viscosity);
	}
	
    return result;
}

double ret_phase_trans_heat(int direction, double temperature)
{
    /*
    directions:
    0:  gas to liquid
    1:  gas to solid
    2: liquid to solid
    */
    double result;
    if (direction == 0)
    {
        result = 2257000;
	}
    if (direction == 1)
    {
        result = 2257000 + 333500;
	}
    if (direction == 2)
    {
        result = 333500;
	}
    return result;
}

int calc_h2otracers_source_rates(double mass_source_rates[], double heat_source_rates[], double densities[], double tracer_temperature_densities[], double temperature[], int number_of_scalars, double delta_t)
{
	/*
	It assumes the following order for the constituents:
	ice - liquid water - dry air - water vapour
	*/
    double diff_density, phase_trans_density, saturation_pressure, water_vapour_pressure, solid_temperature, liquid_temperature;
    // loop over all grid boxes
    for (int i = 0; i < number_of_scalars; ++i)
    {
    	// determining the temperature of the ice
    	if (densities[i] < EPSILON)
    	{
    		solid_temperature = T_0;
		}
    	else
    	{
    		solid_temperature = tracer_temperature_densities[i]/densities[i];
		}
		
		// determining the temperature of the liquid water
    	if (densities[number_of_scalars + i] < EPSILON)
    	{
    		liquid_temperature = T_0;
		}
    	else
    	{
    		liquid_temperature = tracer_temperature_densities[number_of_scalars + i]/densities[number_of_scalars + i];
		}
		
		// determining the saturation pressure
        if (temperature[i] >= T_0)
    	{
            saturation_pressure = saturation_pressure_over_water(temperature[i]);
		}
        else
    	{
            saturation_pressure = saturation_pressure_over_ice(temperature[i]);
		}
		
		// determining the water vapour pressure
        water_vapour_pressure = densities[3*number_of_scalars + i]*R_V*temperature[i];
        
    	// the amount of water vapour that the air can still take 
        diff_density = (saturation_pressure - water_vapour_pressure)/(R_V*temperature[i]);
        
        // the case where the air is not over-saturated
        if (diff_density >= 0)
        {
            // temperature >= 0 °C
            if (temperature[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                mass_source_rates[i] = -densities[i]/delta_t;
                
                // The amount of liquid water per volume that will evaporate.
                // In case the air cannot take all the water, not everything will evaporate.
                phase_trans_density = fmin(densities[number_of_scalars + i], diff_density);
                
                // The source rate for the liquid water consists of two terms:
                // 1.) the evaporation
                // 2.) the melting of ice
                mass_source_rates[number_of_scalars + i] = (densities[i] - phase_trans_density)/delta_t;
                
                // the tendency for the water vapour
                mass_source_rates[3*number_of_scalars + i] = phase_trans_density/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[i] = mass_source_rates[i]*ret_phase_trans_heat(2, solid_temperature);
                
                // the heat source rates acting on the liquid water
                heat_source_rates[number_of_scalars + i] =
                // the evaporation
                -phase_trans_density*ret_phase_trans_heat(0, T_0)/delta_t;
            }
            // temperature < 0 °C
            else
            {
            	// Everything that can sublimate will sublimate.
                phase_trans_density = fmin(densities[i], diff_density);
                
                // the tendency for the ice contains two terms
                // 1.) the freezing
                // 2.) the phase transition through sublimation
                mass_source_rates[i] = (densities[number_of_scalars + i] - phase_trans_density)/delta_t;
                
            	// It is assumed that the still present liquid water vanishes within one time step.
                mass_source_rates[number_of_scalars + i] = -densities[number_of_scalars + i]/delta_t;
                
                // the tendency for the water vapour
                mass_source_rates[3*number_of_scalars + i] = phase_trans_density/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[i] = (
                // the freezing
                densities[number_of_scalars + i]*ret_phase_trans_heat(2, solid_temperature)
                // the sublimation
                - phase_trans_density*ret_phase_trans_heat(1, solid_temperature))/delta_t;
                
                // the heat source rates acting on the liquid water
                heat_source_rates[number_of_scalars + i] = 0;
            }
        }
        // the case where the air is over-saturated
        else
        {
        	// the vanishing of water vapour through the phase transition
            mass_source_rates[3*number_of_scalars + i] = diff_density/delta_t;
            // temperature >= 0 °C
            if (temperature[i] >= T_0)
            {
            	// It is assumed that the still present ice vanishes within one time step.
                mass_source_rates[i] = -densities[i]/delta_t;
                
                // The source rate for the liquid water consists of two terms:
                // 1.) the condensation
                // 2.) the melting of ice
                mass_source_rates[number_of_scalars + i] = (-diff_density + densities[i])/delta_t;
                
                // the tendency for water vapour (through condensation)
                mass_source_rates[3*number_of_scalars + i] = diff_density/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[i] = -densities[i]*ret_phase_trans_heat(2, solid_temperature)/delta_t;
                
                // the heat source rates acting on the liquid water
                heat_source_rates[number_of_scalars + i] =
                // it is only affected by the condensation
                -diff_density*ret_phase_trans_heat(1, liquid_temperature)/delta_t;
            }
            // temperature < 0 °C
            else
            {
                // The source rate for the ice consists of two terms:
                // 1.) the resublimation
                // 2.) the melting of ice
                mass_source_rates[i] = (-diff_density + densities[number_of_scalars + i])/delta_t;
                
                // It is assumed that the liquid water disappears within one time step.
                mass_source_rates[number_of_scalars + i] = -densities[number_of_scalars + i]/delta_t;
                
                // the tendency for water vapour (resublimation)
                mass_source_rates[3*number_of_scalars + i] = diff_density/delta_t;
                
                // the heat source rates acting on the ice
                heat_source_rates[i] =
                // the component through the resublimation
                (-diff_density*ret_phase_trans_heat(1, solid_temperature)
                // the component through freezing
                + densities[number_of_scalars + i]*ret_phase_trans_heat(2, solid_temperature))/delta_t;
                
                // the heat source rates acting on the liquid water
                heat_source_rates[number_of_scalars + i] = 0;
            }
        }
        
        // The latent heat does not affect the gas phase.
        heat_source_rates[2*number_of_scalars + i] = 0;
        heat_source_rates[3*number_of_scalars + i] = 0;
    }
    return 0;
}

double water_vapour_density_from_rel_humidity(double rel_humidity, double temperature, double density)
{
    double water_vapour_density = rel_humidity*saturation_pressure_over_water(temperature)/(R_V*temperature);
    return water_vapour_density;
}

double rel_humidity(double abs_humidity, double temperature)
{
	double vapour_pressure = abs_humidity*R_V*temperature;
	double saturation_pressure;
	if (temperature > T_0)
	{
		saturation_pressure = saturation_pressure_over_water(temperature);
	}
	if (temperature <= T_0)
	{
		saturation_pressure = saturation_pressure_over_ice(temperature);
	}
	double result = vapour_pressure/saturation_pressure;
	return result;
}

double saturation_pressure_over_water(double temperature)
{
    double temp_c = temperature - T_0;
    double result = 100*6.112*exp(17.62*temp_c/(243.12 + temp_c));
    return result;
}

double saturation_pressure_over_ice(double temperature)
{
    double temp_c = temperature - T_0;
    double result = 100*6.112*exp(22.46*temp_c/(272.62 + temp_c));
    return result;
}


double ret_c_p_cond(int solid_or_liquid, int subcategory, double temp)
{
    double result;
    if (solid_or_liquid == 0)
    {
        result = 2060;
    }
    if (solid_or_liquid == 1)
    {
        result = 4184;
    }
    return result;
}

double ret_c_v_cond(int solid_or_liquid, int subcategory, double temp)
{
    double result;
    if (solid_or_liquid == 0)
    {
        result = 2060;
    }
    if (solid_or_liquid == 1)
    {
        result = 4184;
    }
    return result;
}


// gas quantities
// ------------------------

/*
gaseous constituents IDs:
0: dry air
1: water vapour
*/

double entropy_constants_gas_lookup(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 2429487178047751925300627872548148580712448.000000; // ((MEAN_MASS_D*exp(5.0/3))/(3*M_PI*H_BAR*H_BAR))
	}
	if (gas_constituent_id == 1)
	{
		result = 1511084890012154487904341578321985168998400.000000; // ((MEAN_MASS_V*exp(5.0/3))/(3*M_PI*H_BAR*H_BAR))
	}
	return result;
}

double mean_particle_masses_gas_lookup(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 0.004810e-23;
	}
	if (gas_constituent_id == 1)
	{
		result = 0.002991e-23;
	}
	return result;
}

double spec_heat_capacities_v_gas_lookup(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 717.942189;
	}
	if (gas_constituent_id == 1)
	{
		result = 1396.475121;
	}
	return result;
}

double spec_heat_capacities_p_gas_lookup(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 1005.0;
	}
	if (gas_constituent_id == 1)
	{
		result = 1858.0;
	}
	return result;
}

double specific_gas_constants_lookup(int gas_constituent_id)
{
	double result = 0;
	if (gas_constituent_id == 0)
	{
		result = 287.057811;
	}
	if (gas_constituent_id == 1)
	{
		result = 461.524879;
	}
	return result;
}






