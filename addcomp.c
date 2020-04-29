#include "/lib/surface/include/surface.h"
#include "addcomp.h"
#include <math.h>
#include <stdlib.h>
#define N_A (6.0221409e23)
#define K_B (1.380649e-23)
#define M_V 0.01801527
#define R (N_A*K_B)
#define R_V (R/M_V)
#define T_0 273.15

double ret_sink_velocity(long solid_or_liquid, long subcategory, double radius)
{
    return -0.1;
}

double ret_phase_trans_heat(short direction, double temperature)
{
    /*
    directions:
    0:  gas to liquid
    1:  gas to solid
    2: liquid to solid
    */
    double result;
    if (direction == 0)
        result = 2257000;
    if (direction == 1)
        result = 2257000 + 333500;
    if (direction == 2)
        result = 333500;
    return result;
}

int calc_add_comp_source_rates(double mass_source_rates[], double heat_source_rates[], double densities[], double comp_temps[], double temperature[], long number_of_add_comps, long number_of_scalars, double delta_t)
{
    double *saturation_pressure = malloc(number_of_scalars*sizeof(double));
    double *water_vapour_pressure = malloc(number_of_scalars*sizeof(double));
    double diff_density, phase_trans_density;
    for (int i = 0; i < number_of_scalars; ++i)
    {
        if (temperature[i] > T_0)
            saturation_pressure[i] = saturation_pressure_over_water(temperature[i]);
        else
            saturation_pressure[i] = saturation_pressure_over_ice(temperature[i]);
        water_vapour_pressure[i] = densities[2*number_of_scalars + i]*R_V*temperature[i];
        if (water_vapour_pressure[i] <= saturation_pressure[i])
        {
            diff_density = (saturation_pressure[i] - water_vapour_pressure[i])/(R_V*temperature[i]);
            if (temperature[i] >= T_0)
            {
                mass_source_rates[i] = -densities[i]/delta_t;
                phase_trans_density = fmin(densities[number_of_scalars + i], diff_density);
                mass_source_rates[number_of_scalars + i] = (densities[i] - phase_trans_density)/delta_t;
                mass_source_rates[2*number_of_scalars + i] = phase_trans_density/delta_t;
                heat_source_rates[i] = mass_source_rates[i]*ret_phase_trans_heat(1, comp_temps[i]);
                heat_source_rates[number_of_scalars + i] = mass_source_rates[number_of_scalars + i]*ret_phase_trans_heat(1, T_0);
                heat_source_rates[2*number_of_scalars + i] = 0;
            }
            else
            {
                phase_trans_density = fmin(densities[i], diff_density);
                mass_source_rates[i] = (densities[number_of_scalars + i] - phase_trans_density)/delta_t;
                mass_source_rates[number_of_scalars + i] = -densities[number_of_scalars + i]/delta_t;
                mass_source_rates[2*number_of_scalars + i] = phase_trans_density/delta_t;
                heat_source_rates[i] = (densities[number_of_scalars + i]*ret_phase_trans_heat(2, comp_temps[i]) - phase_trans_density*ret_phase_trans_heat(1, comp_temps[i]))/delta_t;
                heat_source_rates[number_of_scalars + i] = mass_source_rates[number_of_scalars + i]*ret_phase_trans_heat(1, comp_temps[number_of_scalars + i]);
                heat_source_rates[2*number_of_scalars + i] = 0;
            }
        }
        else
        {
            diff_density = (water_vapour_pressure[i] - saturation_pressure[i])/(R_V*temperature[i]);
            mass_source_rates[2*number_of_scalars + i] = -diff_density/delta_t;
            if (temperature[i] >= T_0)
            {
                mass_source_rates[i] = -densities[i]/delta_t;
                mass_source_rates[number_of_scalars + i] = (diff_density + densities[i])/delta_t;
                mass_source_rates[2*number_of_scalars + i] = -diff_density/delta_t;
                heat_source_rates[i] = -densities[i]*ret_phase_trans_heat(2, comp_temps[i])/delta_t;
                heat_source_rates[number_of_scalars + i] = mass_source_rates[number_of_scalars + i]*ret_phase_trans_heat(1, comp_temps[number_of_scalars + i]);
                heat_source_rates[2*number_of_scalars + i] = 0;
            }
            else
            {
                mass_source_rates[i] = (diff_density + densities[number_of_scalars + i])/delta_t;
                mass_source_rates[number_of_scalars + i] = -densities[number_of_scalars + i]/delta_t;
                mass_source_rates[2*number_of_scalars + i] = -diff_density/delta_t;
                heat_source_rates[i] = (diff_density*ret_phase_trans_heat(1, comp_temps[i]) + densities[number_of_scalars + i]*ret_phase_trans_heat(2, comp_temps[i]))/delta_t;
                heat_source_rates[number_of_scalars + i] = mass_source_rates[number_of_scalars + i]*ret_phase_trans_heat(1, comp_temps[number_of_scalars + i]);
                heat_source_rates[2*number_of_scalars + i] = 0;
            }
        }
    }
    free(saturation_pressure);
    free(water_vapour_pressure);
    for (int i = 0; i < number_of_scalars*number_of_add_comps; ++i)
        heat_source_rates[i] = 0;
    return 0;
}

double water_vapour_density_from_rel_humidity(double rel_humidity, double temperature, double density)
{
    double water_vapour_density = rel_humidity*saturation_pressure_over_water(temperature)/(R_V*temperature);
    return water_vapour_density;
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


double ret_c_p_cond(long solid_or_liquid, long subcategory, double temp)
{
    double result;
    if (solid_or_liquid == 0)
        result = 2060;
    if (solid_or_liquid == 1)
        result = 4184;
    return result;
}







