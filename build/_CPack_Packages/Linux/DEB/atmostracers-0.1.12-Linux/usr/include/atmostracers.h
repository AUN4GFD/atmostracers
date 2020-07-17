/*
This source file is part of the atmostracers library, which is released under the MIT license.
Github repository: https://github.com/MHBalsmeier/atmostracers
*/

double ret_sink_velocity(int, int, double);
int calc_h2otracers_source_rates(double[], double[], double[], double[], double[], int, int, double);
double water_vapour_density_from_rel_humidity(double, double, double);
double saturation_pressure_over_water(double);
double saturation_pressure_over_ice(double);
double ret_c_p_cond(int, int, double);
double ret_c_v_cond(int, int, double);
double ret_phase_trans_heat(int, double);
double rel_humidity(double, double);
