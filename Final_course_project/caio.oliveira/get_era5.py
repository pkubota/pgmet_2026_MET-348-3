#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Baixa ERA5 (niveis de pressao + superficie) para o caso GOAMAZON/SINGLE.
Requer ~/.cdsapirc configurado (pip install cdsapi)."""
import cdsapi

# Caixa de ~2 graus centrada em (-8, -63): [N, W, S, E]
area = [-6., -65., -10., -61.]
dates = dict(year='2014', month='10', day=['05', '06'],
             time=[f'{h:02d}:00' for h in range(24)])

c = cdsapi.Client()

c.retrieve('reanalysis-era5-pressure-levels', dict(
    product_type='reanalysis', data_format='netcdf', area=area, **dates,
    variable=['geopotential', 'temperature', 'specific_humidity',
              'u_component_of_wind', 'v_component_of_wind', 'vertical_velocity'],
    pressure_level=['1000', '975', '950', '925', '900', '875', '850', '825', '800',
                    '775', '750', '700', '650', '600', '550', '500', '450', '400',
                    '350', '300', '250', '225', '200', '175', '150', '125', '100',
                    '70', '50', '30', '20', '10', '7', '5', '3', '2', '1'],
), 'era5_pl.nc')

c.retrieve('reanalysis-era5-single-levels', dict(
    product_type='reanalysis', data_format='netcdf', area=area, **dates,
    variable=['surface_pressure', 'skin_temperature', 'geopotential',
              'surface_sensible_heat_flux', 'surface_latent_heat_flux'],
), 'era5_sl.nc')

# Obs.: na API nova do CDS o resultado pode vir como .zip; se for o caso, descompacte
# e ajuste os nomes dos arquivos (era5_pl.nc / era5_sl.nc).
