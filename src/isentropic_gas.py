#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 18:02:20 2017

Funciones de corriente isoentrópica:
    
Contenido:
----------
    +Mach: Número de Mach de la corriente

    +Tt_vel: Temperatura de remanso en función de la velocidad
             y la temperatura estática

    +Tt_Mach: Temperatura de remanso en función de la velocidad
             y la temperatura estática
            
    +Tt_Test: Relación temperatura de remanso y estática

    +Pt_Pest: Relación presión de remanso y estática
    
    +Pt_Mach: Presión de remanso en función de la velocidad
              y la presión estática
              
Dependencias:
-------------
    + Numpy

@author: marcosrodriguezjimenez
"""

import numpy as np


def configuracion(modelo_gas):
    if modelo_gas is 'ideal':
        cp_a = aire.gas_ideal.cp_a
        gamma_gas = aire.gas_ideal.gamma_a
        return cp_a, gamma_gas
    elif modelo_gas is 'real':
        aire.gas_real()
    else:
        return 'Modelo desconocido'


def Mach(vel, T_est):
    """
    Función Mach:
    -------------
    
    Obtención del número de Mach.
    
    +Entradas:
        vel: Velocidad de la corriente en (m/s)
        T_est: Temperatura Estática (K)
    + Salidas:
        Mach: Número de Mach de a corriente (Adim)
        
    """
    R_g = 287 #J/kg/k
    gamma_gas = 1.4  #Asumimos aire
    Mach = vel /  np.sqrt(gamma_gas * T_est * R_g)
    return Mach


def Tt_vel(T_est, vel):
    """
    Función Tt_vel:
    ---------------
    
    Temperatura de remanso en función de la velocidad 
    y temperatura estática.
    
    + Entradas:
        T_est: Temperatura Estática (K)
        vel: Velocidad de la corriente en (m/s)
    + Salidas:
        Temperatura de remanso de la corriente (K)
    """
    cp = 1004 #J/kg/K
    return T_est + vel**2/cp/2
    
    
def Tt_Mach(T_est, Mach):
    """
    Función Tt_Mach:
    ---------------
    
    Temperatura de remanso en función del número de Mach
    de la corriente y temperatura estática
    
    + Entradas:
        T_est: Temperatura Estática (K)
        Mach: Mach de corriente (Adim)
    + Salidas:
        Temperatura de remanso de la corriente (K)
    """
    gamma_gas = 1.4  #Asumimos aire
    return T_est*(1+(gamma_gas-1)/2*Mach**2)


def Tt_Test(Mach):
    """
    Función Tt_Test:
    ---------------
    
    Relación de temperatura de remanso y estática
    en función del número de Mach de la corriente
    
    + Entradas:
        Mach: Mach de corriente (Adim)
    + Salidas:
        Relación de temperatura remanso y estática
    """
    gamma_gas = 1.4  #Asumimos aire
    return (1+(gamma_gas-1)/2*Mach**2)


def Pt_Pest(Mach):
    """
    Función Pt_Pest:
    ----------------
    
    Relación entre presión de remanso y estática
    en función del número de Mach de la corriente.
    
    + Entradas:
        Mach: Mach de corriente (Adim)
    + Salidas:
        Relación de presión remanso y estática
    """
    gamma_gas = 1.4 
    return Tt_Test(Mach)**(gamma_gas/(gamma_gas - 1))


def Pt_Mach(Pest, Mach):
    """
    Función Pt_Mach:
    ----------------
    
    Presión de remanso en función de la presión estática
    y el número de Mach de la corriente.
    
    + Entradas:
        Pest: Presión estática de la corriente (Pa)
        Mach: Mach de corriente (Adim)
    + Salidas:
        Presión de remanso de la corriente (Pa)
    """
    gamma_gas = 1.4 
    return Pest * Tt_Test(Mach)**(gamma_gas/(gamma_gas - 1))
    
    
    
    
    
    