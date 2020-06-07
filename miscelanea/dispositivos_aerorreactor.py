# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:54:28 2020

@author: marco
"""

import numpy as np




# Condiciones de iteración
iter_max = 100
tolerancia = 1e-6


def compresor_tempsalida_rend_adiab(eta, Tet, rel_compresion, gamma, verboso=False):
    """
    Solución de la temperatura de remanso a la salida de un compresor (ventilador):
    Entradas:
        - eta: Rendimiento adiabático del compresor [-]
        - Tet: Temperatura de remanso en la entrada [K]
        - rel_compresion: Relación de compresión [-]
        - Coeficiente de dilatación adiabático (gamma) [-]
            - Valor discreto
            - Función de la temperatura
        - verboso: Si es verdadero: print con numero de iteracion, temperatura y residuo para cada iteracion
    
    Si la gamma es función de la temperatura la solución se itera.
    
    Salida:
        - Temperatura de remanso a la salida [K]
        - Número de iteraciones necesarias para converger
        - Residuo de la temperatura (adimensional)
    """
    if callable(gamma):
        # Itera
        sigue_iterando = True
        Tst = Tet
        niteraciones = 0
        while sigue_iterando:
            niteraciones += 1
            
            Tst_ = Tst
            Tmedia = (Tst + Tet)/2
            gamma_ = gamma(Tmedia)
            
            numerador = rel_compresion ** ((gamma_ - 1)/gamma_) - 1
            Tst = (numerador / eta + 1) * Tet
            
            residuo = np.abs((Tst - Tst_)/Tst_)
            
            if verboso:
                print(niteraciones, Tst, residuo)
            if np.abs(residuo)<tolerancia:
                sigue_iterando = False
                return Tst, niteraciones, residuo
            elif niteraciones==iter_max:
                print('ATENCION: NUMERO ITERACIONES MAXIMAS ALCANZADAS')
                sigue_iterando = False
                return Tst, niteraciones, residuo
                
            
    else:
        numerador = rel_compresion ** ((gamma - 1)/gamma) - 1
        Tst = (numerador / eta + 1) * Tet
        return Tst
    
    
def combustor_tempsalida_dosado(eta, Tet, f, L, cp, verboso=False):
    """
    Solución de la temperatura de remanso a la salida de un combustor:
    Entradas:
        - eta: Rendimiento de la combustión [-]
        - Tet: Temperatura de remanso en la entrada [K]
        - f: dosado [-]
        - L: Poder calorífico del combustible [J/kg]
        - cp: Calor específico a presión constante del aire
            - Valor discreto
            - Función de la temperatura
        - verboso: Si es verdadero: print con numero de iteracion, temperatura y residuo para cada iteracion
    
    Si el cp es función de la temperatura la solución se itera.
    
    Salida:
        - Temperatura de remanso a la salida [K]
        - Número de iteraciones necesarias para converger
        - Residuo de la temperatura (adimensional)
    """
    if callable(cp):
        # Itera
        sigue_iterando = True
        Tst = Tet
        niteraciones = 0
        while sigue_iterando:
            niteraciones += 1
            
            Tst_ = Tst
            cp_ = cp(Tst)
            
            Tst = (eta*f*L + cp_*Tet)/((1+f)*cp_)
            
            residuo = np.abs((Tst - Tst_)/Tst_)
            
            if verboso:
                print(niteraciones, Tst, residuo)
            if np.abs(residuo)<tolerancia:
                sigue_iterando = False
                return Tst, niteraciones, residuo
            elif niteraciones==iter_max:
                print('ATENCION: NUMERO ITERACIONES MAXIMAS ALCANZADAS')
                sigue_iterando = False
                return Tst, niteraciones, residuo
                    
                
    else:
        Tst = (eta*f*L + cp*Tet)/((1+f)*cp)
        return Tst
        

def turbina_tempsalida_acoplamiento(eta_mec, Wconsumido, GTurbina, cp, Het, Tet=500, verboso=False):
    """
    Solución de la temperatura de remanso a la salida de una turbina dada la 
    potencia que debe suministrar, dadas unas pérdidas mecánicas en el eje de
    acoplamiento mecánico y dado un cp:
        -(Wconsumido + Wperdido) = WTurbina
        
    Entradas:
        - eta_mec: Tanto por 1 de pérdidas en eje mecánico [-]
        - Wconsumido: Potencia mecánica consumida [W]
        - GTurbina: Gasto circulante por la turbina [kg/s]
        - cp: Calor específico a presión constante del aire
            - Valor discreto
            - Función de la temperatura
        - Het: entalpía total entrante en la turbina [W]
        Tet: Estimación inicial de la temperatura en entrada de la turbina [K]
            Por defecto se considera temperatura de 500K
        - verboso: Si es verdadero: print con numero de iteracion, temperatura y residuo para cada iteracion
                Por defecto se considera Falso.
    
    Si el cp es función de la temperatura la solución se itera.
    
    Salida:
        - Temperatura de remanso a la salida [K]
        - Número de iteraciones necesarias para converger
        - Residuo de la temperatura (adimensional)
    """
    
    Wconsumido = np.abs(Wconsumido)
    WTurbina = -Wconsumido/(1+eta_mec) # Potencia mecánica requerida a la turbina
    hst = (WTurbina + Het)/GTurbina    # Potencia específica a la salida de turbina
    
    
    if callable(cp):
        # Itera
        sigue_iterando = True
        Tst = Tet
        niteraciones = 0
        
        while sigue_iterando:
            niteraciones += 1
            
            Tst_ = Tst
            Tturb = (Tst + Tet)/2
            cp_ = cp(Tturb)
            
            Tst = hst/cp_
            
            residuo = np.abs((Tst - Tst_)/Tst_)
            if verboso:
                print(niteraciones, Tst, residuo)
            if np.abs(residuo)<tolerancia:
                sigue_iterando = False
                return Tst, niteraciones, residuo
            elif niteraciones==iter_max:
                print('ATENCION: NUMERO ITERACIONES MAXIMAS ALCANZADAS')
                sigue_iterando = False
                return Tst, niteraciones, residuo
                    
                
    else:
        Tst = hst/cp
        return Tst
    
    
def turbina_presionsalida_rend_adiab(eta, pet, rel_temperaturas, gamma, Tturb=500, verboso=False):
    """
    Presión de remanso a la salida de una turbina dada la relación de temperaturas
    salida/entrada y la presión de entrada.

    Entradas:
        - eta: Rendimiento adiabático de turbina [-]
        - pet: Presión de remanso a la entrada [Pa]
        - Coeficiente de dilatación adiabático (gamma) [-]
        - Tturb: Temperatura promedio de la turbina p.ej. (Tet+Tst)/2 [K]
            Por defecto 500K
    

    Salida:
        - Temperatura de remanso a la salida [K]
    """
    if callable(gamma):
        gamma_=gamma(Tturb)
    else:
        gamma_=gamma
        
    pst = pet*((rel_temperaturas - 1)/eta + 1)**(gamma_/(gamma_-1))
    return pst


def tobera_temperaturaestatica_rend_adiab(eta, ps, pet, Tet, gamma):
    
    if callable(gamma):
        gamma_=gamma(Tet)
    else:
        gamma_=gamma
        
    factor = (ps/pet)**((gamma_-1)/gamma_) - 1
    factor = factor * eta + 1
    
    
    Ts = Tet * factor
    return Ts




def test(T4t, cp_air, G5, Wca, H4t):
    sigue_iterando = True
    ii = 0
    T5t = 500
    while sigue_iterando:
        ii += 1
        T5t_ = T5t
        Tturb = (T4t+T5t)/2
        
        T5t = (-Wca + H4t)/G5/cp_air(Tturb)
        
        residuo = np.abs((T5t-T5t_)/T5t_)
        print(ii, T5t, T5t_, residuo, cp_air(Tturb), Tturb)
        if residuo<1e-6:
            sigue_iterando=False
            return T5t
        