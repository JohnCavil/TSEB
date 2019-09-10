# This file is part of pyTSEB for running different TSEB models
# Copyright 2016 Hector Nieto and contributors listed in the README.md file.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on Apr 6 2015
@author: Hector Nieto (hnieto@ias.csic.es)

Modified on Jan 27 2016
@author: Hector Nieto (hnieto@ias.csic.es)

DESCRIPTION
===========
This package contains the main routines inherent of Two Source Energy Balance `TSEB` models.
Additional functions needed in TSEB, such as computing of net radiation or estimating the
resistances to heat and momentum transport are imported.

* :doc:`netRadiation` for the estimation of net radiation and radiation partitioning.
* :doc:`ClumpingIndex` for the estimatio of canopy clumping index.
* :doc:`meteoUtils` for the estimation of meteorological variables.
* :doc:`resistances` for the estimation of the resistances to heat and momentum transport.
* :doc:`MOsimilarity` for the estimation of the Monin-Obukhov length and MOST-related variables.

PACKAGE CONTENTS
================

TSEB models
-----------
* :func:`TSEB_2T` TSEB using derived/measured canopy and soil component temperatures.
* :func:`TSEB_PT` Priestley-Taylor TSEB using a single observation of composite radiometric temperature.
* :func:`DTD` Dual-Time Differenced TSEB using composite radiometric temperatures at two times: early morning and near afternoon.

OSEB models
-----------
* :func:`OSEB`. One Source Energy Balance Model.
* :func:`OSEB_BelowCanopy`. One Source Energy Balance Model of baresoil/understory beneath a canopy.
* :func:`OSEB_Canopy`. One Source Energy Balance Model of a very dense canopy, i.e. `Big-leaf` model.

Ancillary functions
-------------------
* :func:`CalcFthetaCampbell`. Gap fraction estimation.
* :func:`CalcG_TimeDiff`. Santanello & Friedl (2003) [Santanello2003]_ soil heat flux model.
* :func:`CalcG_Ratio`. Soil heat flux as a fixed fraction of net radiation [Choudhury1987]_.
* :func:`CalcH_C`. canopy sensible heat flux in a parallel resistance network.
* :func:`CalcH_C_PT`. Priestley- Taylor Canopy sensible heat flux.
* :func:`CalcH_DTD_parallel`. Priestley- Taylor Canopy sensible heat flux for DTD and resistances in parallel.
* :func:`CalcH_DTD_series`. Priestley- Taylor Canopy sensible heat flux for DTD and resistances in series.
* :func:`CalcH_S`. Soil heat flux with resistances in parallel.
* :func:`CalcT_C`. Canopy temperature form composite radiometric temperature.
* :func:`CalcT_C_Series.` Canopy temperature from canopy sensible heat flux and resistance in series.
* :func:`CalcT_CS_Norman`. Component temperatures from dual angle composite radiometric tempertures.
* :func:`CalcT_CS_4SAIL`. Component temperatures from dual angle composite radiometric tempertures. Using 4SAIl for the inversion.
* :func:`Get4SAILEmissionParam`. Effective surface reflectance, and emissivities for soil and canopy using 4SAIL.
* :func:`CalcT_S`. Soil temperature from form composite radiometric temperature.
* :func:`CalcT_S_Series`. Soil temperature from soil sensible heat flux and resistance in series.
'''

import meteoUtils as met
import resistances as res
import MOsimilarity as MO
import netRadiation as rad
import ClumpingIndex as CI
#==============================================================================
# List of constants used in TSEB model and sub-routines   
#==============================================================================
#Change threshold in  Monin-Obukhov lengh to stop the iterations
L_thres=0.00001
# Change threshold in  friction velocity to stop the iterations
u_thres=0.00001
# mimimun allowed friction velocity    
u_friction_min=0.01;
#Maximum number of interations
ITERATIONS=10
# kB coefficient
kB=0.0
#Stephan Boltzmann constant (W m-2 K-4)
sb=5.670373e-8 

def TSEB_2T(Tc,Ts,Ta_K,u,ea,p,Sn_C,Sn_S,Lsky,
            LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
            leaf_width=0.1,z0_soil=0.01,alpha_PT=1.26,x_LAD=1.0,f_c=1.0,f_g=1.0,wc=1.0,
            Resistance_flag=0,calcG_params=[[1],0.35], UseL=False):
                
    ''' TSEB using component canopy and soil temperatures.

    Calculates the turbulent fluxes by the Two Source Energy Balance model 
    using canopy and soil component temperatures that were derived or measured
    previously.
    
    Parameters
    ----------
    Ts : float
        Soil Temperature (Kelvin).
    Tc : float
        Canopy Temperature (Kelvin).
    Ta_K : float 
        Air temperature (Kelvin).
    u : float 
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    Sn_C : float
        Canopy net shortwave radiation (W m-2).
    Sn_S : float
        Soil net shortwave radiation (W m-2).
    Lsky : float
        Downwelling longwave radiation (W m-2)
    LAI : float
        Effective Leaf Area Index (m2 m-2).
    hc : float
        Canopy height (m).
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    zu : float
        Height of measurement of windspeed (m).
    zt : float
        Height of measurement of air temperature (m).
    leaf_width : float, optional
        average/effective leaf width (m).
    z0_soil : float, optional
        bare soil aerodynamic roughness length (m).
    alpha_PT : float, optional
        Priestley Taylor coeffient for canopy potential transpiration, 
        use 1.26 by default.
    Resistance_flag : int, optional
        Flag to determine which Resistances R_x, R_s model to use.
        
            * 0 [Default] Norman et al 1995 and Kustas et al 1999. 
            * 1 : Choudhury and Monteith 1988.
            * 2 : McNaughton and Van der Hurk 1995.   

    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
        
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_soil, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.CalcG_TimeDiff`).
    UseL : float or None, optional
        If included, its value will be used to force the Moning-Obukhov stability length.
    
    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    T_AC : float
        Air temperature at the canopy interface (Kelvin).
    LE_C : float
        Canopy latent heat flux (W m-2).
    H_C : float
        Canopy sensible heat flux (W m-2).
    LE_S : float
        Soil latent heat flux (W m-2).
    H_S : float
        Soil sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    R_s : float
        Soil aerodynamic resistance to heat transport (s m-1).
    R_x : float
        Bulk canopy aerodynamic resistance to heat transport (s m-1).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    u_friction : float
        Friction velocity (m s-1).
    L : float
        Monin-Obuhkov length (m).
    n_iterations : int
        number of iterations until convergence of L.
        
    References
    ----------
    .. [Kustas1997] Kustas, W. P., and J. M. Norman (1997), A two-source approach for estimating
        turbulent fluxes using multiple angle thermal infrared observations,
        Water Resour. Res., 33(6), 1495-1508,
        http://dx.doi.org/10.1029/97WR00704.
    '''
    
    import numpy as np

    # Convert float scalars into numpy arrays and check parameters size
    Tc = np.asarray(Tc)
    [Ts,Ta_K,u,ea,p,Sn_C,Sn_S,Lsky,LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
         leaf_width,z0_soil,alpha_PT,x_LAD,f_c,f_g,wc,calcG_array] = map(_CheckDefaultParameterSize,
         [Ts,Ta_K, u,ea,p,Sn_C,Sn_S,Lsky,LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
            leaf_width,z0_soil,alpha_PT,x_LAD,f_c,f_g,wc,calcG_params[1]], [Tc]*24)
    # Create the output variables
    [flag, T_ac, Ln_S,Ln_C, LE_c,H_c,LE_s,H_s,G,R_s,R_x,R_a] = [np.zeros(Ts.shape) for i in range(12)]
         
    # iteration of the Monin-Obukhov length
    if type(UseL)==bool:
        L = np.asarray(np.zeros(Ts.shape)+np.inf)     # Initially assume stable atmospheric conditions and set variables for 
        max_iterations=ITERATIONS
    else: # We force Monin-Obukhov lenght to the provided array/value
        L=np.asarray(np.ones(Ts.shape)*UseL)
        max_iterations=1 # No iteration
  
    # Calculate the general parameters
    rho_a= met.CalcRho(p, ea, Ta_K)  # Air density
    Cp = met.CalcC_p(p, ea)  # Heat capacity of air
    z_0H=res.CalcZ_0H(z_0M, kB=kB) # Roughness length for heat transport
    
    # Calculate LAI dependent parameters for dataset where LAI > 0
    F = np.asarray(LAI/f_c) # Real LAI    
    # Calculate LAI dependent parameters for dataset where LAI > 0
    omega0=CI.CalcOmega0_Kustas(LAI, f_c,x_LAD=x_LAD, isLAIeff=True)
                
    # And the net longwave radiation
    Ln_C, Ln_S = rad.CalcLnKustas(Tc, Ts, Lsky, LAI, emisVeg, emisGrd)
    
    #Compute Net Radiation
    Rn_soil = Sn_S + Ln_S
    Rn_veg = Sn_C + Ln_C
    Rn = Rn_soil+Rn_veg
    
    #Compute Soil Heat Flux
    i=np.ones(Rn_soil.shape,dtype=bool)
    G[i] = CalcG([calcG_params[0], calcG_array], Rn_soil,i)
    # iteration of the Monin-Obukhov length
    u_friction = MO.CalcU_star(u, zu, L, d_0, z_0M)
    u_friction = np.asarray(np.maximum(u_friction_min, u_friction))
    L_old = np.ones(Tc.shape)
    L_diff = np.asarray(np.ones(Tc.shape)*np.inf)

    # Outer loop for estimating stability. 
    # Stops when difference in consecutives L is below a given threshold
    for n_iterations in range(max_iterations):
        flag[np.isnan(Tr_K)]=255
        i = flag!=255
        if np.all(L_diff < L_thres):
            print("Finished interation with a max. L diff: "+str(np.max(L_diff)))
            break
        
        print("Iteration "+str(n_iterations)+", max. L diff: "+str(np.max(L_diff)))
        
        i = np.logical_and(L_diff >= L_thres, flag!=255) 
        flag[i] = 0
        
        # Calculate the aerodynamic resistance
        R_a[i] = res.CalcR_A(zt[i], u_friction[i], L[i], d_0[i], z_0H[i])
        # Calculate soil and canopy resistances
        U_C = MO.CalcU_C_star (u_friction[i], hc[i], d_0[i], z_0M[i], L=L[i])
        if Resistance_flag == 0:
            u_S=MO.CalcU_Goudriaan (U_C, hc[i], omega0[i]*F[i], leaf_width[i], z0_soil[i]) # Clumped vegetation enhanced wind speed for the soil surface 
            u_d_zm=MO.CalcU_Goudriaan (U_C, hc[i], F[i], leaf_width[i], d_0[i]+z_0M[i]) # Wind speed is highly attenuated within the canopy volume
            R_x[i]=res.CalcR_X_Norman(LAI[i], leaf_width[i], u_d_zm) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
            R_s[i]=res.CalcR_S_Kustas(u_S, Ts[i]-Tc[i])
        elif Resistance_flag ==1:
            R_x[i]=res.CalcR_X_Choudhury(U_C, LAI[i],leaf_width[i]) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
            R_s[i]=res.CalcR_S_Choudhury(u_friction[i],hc[i],z_0M[i],d_0[i],zu[i],z0_soil[i])
        elif Resistance_flag ==2:
            R_x[i]=res.CalcR_X_McNaughton(LAI[i], leaf_width[i], u_friction[i]) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
            R_s[i]=res.CalcR_S_McNaughton(u_friction[i])
        elif Resistance_flag ==3:
            alpha_k=MO.CalcA_Goudriaan(hc[i], omega0[i]*F[i], leaf_width[i]) # Clumped vegetation enhanced wind speed for the soil surface 
            alpha_prime=MO.CalcA_Goudriaan(hc[i], F[i], leaf_width[i]) # Wind speed is highly attenuated within the canopy volume
            R_x[i]=res.CalcR_X_Choudhury(U_C, LAI[i],leaf_width[i],alpha_prime=alpha_prime) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
            R_s[i]=res.CalcR_S_Choudhury(u_friction[i],hc[i],z_0M[i],d_0[i],zu,z0_soil[i],alpha_k=alpha_k)
        else:
            u_S=MO.CalcU_Goudriaan (U_C, hc[i], omega0[i]*F[i], leaf_width[i], z0_soil[i]) # Clumped vegetation enhanced wind speed for the soil surface 
            u_d_zm=MO.CalcU_Goudriaan (U_C, hc[i], F[i], leaf_width[i], d_0[i]+z_0M[i]) # Wind speed is highly attenuated within the canopy volume
            R_x[i]=res.CalcR_X_Norman(LAI[i], leaf_width[i], u_d_zm) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
            R_s[i]=res.CalcR_S_Kustas(u_S, Ts[i]-Tc[i])
        R_s=np.asarray(np.maximum( 1e-3,R_s))
        R_x=np.asarray(np.maximum( 1e-3,R_x))
        R_a=np.asarray(np.maximum( 1e-3,R_a))
        
        # Compute air temperature at the canopy interface
        T_ac[i] = ((Ta_K[i]/R_a[i] + Ts[i]/R_s[i] + Tc[i]/R_x[i])
                    /(1/R_a[i] + 1/R_s[i] + 1/R_x[i]))
        T_ac = np.asarray(np.maximum(1e-3,T_ac))
        
        # Calculate canopy sensible heat flux (Norman et al 1995)
        H_c[i] = rho_a[i]*Cp[i]*(Tc[i]-T_ac[i])/R_x[i]
        # Assume no condensation in the canopy (LE_c<0)
        noC = np.logical_and(i, H_c > Rn_veg)
        H_c[noC] = Rn_veg[noC]
        flag[noC] = 1
        # Assume no thermal inversion in the canopy
        noI = np.logical_and.reduce((i, H_c < CalcH_C_PT(Rn_veg, 1.0, Ta_K, p, Cp, alpha_PT), Rn_veg > 0))
        H_c[noI] = 0
        flag[noI] = 2
            
        # Calculate soil sensible heat flux (Norman et al 1995)
        H_s[i] = rho_a[i]*Cp[i]*(Ts[i]-T_ac[i])/R_s[i]
        # Assume that there is no condensation in the soil (LE_s<0)
        noC = np.logical_and.reduce((i, H_s > Rn_soil-G, (Rn_soil-G) > 0))        
        H_s[noC] = Rn_soil[noC] - G[noC]
        flag[noC] = 3
        # Assume no thermal inversion in the soil
        noI = np.logical_and.reduce((i, H_s < 0, Rn_soil-G > 0))        
        H_s[noI] = 0
        flag[noI] = 4     
                            
        # Evaporation Rate (Kustas and Norman 1999)
        H = np.asarray(H_s+H_c)
        LE = np.asarray(Rn-G-H)
        # Now L can be recalculated and the difference between iterations derived
        if type(UseL)==bool:  
            L[i] = MO.CalcL (u_friction[i], Ta_K[i], rho_a[i], Cp[i], H[i], LE[i])
            L_diff = np.asarray(np.fabs(L-L_old)/np.fabs(L_old))
            L_diff[np.isnan(L_diff)] = np.inf
            L_old=np.array(L)
            L_old[L_old==0] = 1e-36 
                
            # Calculate again the friction velocity with the new stability correctios        
            u_friction[i] = MO.CalcU_star(u[i], zu[i], L[i], d_0[i], z_0M[i])
            u_friction = np.asarray(np.maximum(u_friction_min, u_friction))
        
    # Compute soil and canopy latent heat fluxes
    LE_s = Rn_soil-G-H_s
    LE_c = Rn_veg-H_c
    
    (flag, T_ac,Ln_S, Ln_C, LE_c,H_c,LE_s,H_s,G,R_s,R_x,R_a,u_friction, L,
         n_iterations)=map(np.asarray,(flag, T_ac,Ln_S,Ln_C, LE_c,H_c,LE_s,H_s,
                            G,R_s,R_x,R_a,u_friction, L,n_iterations))
    
    return flag, T_ac,Ln_S,Ln_C, LE_c,H_c,LE_s,H_s,G,R_s,R_x,R_a,u_friction, L,n_iterations

def  TSEB_PT(Tr_K,vza,Ta_K,u,ea,p,Sn_C,Sn_S,Lsky,
            LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
            leaf_width=0.1,z0_soil=0.01,alpha_PT=1.26,x_LAD=1,f_c=1.0,f_g=1.0,wc=1.0,
            Resistance_flag=0,calcG_params=[[1],0.35], UseL=False):
    '''Priestley-Taylor TSEB

    Calculates the Priestley Taylor TSEB fluxes using a single observation of
    composite radiometric temperature and using resistances in series.
    
    Parameters
    ----------
    Tr_K : float
        Radiometric composite temperature (Kelvin).
    vza : float
        View Zenith Angle (degrees).
    Ta_K : float 
        Air temperature (Kelvin).
    u : float 
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    Sn_C : float
        Canopy net shortwave radiation (W m-2).
    Sn_S : float
        Soil net shortwave radiation (W m-2).
    Lsky : float
        Downwelling longwave radiation (W m-2).
    LAI : float
        Effective Leaf Area Index (m2 m-2).
    hc : float
        Canopy height (m).
    emisVeg : float
        Leaf emissivity.
    emisGrd : flaot
        Soil emissivity.
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    zu : float
        Height of measurement of windspeed (m).
    zt : float
        Height of measurement of air temperature (m).
    leaf_width : float, optional
        average/effective leaf width (m).
    z0_soil : float, optional
        bare soil aerodynamic roughness length (m).
    alpha_PT : float, optional
        Priestley Taylor coeffient for canopy potential transpiration, 
        use 1.26 by default.
    x_LAD : float, optional
        Campbell 1990 leaf inclination distribution function chi parameter.
    f_c : float, optional
        Fractional cover.
    f_g : float, optional
        Fraction of vegetation that is green.
    wc : float, optional
        Canopy width to height ratio.
    Resistance_flag : int, optional
        Flag to determine which Resistances R_x, R_s model to use.
        
            * 0 [Default] Norman et al 1995 and Kustas et al 1999. 
            * 1 : Choudhury and Monteith 1988.
            * 2 : McNaughton and Van der Hurk 1995.   

    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
        
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_soil, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.CalcG_TimeDiff`).
    UseL : float or None, optional
        If included, its value will be used to force the Moning-Obukhov stability length.
    
    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    Ts : float
        Soil temperature  (Kelvin).
    Tc : float
        Canopy temperature  (Kelvin).
    T_AC : float
        Air temperature at the canopy interface (Kelvin).
    L_nS : float
        Soil net longwave radiation (W m-2)
    L_nC : float
        Canopy net longwave radiation (W m-2)
    LE_C : float
        Canopy latent heat flux (W m-2).
    H_C : float
        Canopy sensible heat flux (W m-2).
    LE_S : float
        Soil latent heat flux (W m-2).
    H_S : float
        Soil sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    R_s : float
        Soil aerodynamic resistance to heat transport (s m-1).
    R_x : float
        Bulk canopy aerodynamic resistance to heat transport (s m-1).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    u_friction : float
        Friction velocity (m s-1).
    L : float
        Monin-Obuhkov length (m).
    n_iterations : int
        number of iterations until convergence of L.

    References
    ----------
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293,
        http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    .. [Kustas1999] William P Kustas, John M Norman, Evaluation of soil and vegetation heat
        flux predictions using a simple two-source model with radiometric temperatures for
        partial canopy cover, Agricultural and Forest Meteorology, Volume 94, Issue 1,
        Pages 13-29,
        http://dx.doi.org/10.1016/S0168-1923(99)00005-2.
    '''
    
    import numpy as np   
    print(alpha_PT)
    # Convert input float scalars to arrays and parameters size
    Tr_K = np.asarray(Tr_K)    
    (vza,Ta_K,u,ea,p,Sn_C,Sn_S,Lsky,LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
            leaf_width,z0_soil,alpha_PT,x_LAD,f_c,f_g,wc,calcG_array)=map(_CheckDefaultParameterSize,[
            vza,Ta_K,u,ea,p,Sn_C,Sn_S,Lsky,LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
            leaf_width,z0_soil,alpha_PT,x_LAD,f_c,f_g,wc,calcG_params[1]], [Tr_K]*24)
    # Create the output variables
    [flag, Ts, Tc, T_AC, Ln_S, Ln_C, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a]=[np.zeros(Tr_K.shape) for i in range(14)]
     
    # iteration of the Monin-Obukhov length
    if type(UseL)==bool:
        L = np.asarray(np.zeros(Ts.shape)+np.inf)     # Initially assume stable atmospheric conditions and set variables for 
        max_iterations=ITERATIONS
    else: # We force Monin-Obukhov lenght to the provided array/value
        L=np.asarray(np.ones(Ts.shape)*UseL)
        max_iterations=1 # No iteration
    # Calculate the general parameters
    rho= met.CalcRho(p, ea, Ta_K)  # Air density
    c_p = met.CalcC_p(p, ea)  # Heat capacity of air
    z_0H=res.CalcZ_0H(z_0M, kB=kB) # Roughness length for heat transport
    
    # Calculate LAI dependent parameters for dataset where LAI > 0
    omega0=CI.CalcOmega0_Kustas(LAI, f_c,x_LAD=x_LAD, isLAIeff=True)
    F = np.asarray(LAI/f_c) # Real LAI    
    f_theta = CalcFthetaCampbell(vza, F, wc=wc, Omega0=omega0,x_LAD=x_LAD)   # Fraction of vegetation observed by the sensor

    # Initially assume stable atmospheric conditions and set variables for 
    # iteration of the Monin-Obukhov length
    u_friction = MO.CalcU_star(u, zu, L, d_0, z_0M)
    u_friction = np.asarray(np.maximum(u_friction_min, u_friction))
    L_old = np.ones(Tr_K.shape)
    L_diff = np.asarray(np.ones(Tr_K.shape)*float('inf'))

    # First assume that canopy temperature equals the minumum of Air or radiometric T
    Tc = np.asarray(np.minimum(Tr_K, Ta_K))
    flag,Ts = CalcT_S(Tr_K, Tc, f_theta)

    # Outer loop for estimating stability. 
    # Stops when difference in consecutives L is below a given threshold
    for n_iterations in range(max_iterations):
        flag[np.isnan(Tr_K)]=255
        i = flag!=255
        if np.all(L_diff[i] < L_thres): 
            print("Finished interation with a max. L diff: "+str(np.max(L_diff)))
            break
        print("Iteration "+str(n_iterations)+", max. L diff: "+str(np.max(L_diff)))
         
        # Inner loop to iterativelly reduce alpha_PT in case latent heat flux 
        # from the soil is negative. The initial assumption is of potential 
        # canopy transpiration.
        flag[np.logical_and(L_diff >= L_thres, flag!=255)] = 0
        LE_S[np.logical_and(L_diff >= L_thres, flag!=255)] = -1
        alpha_PT_rec = np.asarray(alpha_PT + 0.1)         
        while np.any(LE_S[i] < 0): 
            i = np.logical_and.reduce((LE_S < 0, L_diff >= L_thres, flag != 255))            
            
            alpha_PT_rec[i] -= 0.1 
            
            # There cannot be negative transpiration from the vegetation 
            alpha_PT_rec[alpha_PT_rec <= 0.0] = 0.0 
            flag[np.logical_and(i,alpha_PT_rec == 0.0)] = 5 
             
            flag[np.logical_and.reduce((i, alpha_PT_rec < alpha_PT, alpha_PT_rec > 0.0))] = 3
            
            # Calculate the aerodynamic resistance
            R_a[i] = res.CalcR_A(zt[i], u_friction[i], L[i], d_0[i], z_0H[i])
            # Calculate soil and canopy resistances
            U_C = MO.CalcU_C_star (u_friction[i], hc[i], d_0[i], z_0M[i], L=L[i])
            if Resistance_flag == 0:
                u_S=MO.CalcU_Goudriaan (U_C, hc[i], omega0[i]*F[i], leaf_width[i], z0_soil[i]) # Clumped vegetation enhanced wind speed for the soil surface 
                u_d_zm=MO.CalcU_Goudriaan (U_C, hc[i], F[i], leaf_width[i], d_0[i]+z_0M[i]) # Wind speed is highly attenuated within the canopy volume
                R_x[i]=res.CalcR_X_Norman(LAI[i], leaf_width[i], u_d_zm) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_Kustas(u_S, Ts[i]-Tc[i])
            elif Resistance_flag ==1:
                R_x[i]=res.CalcR_X_Choudhury(U_C, LAI[i],leaf_width[i]) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_Choudhury(u_friction[i],hc[i],z_0M[i],d_0[i],zu[i],z0_soil[i])
            elif Resistance_flag ==2:
                R_x[i]=res.CalcR_X_McNaughton(LAI[i], leaf_width[i], u_friction[i]) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_McNaughton(u_friction[i])
            elif Resistance_flag ==3:
                alpha_k = MO.CalcA_Goudriaan(hc[i], omega0[i]*F[i], leaf_width[i]) # Clumped vegetation enhanced wind speed for the soil surface 
                alpha_prime = MO.CalcA_Goudriaan(hc[i], F[i], leaf_width[i]) # Wind speed is highly attenuated within the canopy volume
                R_x[i]=res.CalcR_X_Choudhury(U_C, LAI[i],leaf_width[i],alpha_prime=alpha_prime) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_Choudhury(u_friction[i],hc[i],z_0M[i],d_0[i],zu,z0_soil[i],alpha_k=alpha_k)
            else:
                u_S=MO.CalcU_Goudriaan (U_C, hc[i], omega0[i]*F[i], leaf_width[i], z0_soil[i]) # Clumped vegetation enhanced wind speed for the soil surface 
                u_d_zm=MO.CalcU_Goudriaan (U_C, hc[i], F[i], leaf_width[i], d_0[i]+z_0M[i]) # Wind speed is highly attenuated within the canopy volume
                R_x[i]=res.CalcR_X_Norman(LAI[i], leaf_width[i], u_d_zm) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_Kustas(u_S, Ts[i]-Tc[i])
            R_s=np.asarray(np.maximum( 0.001,R_s)) #CHANGE------------!!!!!!!!
            R_x=np.asarray(np.maximum( 0.001,R_x)) #CHANGE------------!!!!!!!!
            R_a=np.asarray(np.maximum( 0.001,R_a)) #CHANGE------------!!!!!!!!
            
            R_s=np.asarray(np.minimum( 200,R_s)) #CHANGE------------!!!!!!!!
            R_x=np.asarray(np.minimum( 80,R_x)) #CHANGE------------!!!!!!!!
            R_a=np.asarray(np.minimum( 80,R_a)) #CHANGE-----------!!!!!!!!!

            # Calculate net longwave radiation with current values of Tc and Ts
            Ln_C[i], Ln_S[i] = rad.CalcLnKustas (Tc[i], Ts[i], Lsky[i], LAI[i], emisVeg[i], emisGrd[i])
            delta_Rn = Sn_C + Ln_C 
            Rn_soil=Sn_S+Ln_S
            
            # Calculate the canopy and soil temperatures using the Priestley Taylor appoach
            H_C[i] = CalcH_C_PT(delta_Rn[i], f_g[i], Ta_K[i], p[i], c_p[i], alpha_PT_rec[i])
            Tc[i] = CalcT_C_Series(Tr_K[i], Ta_K[i], R_a[i], R_x[i], R_s[i], f_theta[i], H_C[i], rho[i], c_p[i])
            
            # Calculate soil temperature
            flag_t = np.zeros(flag.shape)            
            flag_t[i], Ts[i] = CalcT_S(Tr_K[i], Tc[i], f_theta[i])
            flag[flag_t==255] = 255
            LE_S[flag_t==255] = 0
            
            
            # Recalculate soil resistance using new soil temperature
            if Resistance_flag == 0:
                R_s[i] = res.CalcR_S_Kustas(u_S, Ts[i]-Tc[i])
                R_s = np.asarray(np.maximum( 0.001,R_s)) #CHANGE------------!!!!!!!!
                R_s=np.asarray(np.minimum( 200,R_s))   #CHANGE------------!!!!!!!!
            
            i = np.logical_and.reduce((LE_S < 0, L_diff >= L_thres, flag != 255))            
            
            # Get air temperature at canopy interface
            T_AC[i] = (( Ta_K[i]/R_a[i] + Ts[i]/R_s[i] + Tc[i]/R_x[i] )
                /(1.0/R_a[i] + 1.0/R_s[i] + 1.0/R_x[i]))
            
            # Calculate soil fluxes
            H_S[i] =  rho[i] * c_p[i] * (Ts[i] - T_AC[i])/ R_s[i]
            
            #Compute Soil Heat Flux Ratio
            G[i] = CalcG([calcG_params[0], calcG_array], Rn_soil, i)

            # Estimate latent heat fluxes as residual of energy balance at the
            # soil and the canopy            
            LE_S[i] = Rn_soil[i] - G[i] - H_S[i]
            LE_C[i] = delta_Rn[i] - H_C[i]        

            # Special case if there is no transpiration from vegetation. 
            # In that case, there should also be no evaporation from the soil
            # and the energy at the soil should be conserved.
            # See end of appendix A1 in Guzinski et al. (2015).  
            noT = np.logical_and(i, LE_C == 0)                    
            H_S[noT] = np.minimum(H_S[noT], Rn_soil[noT] - G[noT])
            G[noT] = np.maximum(G[noT], Rn_soil[noT] - H_S[noT])
            LE_S[noT] = 0

            # Calculate total fluxes
            H = np.asarray(H_C + H_S)
            LE = np.asarray(LE_C + LE_S)
            # Now L can be recalculated and the difference between iterations derived
            if type(UseL)==bool:        
                L[i] = MO.CalcL (u_friction[i], Ta_K[i], rho[i], c_p[i], H[i], LE[i])
                # Calculate again the friction velocity with the new stability correctios        
                u_friction[i] = MO.CalcU_star(u[i], zu[i], L[i], d_0[i], z_0M[i])
                u_friction = np.asarray(np.maximum(u_friction_min, u_friction))
                
            
    
        if type(UseL)==bool:        
            L_diff = np.asarray(np.fabs(L-L_old)/np.fabs(L_old))
            L_diff[np.isnan(L_diff)] = float('inf')
            L_old=np.array(L)
            L_old[L_old==0] = 1e-36 

    (flag, Ts, Tc, T_AC,L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,u_friction, 
         L,n_iterations)=map(np.asarray,(flag, Ts, Tc, T_AC,Ln_S, Ln_C, LE_C,H_C,
                LE_S,H_S,G,R_s,R_x,R_a,u_friction, L,n_iterations))    
    
    
    
    
    
#    print("HC is: ")
#    print(H_C[0,4])
#    print("R_a is: ")
#    print(R_a[0,4])
#    print("R_s is: ")
#    print(R_s[0,4])
#    print("R_x is: ")
#    print(R_x[0,4]) 
#    print("Tc is: ")
#    print(Tc[0,4])
#    print("Ta is: ")
#    print(Ta_K[0,4]) 
#    print("Ts is: ")
#    print(Ts[0,4]) 
#    
    return flag, Ts, Tc, T_AC,L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,u_friction, L,n_iterations
    
def  DTD(Tr_K_0,Tr_K_1,vza,Ta_K_0,Ta_K_1,u,ea,p,Sn_C,Sn_S,
             Lsky,LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,
             leaf_width=0.1,z0_soil=0.01,alpha_PT=1.26,x_LAD=1,f_c=1.0,f_g=1.0,wc=1.0,
             Resistance_flag=0,calcG_params=[[1],0.35], UseRi=False):
    ''' Calculate daytime Dual Time Difference TSEB fluxes
    
    Parameters
    ----------
    Tr_K_0 : float
        Radiometric composite temperature around sunrise(Kelvin).
    Tr_K_1 : float
        Radiometric composite temperature near noon (Kelvin).
    vza : float
        View Zenith Angle near noon (degrees).
    Ta_K_0 : float 
        Air temperature around sunrise (Kelvin).
    Ta_K_1 : float 
        Air temperature near noon (Kelvin).
    u : float 
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    Sn_C : float
        Canopy net shortwave radiation (W m-2).
    Sn_S : float
        Soil net shortwave radiation (W m-2).
    Lsky : float
        Downwelling longwave radiation (W m-2).
    LAI : float
        Effective Leaf Area Index (m2 m-2).
    hc : float
        Canopy height (m).
    emisVeg : float
        Leaf emissivity.
    emisGrd : flaot
        Soil emissivity.
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    zu : float
        Height of measurement of windspeed (m).
    zt : float
        Height of measurement of air temperature (m).
    leaf_width : Optional[float]
        average/effective leaf width (m).
    z0_soil : Optional[float]
        bare soil aerodynamic roughness length (m).
    alpha_PT : Optional[float]
        Priestley Taylor coeffient for canopy potential transpiration, 
        use 1.26 by default.
    x_LAD : Optional[float]
        Campbell 1990 leaf inclination distribution function chi parameter.
    f_c : Optiona;[float]
        Fractional cover.
    f_g : Optional[float]
        Fraction of vegetation that is green.
    wc : Optional[float]
        Canopy width to height ratio.
    Resistance_flag : int, optional
        Flag to determine which Resistances R_x, R_s model to use.
        
            * 0 [Default] Norman et al 1995 and Kustas et al 1999. 
            * 1 : Choudhury and Monteith 1988.
            * 2 : McNaughton and Van der Hurk 1995.   

    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
        
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_soil, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.CalcG_TimeDiff`).
    UseRi : float or None, optional
        If included, its value will be used to force the Richardson Number.
    
    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    Ts : float
        Soil temperature  (Kelvin).
    Tc : float
        Canopy temperature  (Kelvin).
    T_AC : float
        Air temperature at the canopy interface (Kelvin).
    L_nS : float
        Soil net longwave radiation (W m-2).
    L_nC : float
        Canopy net longwave radiation (W m-2).
    LE_C : float
        Canopy latent heat flux (W m-2).
    H_C : float
        Canopy sensible heat flux (W m-2).
    LE_S : float
        Soil latent heat flux (W m-2).
    H_S : float
        Soil sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    R_s : float
        Soil aerodynamic resistance to heat transport (s m-1).
    R_x : float
        Bulk canopy aerodynamic resistance to heat transport (s m-1).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    u_friction : float
        Friction velocity (m s-1).
    L : float
        Monin-Obuhkov length (m).
    Ri : float
        Richardson number.
    n_iterations : int
        number of iterations until convergence of L.

    References
    ----------
    .. [Norman2000] Norman, J. M., W. P. Kustas, J. H. Prueger, and G. R. Diak (2000),
        Surface flux estimation using radiometric temperature: A dual-temperature-difference
        method to minimize measurement errors, Water Resour. Res., 36(8), 2263-2274,
        http://dx.doi.org/10.1029/2000WR900033.
    .. [Guzinski2015] Guzinski, R., Nieto, H., Stisen, S., and Fensholt, R. (2015) Inter-comparison
        of energy balance and hydrological models for land surface energy flux estimation over
        a whole river catchment, Hydrol. Earth Syst. Sci., 19, 2017-2036,
        http://dx.doi.org/10.5194/hess-19-2017-2015.
    '''

    import numpy as np   
    
    # Convert input scalars to numpy arrays and parameters size
    Tr_K_0 = np.asarray(Tr_K_0)    
    (Tr_K_1,vza,Ta_K_0,Ta_K_1,u,ea,p,Sn_C,Sn_S,Lsky,LAI,hc,emisVeg,
         emisGrd,z_0M,d_0,zu,zt,leaf_width,z0_soil,alpha_PT,f_c,f_g,wc,
         calcG_array)=map(_CheckDefaultParameterSize,[Tr_K_1,vza,Ta_K_0,Ta_K_1,u,ea,p,Sn_C,
        Sn_S,Lsky,LAI,hc,emisVeg,emisGrd,z_0M,d_0,zu,zt,leaf_width,z0_soil,
        alpha_PT,f_c,f_g,wc,calcG_params[1]], [Tr_K_0]*25)

    # Create the output variables
    [flag, Ts, Tc, T_AC, Ln_S, Ln_C, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,H]=[np.zeros(Tr_K_1.shape) for i in range(15)]
    
    # Calculate the general parameters
    rho= met.CalcRho(p, ea, Ta_K_1)  # Air density
    c_p = met.CalcC_p(p, ea)  # Heat capacity of air 
    
    # Calculate LAI dependent parameters for dataset where LAI > 0
    omega0 = CI.CalcOmega0_Kustas(LAI, f_c, isLAIeff=True,x_LAD=x_LAD) # Clumping factor at nadir
    F = np.asarray(LAI/f_c) # Real LAI    
    f_theta = CalcFthetaCampbell(vza, F, wc=wc ,Omega0=omega0,x_LAD=x_LAD)   # Fraction of vegetation observed by the sensor

    # L is not used in the DTD, since Richardson number is used instead to
    # avoid dependance on non-differential temperatures. But it is still saved
    # in the output for testing purposes.    
    if type(UseRi)==bool:
        # Calculate the Richardson number    
        Ri = MO.CalcRichardson (u, zu, d_0, Tr_K_0, Tr_K_1, Ta_K_0, Ta_K_1)
    else: # We force Monin-Obukhov lenght to the provided array/value
        Ri=np.asarray(np.ones(Tr_K_1.shape)*UseRi)
        
    # calculate the resistances
    # First calcualte u_S, wind speed at the soil surface
    u_friction = MO.CalcU_star(u, zu, Ri, d_0, z_0M, useRi=True)
    u_friction = np.asarray(np.maximum(u_friction_min, u_friction))    
    z_0H=res.CalcZ_0H(z_0M,kB=kB)
    # First assume that canpy temperature equals the minumum of Air or radiometric T
    Tc = np.asarray(np.minimum(Tr_K_1, Ta_K_1))
    flag, Ts = CalcT_S(Tr_K_1, Tc, f_theta)
    
    # Outer loop until canopy and soil temperatures have stabilised 
    Tc_prev = np.zeros(Tr_K_1.shape)
    Tc_thres = 0.1
    Tc_diff = np.fabs(Tc - Tc_prev)
    for n_iterations in range(ITERATIONS):
        i = flag!=255
        if np.all(Tc_diff[i] < Tc_thres): 
            print("Finished interation with a max. Tc diff: "+str(np.max(Tc_diff)))
            break
        print("Iteration "+str(n_iterations)+", maximum Tc difference between iterations: "+str(np.max(Tc_diff))) 

        # Inner loop to iterativelly reduce alpha_PT in case latent heat flux 
        # from the soil is negative. The initial assumption is of potential 
        # canopy transpiration.
        flag[np.logical_and(Tc_diff >= Tc_thres, flag!=255)] = 0
        LE_S[np.logical_and(Tc_diff >= Tc_thres, flag!=255)] = -1
        alpha_PT_rec = np.asarray(alpha_PT + 0.1)         
        while np.any(LE_S < 0): 
            i = np.logical_and.reduce((LE_S < 0, Tc_diff >= Tc_thres, flag != 255))            
            
            alpha_PT_rec[i] -= 0.1 
            
            # There cannot be negative transpiration from the vegetation 
            alpha_PT_rec[alpha_PT_rec <= 0.0] = 0.0 
            flag[np.logical_and(i,alpha_PT_rec == 0.0)] = 5 
             
            flag[np.logical_and.reduce((i, alpha_PT_rec < alpha_PT, alpha_PT_rec > 0.0))] = 3

            # Calculate wind profiles
            U_C = MO.CalcU_C(u_friction[i], hc[i], d_0[i], z_0M[i])
            R_a[i] = res.CalcR_A (zu[i], u_friction[i], Ri[i], d_0[i], z_0H[i], useRi=True)
            if Resistance_flag == 0:
                u_S=MO.CalcU_Goudriaan (U_C, hc[i], omega0[i]*F[i], leaf_width[i], z0_soil[i]) # Clumped vegetation enhanced wind speed for the soil surface 
                u_d_zm=MO.CalcU_Goudriaan (U_C, hc[i], F[i], leaf_width[i], d_0[i]+z_0M[i]) # Wind speed is highly attenuated within the canopy volume
                R_x[i]=res.CalcR_X_Norman(LAI[i], leaf_width[i], u_d_zm) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                deltaT=np.asarray((Tr_K_1[i] - Tr_K_0[i]) - (Ta_K_1[i]- Ta_K_0[i]))#based on equation from Guzinski et. al., 2014 
                R_s[i]=res.CalcR_S_Kustas(u_S, deltaT)
            elif Resistance_flag ==1:
                R_x[i]=res.CalcR_X_Choudhury(U_C, LAI[i],leaf_width[i]) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_Choudhury(u_friction[i],hc[i],z_0M[i],d_0[i],zu[i],z0_soil[i])
            elif Resistance_flag ==2:
                R_x[i]=res.CalcR_X_McNaughton(LAI[i], leaf_width[i], u_friction[i]) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_McNaughton(u_friction[i])
            elif Resistance_flag ==3:
                alpha_k=MO.CalcA_Goudriaan(hc[i],omega0[i]*F[i],leaf_width[i]) # Clumped vegetation enhanced wind speed for the soil surface 
                alpha_prime=MO.CalcA_Goudriaan(hc[i],F[i],leaf_width[i])# Wind speed is highly attenuated within the canopy volume
                R_x[i]=res.CalcR_X_Choudhury(U_C, LAI[i],leaf_width[i],alpha_prime=alpha_prime) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                R_s[i]=res.CalcR_S_Choudhury(u_friction[i],hc[i],z_0M[i],d_0[i],zu,z0_soil[i],alpha_k=alpha_k)
            else:
                u_S=MO.CalcU_Goudriaan (U_C, hc[i], omega0[i]*F[i], leaf_width[i], z0_soil[i]) # Clumped vegetation enhanced wind speed for the soil surface 
                u_d_zm=MO.CalcU_Goudriaan (U_C, hc[i], F[i], leaf_width[i], d_0[i]+z_0M[i]) # Wind speed is highly attenuated within the canopy volume
                R_x[i]=res.CalcR_X_Norman(LAI[i], leaf_width[i], u_d_zm) # Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
                deltaT=np.asarray((Tr_K_1[i] - Tr_K_0[i]) - (Ta_K_1[i]- Ta_K_0[i]))#based on equation from Guzinski et. al., 2014 
                R_s[i]=res.CalcR_S_Kustas(u_S, deltaT)
            R_s=np.asarray(np.maximum( 1e-3,R_s))
            R_x=np.asarray(np.maximum( 1e-3,R_x))
            R_a=np.asarray(np.maximum( 1e-3,R_a))
   
            # Calculate net longwave radiation with current values of Tc and Ts
            Ln_C[i], Ln_S[i] = rad.CalcLnKustas (Tc[i], Ts[i], Lsky[i], LAI[i], emisVeg[i], emisGrd[i])
            
            # Calculate total net radiation of soil and canopy        
            delta_Rn = Sn_C + Ln_C 
            Rn_soil = Sn_S + Ln_S
            
            # Calculate sensible heat fluxes at time t1
            H_C[i] = CalcH_C_PT(delta_Rn[i], f_g[i], Ta_K_1[i], p[i], c_p[i], alpha_PT_rec[i])
            H[i] = CalcH_DTD_series(Tr_K_1[i], Tr_K_0[i], Ta_K_1[i], Ta_K_0[i], rho[i], c_p[i], f_theta[i],
                R_s[i], R_a[i], R_x[i], H_C[i])
            H_S[i] = H[i] - H_C[i]
                               
            # Calculate ground heat flux
            G[i] = CalcG([calcG_params[0], calcG_array], Rn_soil, i)
            
            # Calculate latent heat fluxes as residuals
            LE_S[i] = Rn_soil[i] - H_S[i] - G[i]
            LE_C[i] = delta_Rn[i] - H_C[i]

            # Special case if there is no transpiration from vegetation. 
            # In that case, there should also be no evaporation from the soil
            # and the energy at the soil should be conserved.
            # See end of appendix A1 in Guzinski et al. (2015).                      
            noT = np.logical_and(i, LE_C == 0)                    
            H_S[noT] = np.minimum(H_S[noT], Rn_soil[noT] - G[noT])
            G[noT] = np.maximum(G[noT], Rn_soil[noT] - H_S[noT])
            LE_S[noT] = 0
    
            # Recalculate soil and canopy temperatures. They are used only for  
            # estimation of longwave radiation, so the use of non-differential Tr
            # and Ta shouldn't affect the turbulent fluxes much
            Tc[i] = CalcT_C_Series(Tr_K_1[i], Ta_K_1[i], R_a[i], R_x[i], R_s[i], f_theta[i], H_C[i], rho[i], c_p[i])
            flag_t = np.zeros(flag.shape)            
            flag_t[i],Ts[i] = CalcT_S(Tr_K_1[i], Tc[i], f_theta[i])
            flag[flag_t==255] = 255
            LE_S[flag_t==255] = 0  
        
        Tc_diff = np.asarray(np.fabs(Tc - Tc_prev))
        Tc_prev= np.array(Tc)
    
    # L is only calculated for testing purposes
    L = MO.CalcL(u_friction, Ta_K_1, rho, c_p, H, LE_C + LE_S)

    (flag, Ts, Tc, T_AC,L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,u_friction, 
         L,Ri,n_iterations)=map(np.asarray,(flag, Ts, Tc, T_AC,Ln_S, Ln_C, LE_C,
                H_C,LE_S,H_S,G, R_s,R_x,R_a,u_friction, L,Ri,n_iterations))               
    return [flag, Ts, Tc, T_AC,L_nS,L_nC, LE_C,H_C,LE_S,H_S,G,R_s,R_x,R_a,
            u_friction, L,Ri,n_iterations]        

def  OSEB(Tr_K,Ta_K,u,ea,p,Sn,Lsky,emis,z_0M,d_0,zu,zt, calcG_params=[[1],0.35], UseL=False,T0_K = []):
    '''Calulates bulk fluxes from a One Source Energy Balance model

    Parameters
    ----------
    Tr_K : float
        Radiometric composite temperature (Kelvin).
    Ta_K : float 
        Air temperature (Kelvin).
    u : float 
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    S_n : float
        Solar irradiance (W m-2).
    Lsky : float
        Downwelling longwave radiation (W m-2)
    emis : float
        Surface emissivity.
    albedo : float
        Surface broadband albedo.        
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    zu : float
        Height of measurement of windspeed (m).
    zt : float
        Height of measurement of air temperature (m).
    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
        
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_soil, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.CalcG_TimeDiff`).
    UseL : Optional[float]
        If included, its value will be used to force the Moning-Obukhov stability length.
    T0_K: Optional[tuple(float,float)]
        If given it contains radiometric composite temperature (K) at time 0 as 
        the first element and air temperature (K) at time 0 as the second element, 
        in order to derive differential temperatures like is done in DTD
        
    
    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    Ln : float
        Net longwave radiation (W m-2)
    LE : float
        Latent heat flux (W m-2).
    H : float
        Sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    u_friction : float
        Friction velocity (m s-1).
    L : float
        Monin-Obuhkov length (m).
    n_iterations : int
        number of iterations until convergence of L.
    '''
    import numpy as np     
    
    # Convert input scalars to numpy arrays and check parameters size
    Tr_K = np.asarray(Tr_K)
    (Ta_K,u,ea,p,Sn,Lsky,emis,z_0M,d_0,zu,zt, calcG_array)=map(_CheckDefaultParameterSize,
        [Ta_K,u,ea,p,Sn,Lsky,emis,z_0M,d_0,zu,zt,calcG_params[1]], [Tr_K]*12)
    # Create the output variables
    [flag,Ln, LE,H,G,R_a]=[np.zeros(Tr_K.shape) for i in range(6)]
     
    # iteration of the Monin-Obukhov length
    if type(UseL)==bool:
        L = np.zeros(Tr_K.shape)+np.inf     # Initially assume stable atmospheric conditions and set variables for 
        max_iterations=ITERATIONS
    else: # We force Monin-Obukhov lenght to the provided array/value
        L=np.ones(Tr_K.shape)*UseL
        max_iterations=1 # No iteration

    # Check if differential temperatures are to be used  
    if len(T0_K) == 2:
        differentialT = True 
        Tr_K_0 = np.asarray(T0_K[0])
        Ta_K_0 = np.asarray(T0_K[1])
    else:
        differentialT = False
   
    # Initially assume stable atmospheric conditions and set variables for 
    L_old=np.ones(Tr_K.shape)
    # Calculate the general parameters
    rho= met.CalcRho(p, ea, Ta_K)  #Air density
    c_p = met.CalcC_p(p, ea)  #Heat capacity of air

    # With differential temperatures use Richardson number to approximate L,
    # same as is done in DTD    
    if differentialT:
        if type(UseL)==bool:
            Ri = MO.CalcRichardson (u, zu, d_0, Tr_K_0, Tr_K, Ta_K_0, Ta_K)
        else:
            Ri=np.array(L)
        u_friction = MO.CalcU_star(u, zu, Ri, d_0, z_0M, useRi=True) 
    else:
        u_friction = MO.CalcU_star(u, zu, L, d_0, z_0M)
    u_friction = np.maximum(u_friction_min, u_friction)
    L_old = np.ones(Tr_K.shape)
    L_diff = np.ones(Tr_K.shape)*float('inf')
    
    z_0H=res.CalcZ_0H(z_0M,kB=kB)
    
    # Calculate Net radiation
    Ln=emis*Lsky-emis*met.CalcStephanBoltzmann(Tr_K)
    Rn=np.asarray(Sn+Ln)

    #Compute Soil Heat Flux
    i=np.ones(Rn.shape,dtype=bool)
    G[i] = CalcG([calcG_params[0], calcG_array], Rn,i)
    
    # Loop for estimating atmospheric stability. 
    # Stops when difference in consecutive L and u_friction is below a 
    # given threshold
    for n_iterations in range(max_iterations):
        flag = np.zeros(Tr_K.shape)
        #Stop the iteration if differences are below the threshold
        if np.all(L_diff < L_thres):
            break

        # Calculate the aerodynamic resistances
        if differentialT:
            R_a=res.CalcR_A (zu, u_friction, Ri, d_0, z_0H, useRi=True)
        else:
            R_a=res.CalcR_A ( zt, u_friction, L, d_0, z_0H)
        R_a = np.asarray(np.maximum( 1e-3,R_a))
        
        # Calculate bulk fluxes assuming that since there is no vegetation,
        # Tr is the heat source
        if differentialT:
            H =  rho * c_p * ((Tr_K - Tr_K_0) - (Ta_K - Ta_K_0))/ R_a
        else:
            H =  rho * c_p * (Tr_K - Ta_K)/ R_a
        H=np.asarray(H)
        LE = np.asarray(Rn - G - H)
        
        # Avoid negative ET during daytime and make sure that energy is conserved
        flag[LE<0] = 5
        H[LE<0] = np.minimum(H[LE<0], Rn[LE<0] - G[LE<0])
        G[LE<0] = np.maximum(G[LE<0], Rn[LE<0] - H[LE<0])
        LE[LE<0] = 0
        
        if type(UseL)==bool:
            # Now L can be recalculated and the difference between iterations derived
            L=MO.CalcL (u_friction, Ta_K, rho, c_p, H, LE)
            L_diff=np.fabs(L-L_old)/np.fabs(L_old)
            L_old=np.array(L)
            L_old[np.fabs(L_old)==0] = 1e-36
    
            # Calculate again the friction velocity with the new stability correction
            # and derive the change between iterations
            if not differentialT:        
                u_friction=MO.CalcU_star (u, zu, L, d_0, z_0M)
                u_friction = np.maximum(u_friction_min, u_friction)
        
    flag,Ln, LE,H,G,R_a,u_friction, L,n_iterations=map(np.asarray,(flag,Ln, LE,
                        H,G,R_a,u_friction, L,n_iterations))
    
    return flag,Ln, LE,H,G,R_a,u_friction, L,n_iterations

def OSEB_BelowCanopy(Tr_K,Ta_K,u,ea,p,Sn,Lsky,LAI,hc,emisVeg, emisGrd,z_0M,d_0,zu,zt,
                       leaf_width=0.01,z0_soil=0.01,x_LAD=1.0,f_c=1.0,wc=1.0,
                       Resistance_flag=0,calcG_params=[[1],0.35],
                        Tc=None,UseL=None):
    '''Calculates a One Source Energy Balance below a canopy
    
    Parameters
    ----------
    Tr_K : float
        Radiometric composite temperature (Kelvin).
    Ta_K : float 
        Air temperature (Kelvin).
    u : float 
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    Sdn_dir : float
        Beam solar irradiance (W m-2).
    Sdn_dif : float
        Difuse solar irradiance (W m-2).
    fvis : float
        Fraction of shortwave radiation corresponding to the PAR region (400-700nm).
    fnir : float
        Fraction of shortwave radiation corresponding to the NIR region (700-2500nm).
    sza : float
        Solar Zenith Angle (degrees)
    Lsky : float
        Downwelling longwave radiation (W m-2)
    LAI : float
        Effective Leaf Area Index (m2 m-2).
    hc : float
        Canopy height (m).
    emisVeg : float
        Leaf emissivity.
    emisGrd : flaot
        Soil emissivity.
    spectraVeg : dict('rho_leaf_vis', 'tau_leaf_vis', 'rho_leaf_nir','tau_leaf_nir')
        Leaf spectrum dictionary.  
        
        rho_leaf_vis : float
            leaf bihemispherical reflectance in the visible (400-700 nm).
        tau_leaf_vis : float
            leaf bihemispherical transmittance in the visible (400-700nm).
        rho_leaf_nir : float
            leaf bihemispherical reflectance in the optical infrared (700-2500nm),
        tau_leaf_nir : float
            leaf bihemispherical transmittance in the optical  infrared (700-2500nm).
    spectraGrd : dict('rho rsoilv', 'rsoiln')
        Soil spectrum dictonary.
        
            rsoilv : float
                soil bihemispherical reflectance in the visible (400-700 nm).
            rsoiln : float
                soil bihemispherical reflectance in the optical infrared (700-2500nm).
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    zu : float
        Height of measurement of windspeed (m).
    zt : float
        Height of measurement of air temperature (m).
    leaf_width : Optional[float]
        average/effective leaf width (m).
    z0_soil : Optional[float]
        bare soil aerodynamic roughness length (m).
    x_LAD : Optional[float]
        Campbell 1990 leaf inclination distribution function chi parameter.
    f_c : Optiona;[float]
        Fractional cover.
    f_g : Optional[float]
        Fraction of vegetation that is green.
    wc : Optional[float]
        Canopy width to height ratio.
    Resistance_flag : int, optional
        Flag to determine which Resistance R_s model to use.
        
            * 0 [Default] Norman et al 1995 and Kustas et al 1999. 
            * 1 : Choudhury and Monteith 1988.
            * 2 : McNaughton and Van der Hurk 1995.   

    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
        
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_soil, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.CalcG_TimeDiff`).
    Tc : Optional(None,float)
        If provided, the temperature of the overstory.
    UseL : Optional[float]
        If included, its value will be used to force the Moning-Obukhov stability length.
    
    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    Ln : float
        Net longwave radiation of the understory (W m-2)
    LE : float
        Understory latent heat flux (W m-2).
    H : float
        Understory sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    R_s : float
        Soil aerodynamic resistance to heat transport (s m-1).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    u_friction : float
        Friction velocity (m s-1).
    L : float
        Monin-Obuhkov length (m).
    n_iterations : int
        number of iterations until convergence of L.
    '''
    
    import numpy as np
    
    # Convert the input scalars to numpy arrays and check parameters size
    Tr_K = np.asarray(Tr_K)    
    (Ta_K,u,ea,p,Sn,Lsky,LAI,hc,emisVeg, emisGrd,z_0M,d_0,zu,zt,
        leaf_width,z0_soil,x_LAD,f_c,wc,calcG_array)=map(_CheckDefaultParameterSize,
        [Ta_K,u,ea,p,Sn,Lsky,LAI,hc,emisVeg, emisGrd,z_0M,d_0,zu,zt,leaf_width,z0_soil,
        x_LAD,f_c,wc,calcG_params[1]], [Tr_K]*20)

    # Create the output variables
    flag, Ln, LE,H,G,R_s,R_a=[np.zeros(Tr_K.shape) for i in range(7)]
    # iteration of the Monin-Obukhov length
    if type(UseL)==bool:
        L = np.zeros(Tr_K.shape)+np.inf     # Initially assume stable atmospheric conditions and set variables for 
        max_iterations=ITERATIONS
    else: # We force Monin-Obukhov lenght to the provided array/value
        L=np.ones(Tr_K.shape)*UseL
        max_iterations=1 # No iteration

    u_friction = MO.CalcU_star(u, zu, L, d_0,z_0M)

    # Initially assume stable atmospheric conditions and set variables for 
    L_old = np.ones(Tr_K.shape)
    L_diff = np.ones(Tr_K.shape)*float('inf')

    # calculate the general parameters
    rho= met.CalcRho(p, ea, Ta_K)  #Air density
    c_p = met.CalcC_p(p, ea)  #Heat capacity of air  
    z_0H=res.CalcZ_0H(z_0M,kB=kB)

    # Net longwave radiation
    if not Tc:
        # We assume that the temperature of the canpy above equals the minumum of Air or radiometric
        Tc=np.array(Ta_K)
    _, Ln = rad.CalcLnKustas (Tc, Tr_K, Lsky, LAI,emisVeg, emisGrd,x_LAD)
    Rn=Sn+Ln
    #Compute Soil Heat Flux
    i=np.ones(Rn.shape,dtype=bool)
    G[i] = CalcG([calcG_params[0], calcG_array], Rn,i)
    # Calculate LAI dependent parameters for dataset where LAI > 0
    omega0 = CI.CalcOmega0_Kustas(LAI, f_c, isLAIeff=True,x_LAD=x_LAD) # Clumping factor at nadir
    F = np.asarray(LAI/f_c) # Real LAI    
    
    # Loop for estimating atmospheric stability. 
    # Stops when difference in consecutive L and u_friction is below a 
    # given threshold
    for n_iterations in range(max_iterations):
        #Stop the iteration if differences are below the threshold
        if np.all(L_diff < L_thres):
            break

        flag = np.zeros(Tr_K.shape)
        
        # calculate the aerodynamic resistances
        R_a=res.CalcR_A ( zt, u_friction, L, d_0, z_0H)
        # Calculate wind speed at the canopy height
        U_C=MO.CalcU_C_star (u_friction, hc, d_0, z_0M,L=L)
        # Calculate soil resistance
        if Resistance_flag == 0:
            u_S=MO.CalcU_Goudriaan (U_C, hc, F*omega0, leaf_width, z0_soil) # Clumped vegetation enhanced wind speed for the soil surface 
            R_s=res.CalcR_S_Kustas(u_S, Tr_K-Ta_K)
        elif Resistance_flag ==1:
            R_s=res.CalcR_S_Choudhury(u_friction,hc,z_0M,d_0,zu,z0_soil)
        elif Resistance_flag ==2:
            R_s=res.CalcR_S_McNaughton(u_friction)
        elif Resistance_flag ==3:
            alpha_k=MO.CalcA_Goudriaan(hc,F*omega0,leaf_width) # Clumped vegetation enhanced wind speed for the soil surface 
            R_s=res.CalcR_S_Choudhury(u_friction,hc,z_0M,d_0,zu,z0_soil,alpha_k=alpha_k)
        else:
            u_S=MO.CalcU_Goudriaan (U_C, hc, F*omega0, leaf_width, z0_soil) # Clumped vegetation enhanced wind speed for the soil surface 
            R_s=res.CalcR_S_Kustas(u_S, Tr_K-Ta_K)
        R_s=np.asarray(np.maximum( 1e-3,R_s))
        R_a=np.asarray(np.maximum( 1e-3,R_a))

        H =  np.asarray(rho * c_p * (Tr_K - Ta_K)/ (R_s+R_a))
        LE = np.asarray(Rn - G - H)
        
        # Avoid negative ET during daytime and make sure that energy is conserved
        flag[LE<0] = 5
        H[LE<0] = np.minimum(H[LE<0], Rn[LE<0] - G[LE<0])
        G[LE<0] = np.maximum(G[LE<0], Rn[LE<0] - H[LE<0])
        LE[LE<0] = 0
        
        if type(UseL)==bool:
            # Now L can be recalculated and the difference between iterations derived
            L=MO.CalcL (u_friction, Ta_K, rho, c_p, H, LE)
            L_diff=np.fabs(L-L_old)/np.fabs(L_old)
            L_old=np.array(L)
            L_old[np.fabs(L_old)==0] = 1e-36
    
            # Calculate again the friction velocity with the new stability correction
            # and derive the change between iterations
            u_friction=MO.CalcU_star (u, zu, L, d_0, z_0M)
            u_friction = np.maximum(u_friction_min, u_friction)
    
       
    return flag,Ln, LE,H,G,R_s,R_a,u_friction, L,n_iterations

def OSEB_Canopy(Tr_K,Ta_K,u,ea,p,Sn,Lsky,LAI,hc,emis,z_0M,d_0,zu,zt,
                       leaf_width=0.01,f_c=1.0,
                       Resistance_flag=0,calcG_params=[[0],0],UseL=None):
    '''Calculates a One Source Energy Balance of a dense a canopy
    
    Parameters
    ----------
    Tr_K : float
        Radiometric composite temperature (Kelvin).
    Ta_K : float 
        Air temperature (Kelvin).
    u : float 
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    S_n : float
        Net Shortwave radiation (W m-2).
    Lsky : float
        Downwelling longwave radiation (W m-2)
    LAI : float
        Effective Leaf Area Index (m2 m-2).
    hc : float
        Canopy height (m).
    emis : float
        Canopy emissivity.
    spectra : dict('rho_leaf_vis', 'tau_leaf_vis', 'rho_leaf_nir','tau_leaf_nir')
        Canopy spectrum dictionary.
        
            rho_leaf_vis : float
                leaf bihemispherical reflectance in the visible (400-700 nm).
            tau_leaf_vis : float
                leaf bihemispherical transmittance in the visible (400-700nm).
            rho_leaf_nir : float
                leaf bihemispherical reflectance in the optical infrared (700-2500nm),
            tau_leaf_nir : float
                leaf bihemispherical transmittance in the optical  infrared (700-2500nm).
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    zu : float
        Height of measurement of windspeed (m).
    zt : float
        Height of measurement of air temperature (m).
    leaf_width : Optional[float]
        average/effective leaf width (m).
    Resistance_flag : int, optional
        Flag to determine which Resistance R_x model to use.
        
            * 0 [Default] Norman et al 1995 and Kustas et al 1999. 
            * 1 : Choudhury and Monteith 1988.
            * 2 : McNaughton and Van der Hurk 1995.   

    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
        
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_soil, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.CalcG_TimeDiff`).
    UseL : Optional[float]
        If included, its value will be used to force the Moning-Obukhov stability length.
    
    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    L_n : float
        Net longwave radiation of the overstory (W m-2)
    LE : float
        Understory latent heat flux (W m-2).
    H : float
        Understory sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    R_x : float
        Canopy boundary layer resistance to heat transport (s m-1).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    u_friction : float
        Friction velocity (m s-1).
    L : float
        Monin-Obuhkov length (m).
    n_iterations : int
        number of iterations until convergence of L.
    '''  

    import numpy as np
    
    # Convert the input scalars to numpy arrays and check parameters size
    Tr_K = np.asarray(Tr_K)    
    (Ta_K,u,ea,p,Sn,Lsky,LAI,hc,emis,z_0M,d_0,zu,zt,leaf_width,f_c,
         calcG_array)=map(_CheckDefaultParameterSize,[Ta_K,u,ea,p,Sn,Lsky,LAI,hc,emis,z_0M,
        d_0,zu,zt,leaf_width,f_c,calcG_params[1]], [Tr_K]*16)
    # Create the output variables
    flag, Ln, LE,H,G,R_x,R_a=[np.zeros(Tr_K.shape) for i in range(7)]
    # iteration of the Monin-Obukhov length
    if type(UseL)==bool:
        L = np.zeros(Tr_K.shape)+np.inf     # Initially assume stable atmospheric conditions and set variables for 
        max_iterations=ITERATIONS
    else: # We force Monin-Obukhov lenght to the provided array/value
        L=np.ones(Tr_K.shape)*UseL
        max_iterations=1 # No iteration

    u_friction = MO.CalcU_star(u, zu, L, d_0,z_0M)

    # Initially assume stable atmospheric conditions and set variables for 
    L_old = np.ones(Tr_K.shape)
    L_diff = np.ones(Tr_K.shape)*float('inf')

    # calculate the general parameters
    F=LAI/f_c
    rho= met.CalcRho(p, ea, Ta_K)  #Air density
    c_p = met.CalcC_p(p, ea)  #Heat capacity of air  
    z_0H=res.CalcZ_0H(z_0M,kB=kB)
    # Net shortwave radiation
    Ln=emis*Lsky-emis*met.CalcStephanBoltzmann(Tr_K)
    Ln, _ = rad.CalcLnKustas (Tr_K, Tr_K, Lsky, LAI,emis, 0.95)
    Rn=Sn+Ln
    
    #Compute Soil Heat Flux
    i=np.ones(Rn.shape,dtype=bool)
    G[i] = CalcG([calcG_params[0], calcG_array], Rn,i)        
    # Calculate LAI dependent parameters for dataset where LAI > 0
    omega0 = CI.CalcOmega0_Kustas(LAI, f_c, isLAIeff=True,x_LAD=x_LAD) # Clumping factor at nadir
    F = np.asarray(LAI/f_c) # Real LAI    
        
    # loop for estimating stability, stop when difference in consecutives L is below 0.01
    for n_iterations in range(max_iterations):
        #Stop the iteration if differences are below the threshold
        if np.all(L_diff < L_thres):
            break

        flag = np.zeros(Tr_K.shape)
         
        # calculate the aerodynamic resistances
        R_a=res.CalcR_A ( zt, u_friction, L, d_0, z_0H)
        # Calculate wind speed at the canopy height
        U_C=MO.CalcU_C_star (u_friction, hc, d_0, z_0M,L=L)
        # Calculate soil resistance
        if Resistance_flag == 0:
            u_d_zm = MO.CalcU_Goudriaan (U_C, hc, F, leaf_width,d_0+z_0M) # Wind speed is highly attenuated within the canopy volume
            R_x=res.CalcR_X_Norman(LAI, leaf_width, u_d_zm)# Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
        elif Resistance_flag ==1:
            R_x=res.CalcR_X_Choudhury(U_C, LAI,leaf_width)# Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
        elif Resistance_flag ==2:
            R_x=res.CalcR_X_McNaughton(LAI, leaf_width, u_friction)# Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
        elif Resistance_flag ==3:
            alpha_prime=MO.CalcA_Goudriaan(hc,F,leaf_width) # Wind speed is highly attenuated within the canopy volume
            R_x=res.CalcR_X_Choudhury(U_C, LAI,leaf_width,alpha_prime=alpha_prime)# Vegetation in series with soil, i.e. well mixed, so we use the landscape LAI
        else:
            u_d_zm = MO.CalcU_Goudriaan (U_C, hc, F, leaf_width,d_0+z_0M) # Wind speed is highly attenuated within the canopy volume
            R_x=res.CalcR_X_Norman(LAI, leaf_width, u_d_zm)
        R_x=np.asarray(np.maximum( 1e-3,R_x))
        R_a=np.asarray(np.maximum( 1e-3,R_a))
        # calculate surface fluxes
        H =  np.asarray(rho * c_p * (Tr_K - Ta_K)/ (R_x+R_a))
        LE = np.asarray(Rn - G - H)
        # Avoid negative ET during daytime and make sure that energy is conserved
        flag[LE<0] = 5
        H[LE<0] = np.minimum(H[LE<0], Rn[LE<0] - G[LE<0])
        G[LE<0] = np.maximum(G[LE<0], Rn[LE<0] - H[LE<0])
        LE[LE<0] = 0
        
        if type(UseL)==bool:
            # Now L can be recalculated and the difference between iterations derived
            L=MO.CalcL (u_friction, Ta_K, rho, c_p, H, LE)
            L_diff=np.fabs(L-L_old)/np.fabs(L_old)
            L_old=np.array(L)
            L_old[np.fabs(L_old)==0] = 1e-36
    
            # Calculate again the friction velocity with the new stability correction
            # and derive the change between iterations
            u_friction=MO.CalcU_star (u, zu, L, d_0, z_0M)
            u_friction = np.maximum(u_friction_min, u_friction)
    
    (flag, Ln, LE,H,G,R_x,R_a,u_friction, L,n_iterations)=map(np.asarray,(flag, 
            Ln, LE,H,G,R_x,R_a,u_friction, L,n_iterations))
    
    return flag, Ln, LE,H,G,R_x,R_a,u_friction, L,n_iterations
  
def CalcFthetaCampbell(theta,F,wc=1,Omega0=1, x_LAD=1):
    '''Calculates the fraction of vegetatinon observed at an angle.
    
    Parameters
    ----------
    theta : float
        Angle of incidence (degrees).
    F : float
        Real Leaf (Plant) Area Index.
    wc : float
        Ratio of vegetation height versus width, optional (default = 1).
    Omega0 : float
        Clumping index at nadir, optional (default =1).
    x_LAD : float
        Chi parameter for the ellipsoidal Leaf Angle Distribution function, 
        use x_LAD=1 for a spherical LAD.
    
    Returns
    -------
    f_theta : float
        fraction of vegetation obsserved at an angle.
    
    References
    ----------
    .. [Campbell1998] Campbell, G. S. & Norman, J. M. (1998), An introduction to environmental
        biophysics. Springer, New York
        https://archive.org/details/AnIntroductionToEnvironmentalBiophysics.
    .. [Norman1995] J.M. Norman, W.P. Kustas, K.S. Humes, Source approach for estimating
        soil and vegetation energy fluxes in observations of directional radiometric
        surface temperature, Agricultural and Forest Meteorology, Volume 77, Issues 3-4,
        Pages 263-293, http://dx.doi.org/10.1016/0168-1923(95)02265-Y.
    '''

    import numpy as np    
    
    # First calcualte the angular clumping factor Omega based on eq (3) from
    # W.P. Kustas, J.M. Norman,  Agricultural and Forest Meteorology 94 (1999)
    OmegaTheta = Omega0 / (Omega0 + (1.0 - Omega0) * np.exp(-2.2 * np.radians(theta)**(3.8 - 0.46 * wc)))    #CHECK: should theta here be in degrees or radians
    # Estimate the beam extinction coefficient based on a elipsoidal LAD function
    # Eq. 15.4 of Campbell and Norman (1998)
    K_be=rad.CalcKbe_Campbell(theta,x_LAD)
    ftheta=1.0-np.exp(-K_be*OmegaTheta*F)
    return np.asarray(ftheta)

def CalcG(calcG_params, Rn_soil, i = None):
    import numpy as np    
    if i is None:
        i = np.ones(Rn_soil.shape, dtype=bool)
    if calcG_params[0][0]==0:
        G=calcG_params[1][i]
    elif calcG_params[0][0]==1:
        G=CalcG_Ratio(Rn_soil[i], calcG_params[1][i])
    elif calcG_params[0][0]==2:
        G=CalcG_TimeDiff (Rn_soil[i], [calcG_params[1][i], calcG_params[0][1], calcG_params[0][2], calcG_params[0][3]])
    return np.asarray(G)  

def CalcG_TimeDiff (R_n, G_param=[12.0,0.35, 3.0,24.0]):
    ''' Estimates Soil Heat Flux as function of time and net radiation.
    
    Parameters
    ----------
    R_n : float
        Net radiation (W m-2).
    G_param : tuple(float,float,float,float)
        tuple with parameters required (time, Amplitude,phase_shift,shape).
        
            time: float 
                time of interest (decimal hours).
            Amplitude : float 
                maximum value of G/Rn, amplitude, default=0.35.
            phase_shift : float
                shift of peak G relative to solar noon (default 3hrs after noon).
            shape : float
                shape of G/Rn, default 24 hrs.
    
    Returns
    -------
    G : float
        Soil heat flux (W m-2).

    References
    ----------
    .. [Santanello2003] Joseph A. Santanello Jr. and Mark A. Friedl, 2003: Diurnal Covariation in
        Soil Heat Flux and Net Radiation. J. Appl. Meteor., 42, 851-862,
        http://dx.doi.org/10.1175/1520-0450(2003)042<0851:DCISHF>2.0.CO;2.'''
    
    import numpy as np
    # Get parameters
    time=12.0-G_param[0]
    A = G_param[1]
    phase_shift=G_param[2]
    B = G_param[3]
    G_ratio=A*np.cos(2.0*np.pi*(time+phase_shift)/B)
    G = R_n * G_ratio
    return np.asarray(G)

def CalcG_Ratio(Rn_soil,G_ratio=0.35):
    '''Estimates Soil Heat Flux as ratio of net soil radiation.
    
    Parameters
    ----------
    Rn_soil : float
        Net soil radiation (W m-2).
    G_ratio : float, optional
        G/Rn_soil ratio, default=0.35.
    
    Returns
    -------
    G : float
        Soil heat flux (W m-2).

    References
    ----------
    .. [Choudhury1987] B.J. Choudhury, S.B. Idso, R.J. Reginato, Analysis of an empirical model
        for soil heat flux under a growing wheat crop for estimating evaporation by an
        infrared-temperature based energy balance equation, Agricultural and Forest Meteorology,
        Volume 39, Issue 4, 1987, Pages 283-297,
        http://dx.doi.org/10.1016/0168-1923(87)90021-9.
    '''
    import numpy as np
    G= G_ratio*Rn_soil
    return np.asarray(G)

def CalcH_C (T_C, T_A, R_A, rho, c_p):
    '''Calculates canopy sensible heat flux in a parallel resistance network.
    
    Parameters
    ----------
    T_C : float
        Canopy temperature (K).
    T_A : float
        Air temperature (K).
    R_A : float
        Aerodynamic resistance to heat transport (s m-1).
    rho : float
        air density (kg m-3).
    c_p : float
        Heat capacity of air at constant pressure (J kg-1 K-1).
    
    Returns
    -------
    H_C : float
        Canopy sensible heat flux (W m-2).'''
    import numpy as np
    H_C = rho*c_p*(T_C-T_A)/R_A
    return np.asarray(H_C)

def  CalcH_C_PT (delta_R_ni, f_g, T_a_K, P, c_p, alpha):
    '''Calculates canopy sensible heat flux based on the Priestley and Taylor formula.
    
    Parameters
    ----------
    delta_R_ni : float
        net radiation divergence of the vegetative canopy (W m-2).
    f_g : float
        fraction of vegetative canopy that is green.
    T_a_K : float
        air temperature (Kelvin).
    P : float
        air pressure (mb).
    c_p : float
        heat capacity of moist air (J kg-1 K-1).
    alpha : float 
        the Priestley Taylor parameter.
    
    Returns
    -------
    H_C : float
        Canopy sensible heat flux (W m-2).

    References
    ----------
    Equation 14 in [Norman1995]_
    '''  
    import numpy as np
    # slope of the saturation pressure curve (kPa./deg C)
    s = met.CalcDeltaVaporPressure( T_a_K)
    s=s*10 # to mb
    # latent heat of vaporisation (MJ./kg)
    Lambda=met.CalcLambda(T_a_K)
    # psychrometric constant (mb C-1)
    gama=met.CalcPsicr(P,Lambda)
    s_gama = s / (s + gama)
    H_C = delta_R_ni * (1.0 - alpha * f_g * s_gama)
    return np.asarray(H_C)

def CalcH_DTD_parallel (T_R1, T_R0, T_A1, T_A0, rho, c_p, f_theta1, R_S1, R_A1, R_AC1, H_C1):
    '''Calculates the DTD total sensible heat flux at time 1 with resistances in parallel.
    
    Parameters
    ----------
    T_R1 : float
        radiometric surface temperature at time t1 (K).
    T_R0 : float
        radiometric surface temperature at time t0 (K).
    T_A1 : float
        air temperature at time t1 (K).
    T_A0 : float
        air temperature at time t0 (K).
    rho : float
        air density at time t1 (kg m-3).
    cp : float
        heat capacity of moist air (J kg-1 K-1).
    f_theta_1 : float
        fraction of radiometer field of view that is occupied by vegetative cover at time t1.
    R_S1 : float
        resistance to heat transport from the soil surface at time t1 (s m-1).
    R_A1 : float
        resistance to heat transport in the surface layer at time t1 (s m-1).
    R_A1 : float
        resistance to heat transport at the canopy interface at time t1 (s m-1).
    H_C1 : float
        canopy sensible heat flux at time t1 (W m-2).
    
    Returns
    -------
    H : float
        Total sensible heat flux at time t1 (W m-2).

    References
    ----------
    .. [Guzinski2013] Guzinski, R., Anderson, M. C., Kustas, W. P., Nieto, H., and Sandholt, I. (2013)
        Using a thermal-based two source energy balance model with time-differencing to
        estimate surface energy fluxes with day-night MODIS observations,
        Hydrol. Earth Syst. Sci., 17, 2809-2825,
        http://dx.doi.org/10.5194/hess-17-2809-2013.
    '''
    import numpy as np

    #% Ignore night fluxes
    H = (rho*c_p *(((T_R1-T_R0)-(T_A1-T_A0))/((1.0-f_theta1)*(R_A1+R_S1))) +
        H_C1*(1.0-((f_theta1*R_AC1)/((1.0-f_theta1)*(R_A1+R_S1)))))
    return np.asarray(H)   
    
def CalcH_DTD_series(T_R1, T_R0, T_A1, T_A0, rho, c_p, f_theta, R_S, R_A, R_x, H_C):
    '''Calculates the DTD total sensible heat flux at time 1 with resistances in series
    
    Parameters
    ----------
    T_R1 : float
        radiometric surface temperature at time t1 (K).
    T_R0 : float
        radiometric surface temperature at time t0 (K).
    T_A1 : float
        air temperature at time t1 (K).
    T_A0 : float
        air temperature at time t0 (K).
    rho : float
        air density at time t1 (kg m-3).
    cp : float
        heat capacity of moist air (J kg-1 K-1).
    f_theta : float
        fraction of radiometer field of view that is occupied by vegetative cover at time t1.
    R_S : float
        resistance to heat transport from the soil surface at time t1 (s m-1).
    R_A : float
        resistance to heat transport in the surface layer at time t1 (s m-1).
    R_x : float
        Canopy boundary resistance to heat transport at time t1 (s m-1).
    H_C : float
        canopy sensible heat flux at time t1 (W m-2).
    
    Returns
    -------
    H : float
        Total sensible heat flux at time t1 (W m-2).

    References
    ----------
    .. [Guzinski2014] Guzinski, R., Nieto, H., Jensen, R., and Mendiguren, G. (2014)
        Remotely sensed land-surface energy fluxes at sub-field scale in heterogeneous
        agricultural landscape and coniferous plantation, Biogeosciences, 11, 5021-5046,
        http://dx.doi.org/10.5194/bg-11-5021-2014.
    '''
    import numpy as np
    H = rho*c_p*((T_R1-T_R0)-(T_A1-T_A0))/((1.0-f_theta)*R_S + R_A) + \
        H_C*((1.0-f_theta)*R_S - f_theta*R_x)/((1.0-f_theta)*R_S + R_A)
    return np.asarray(H)
     
def CalcH_S (T_S, T_A, R_A, R_S, rho, c_p):
    '''Calculates soil sensible heat flux in a parallel resistance network.
    
    Parameters
    ----------
    T_S : float
        Soil temperature (K).
    T_A : float
        Air temperature (K).
    R_A : float
        Aerodynamic resistance to heat transport (s m-1).
    R_A : float
        Aerodynamic resistance at the soil boundary layer (s m-1).
    rho : float
        air density (kg m-3).
    c_p : float
        Heat capacity of air at constant pressure (J kg-1 K-1).
   
    Returns
    -------
    H_C : float
        Canopy sensible heat flux (W m-2).

    References
    ----------
    Equation 7 in [Norman1995]_
    '''
    import numpy as np

    H_S = rho*c_p*((T_S-T_A)/(R_S+R_A))
    return np.asarray(H_S)
    
def  CalcT_C (T_R, T_S, f_theta):
    '''Estimates canopy temperature from the directional composite radiometric temperature.
    
    Parameters
    ----------
    T_R : float
        Directional Radiometric Temperature (K).
    T_S : float
        Soil Temperature (K).
    f_theta : float
        Fraction of vegetation observed.

    Returns
    -------
    flag : int
        Error flag if inversion not possible (255).
    T_C : float
        Canopy temperature (K).

    References
    ----------
    Eq. 1 in [Norman1995]_
    '''
    
    import numpy as np    
    # Convert input scalars to numpy array
    (T_R, T_S, f_theta)=map(np.asarray,(T_R, T_S, f_theta))
    T_temp = np.asarray(T_R**4 - (1.0 - f_theta)*T_S**4)
    T_C = np.zeros(T_R.shape)
    flag = np.zeros(T_R.shape)     
    
    # Succesfull inversion
    T_C[T_temp>=0] = ( T_temp[T_temp>=0] / f_theta[T_temp>=0])**0.25

    # Unsuccesfull inversion
    T_C[T_temp<0] = 1e-6
    flag[T_temp<0] = 255    

    return np.asarray(flag),np.asarray(T_C)


def CalcT_C_Series(Tr_K,Ta_K, R_a, R_x, R_s, f_theta, H_C, rho, c_p):
    '''Estimates canopy temperature from canopy sensible heat flux and 
    resistance network in series.
    
    Parameters
    ----------
    Tr_K : float
        Directional Radiometric Temperature (K).
    Ta_K : float
        Air Temperature (K).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    R_x : float
        Bulk aerodynamic resistance to heat transport at the canopy boundary layer (s m-1).
    R_s : float
        Aerodynamic resistance to heat transport at the soil boundary layer (s m-1).
    f_theta : float
        Fraction of vegetation observed.
    H_C : float
        Sensible heat flux of the canopy (W m-2).
    rho : float
        Density of air (km m-3).
    c_p : float
        Heat capacity of air at constant pressure (J kg-1 K-1).
    
    Returns
    -------
    T_c : float
        Canopy temperature (K).
    
    References
    ----------
    Eqs. A5-A13 in [Norman1995]_'''
    import numpy as np
    
    T_R_K_4=Tr_K**4
    # equation A7 from Norman 1995, linear approximation of temperature of the canopy
    T_C_lin = (( Ta_K/R_a + Tr_K/(R_s*(1.0-f_theta)) 
        + H_C*R_x/(rho*c_p)*(1.0/R_a + 1.0/R_s + 1.0/R_x)) 
        /(1.0/R_a + 1.0/R_s + f_theta/(R_s*(1.0 - f_theta))))
    # equation A12 from Norman 1995
    T_D = (T_C_lin*(1+R_s/R_a) - H_C*R_x/(rho*c_p)*(1.0 + R_s/R_x + R_s/R_a)
            - Ta_K*R_s/R_a)
    # equation A11 from Norman 1995
    delta_T_C = ((T_R_K_4 - f_theta*T_C_lin**4 - (1.0-f_theta)*T_D**4) 
        / (4.0* (1.0-f_theta)* T_D**3* (1.0 + R_s/R_a) + 4.0*f_theta*T_C_lin**3))
    # get canopy temperature in Kelvin
    Tc = T_C_lin + delta_T_C
    return np.asarray(Tc)
   
def CalcT_CS_Norman (F, vza_n, vza_f, T_n, T_f,wc=1,x_LAD=1, omega0=1):
    '''Estimates canopy and soil temperature by analytical inversion of Eq 1 in [Norman1995]
    of two directional radiometric observations. Ignoring shawows.
    
    Parameters
    ----------
    F : float
        Real Leaf (Plant) Area Index.
    vza_n : float
        View Zenith Angle during the nadir observation (degrees).
    vza_f : float
        View Zenith Angle during the oblique observation (degrees).
    T_n : float
        Radiometric temperature in the nadir obsevation (K).
    T_f : float
        Radiometric temperature in the oblique observation (K).
    wc : float,optional
        Canopy height to width ratio, use wc=1 by default.
    x_LAD : float,optional
        Chi parameter for the ellipsoildal Leaf Angle Distribution function of 
        Campbell 1988 [default=1, spherical LIDF].
    omega0 : float,optional
        Clumping index at nadir, use omega0=1 by default.
    
    Returns
    -------
    Tc : float
        Canopy temperature (K).
    Ts : float
        Soil temperature (K).
    
    References
    ----------
    inversion of Eq. 1 in [Norman1995]_
    '''
    
    import numpy as np
    # Convert  the input scalars to numpy arrays
    F, vza_n, vza_f, T_n, T_f,wc,x_LAD, omega0=map(np.asarray,
                    (F, vza_n, vza_f, T_n, T_f,wc,x_LAD, omega0))
    # Calculate the fraction of vegetation observed by each angle
    f_theta_n=CalcFthetaCampbell(vza_n, F, wc=wc,Omega0=omega0,x_LAD=x_LAD)
    f_theta_f=CalcFthetaCampbell(vza_f, F, wc=wc,Omega0=omega0,x_LAD=x_LAD)
    # Solve the sytem of two unknowns and two equations
    Ts_4=np.asarray((f_theta_f*T_n**4-f_theta_n*T_f**4)/(f_theta_f-f_theta_n))
    Tc_4=np.asarray((T_n**4-(1.0-f_theta_n)*Ts_4)/f_theta_n)
    
    Tc_K = np.zeros(T_n.shape)    
    Ts_K = np.zeros(T_n.shape)    
    
    # Successful inversion
    i = np.logical_and(Tc_4 > 0, Ts_4 > 0)
    Tc_K[i] = Tc_4[i]**0.25   
    Ts_K[i] = Ts_4[i]**0.25

    # Unsuccessful inversion
    Tc_K[~i] = float('nan')    
    Ts_K[~i] = float('nan')    
    
    return np.asarray(Tc_K), np.asarray(Ts_K)

def CalcT_CS_4SAIL (LAI, lidf, hotspot, Eo_n, Eo_f, L_sky, sza_n, sza_f, vza_n, vza_f, psi_n, psi_f, e_v,e_s):
    '''Estimates canopy and soil temperature by analytical inversion of 4SAIL (Eq. 12 in [Verhoef2007]_)
    of two directional radiometric observations. Ignoring shadows.
        
    Parameters
    ----------
    LAI : float
        Leaf (Plant) Area Index.
    lidf : list
        Campbell 1988 Leaf Inclination Distribution Function, default 5 degrees angle step.
    hotspot : float
        hotspot parameters, use 0 to ignore the hotspot effect (turbid medium).
    Eo_n : float
        Surface land Leaving thermal radiance (emitted thermal radiation).
        at the nadir observation (W m-2).
    Eo_f : float
        Surface land Leaving thermal radiance (emitted thermal radiation)
        at the oblique observation (W m-2).
    Lsky : float
        Broadband incoming longwave radiation (W m-2).
    sza_n : float
        Sun Zenith Angle during the nadir observation (degrees).
    sza_f : float
        Sun Zenith Angle during the oblique observation (degrees).
    vza_n : float
        View Zenith Angle during the nadir observation (degrees).
    vza_f : float
        View Zenith Angle during the oblique observation (degrees).
    psi_n : float
        Relative (sensor-sun) Azimuth Angle during the nadir observation (degrees).
    psi_f : float
        Relative (sensor-sun) Azimuth Angle during the oblique observation (degrees).
    e_v : float
        broadband leaf emissivity.
    e_s : float
        broadband soil emissivity.
    
    Returns
    -------
    Tc_K : float
        Canopy temperature (K).
    Ts_K : float
        Soil temperature (K).
    
    References
    ----------
    .. [Verhoef2007] Verhoef, W.; Jia, Li; Qing Xiao; Su, Z., (2007) Unified Optical-Thermal
        Four-Stream Radiative Transfer Theory for Homogeneous Vegetation Canopies,
        IEEE Transactions on Geoscience and Remote Sensing, vol.45, no.6, pp.1808-1822,
        http://dx.doi.org/10.1109/TGRS.2007.895844 based on  in Verhoef et al. (2007)
    '''
    import numpy as np

    # Apply Kirkchoff to work with reflectances instead of emissivities
    r_s=1.-e_s;
    r_v=1.-e_v;
    #Get nadir parameters for the inversion
    [rdot_star_n,emiss_v_eff_n,emiss_s_eff_n,gamma_sot,
         emiss_sot]=Get4SAILEmissionParam(LAI,hotspot,lidf,sza_n,vza_n,psi_n,r_v,r_s)
    # Calculate the total emission of the surface at nadir observation
    L_emiss_n=Eo_n-rdot_star_n*L_sky
    #Get forward parameters for the inversion
    [rdot_star_f,emiss_v_eff_f,emiss_s_eff_f,gamma_sot,
         emiss_sot]=Get4SAILEmissionParam(LAI,hotspot,lidf,sza_f,vza_f,psi_f,r_v,r_s)
    # Calculate the total emission of the surface at oblique observation
    L_emiss_f=Eo_f-rdot_star_f*L_sky
    # Invert 4SAIL to get the BB emission of vegetation and soil
    H_v=((emiss_s_eff_n*L_emiss_f-emiss_s_eff_f*L_emiss_n)/
        (emiss_s_eff_n*emiss_v_eff_f-emiss_s_eff_f*emiss_v_eff_n))
    H_s=(L_emiss_n-emiss_v_eff_n*H_v)/emiss_s_eff_n
    # Invert Stephan Boltzmann to obtain vegetation and soil temperatures
    Tc_K=(H_v/sb)**0.25
    Ts_K=(H_s/sb)**0.25
    return np.asarray(Tc_K), np.asarray(Ts_K)

def Get4SAILEmissionParam(LAI,hotspot,lidf,sza,vza,psi,rho_v,rho_s,tau_v=0.0):
    '''Calculates the effective surface reflectance, and emissivities for 
    soil and canopy using 4SAIL.
    
    Parameters
    ----------
    LAI : float
        Leaf (Plant) Area Index.
    hotspot : float
        hotspot parameters, use 0 to ignore the hotspot effect (turbid medium).
    lidf : list
        Campbell 1988 Leaf Inclination Distribution Function, 5 angle step.
    sza : float
        Sun Zenith Angle during the nadir observation (degrees).
    vza : float
        View Zenith Angle during the nadir observation (degrees).
    psi : float
        Relative (sensor-sun) Azimuth Angle during the nadir observation (degrees).
    psi_f : float
        Relative (sensor-sun) Azimuth Angle during the oblique observation (degrees).
    rho_v : float
        leaf reflectance (1-leaf emissivity).
    rho_s : float
        soil emissivity (1-soil emissivity).
    tau_v : float
        leaf transmittance (default zero transmittance in the TIR).
    
    Returns
    -------
    rdot_star : float
        surface effective reflectance.
    emiss_v_eff : float
        canopy effective emissivity.
    emiss_s_eff : float
        soil effective emissivity.
    gamma_sot : float
        directional canopy absortivity.
    emiss_sot : float
        directional canopy emissivity.

    References
    ----------
    Equations 5, 11, and 13 in [Verhoef2007]_
    '''
    import numpy as np
    from FourSAIL import FourSAIL
    # Run 4 SAIL
    [tss, too, tsstoo, rdd, tdd, rsd, tsd, rdo, tdo, rso, rsos, rsod, rddt, 
        rsdt, rdot, rsodt, rsost, rsot,gamma_sdf,gammas_db,gamma_so]=FourSAIL(LAI,
        hotspot,lidf,sza,vza,psi,rho_v,tau_v,rho_s)
    # Eq. 5 in [Verhoef2007]_    
    gamma_d=1.-rdd-tdd
    gamma_o=1.-rdo-tdo-too 
    # Eq. 13 in [Verhoef2007]_
    dn=1.0-rddt-rdd
    emiss_o=1.-rdot
    emiss_d=1.-rddt
    rdot_star=rdo+(tdd*(rddt*tdo+rdot*too)/dn) 
    # Get the coefficients from Eq 11 in [Verhoef2007]_
    emiss_v_eff=gamma_o+(gamma_d*(rddt*tdo+rdot*too)/dn) # 2nd element in Eq. 11 [Verhoef2007]_
    emiss_s_eff=emiss_o*too+(emiss_d*(tdo+rdd*rdot*too)/dn) # 3rd element in Eq. 11 [Verhoef2007]_
    gamma_sot=gamma_so+(gamma_sdf*(rddt*tdo+rdot*too)/dn) # 4th element in Eq. 11 [Verhoef2007]_
    emiss_sot=emiss_o*tsstoo+tss*(emiss_d*(tdo+rdd*rdot*too)/dn) # 5th element in Eq. 11 [Verhoef2007]_ 
    
    rdot_star,emiss_v_eff,emiss_s_eff,gamma_sot,emiss_sot=map(np.asarray,(rdot_star,
                            emiss_v_eff,emiss_s_eff,gamma_sot,emiss_sot))
    return rdot_star,emiss_v_eff,emiss_s_eff,gamma_sot,emiss_sot

def  CalcT_S (T_R, T_C, f_theta):
    '''Estimates soil temperature from the directional LST.
    
    Parameters
    ----------
    T_R : float
        Directional Radiometric Temperature (K).
    T_C : float
        Canopy Temperature (K).
    f_theta : float
        Fraction of vegetation observed.

    Returns
    -------
    flag : float
        Error flag if inversion not possible (255).
    T_S: float
        Soil temperature (K).
    
    References
    ----------
    Eq. 1 in [Norman1995]_'''
    
    
    import numpy as np    
    # Convert the input scalars to numpy arrays
    T_R, T_C, f_theta=map(np.asarray,(T_R, T_C, f_theta))
    T_temp = T_R**4 - f_theta*T_C**4
    T_S = np.zeros(T_R.shape)
    flag = np.zeros(T_R.shape)     
    
    # Succesfull inversion
    T_S[T_temp>=0] = ( T_temp[T_temp>=0] / (1.0 - f_theta[T_temp>=0]))**0.25

    # Unsuccesfull inversion
    T_S[T_temp<0] = 1e-6
    flag[T_temp<0] = 255

    return np.asarray(flag),np.asarray(T_S)

def CalcT_S_4SAIL (T_R, T_C, rdot_star,emiss_v_eff,emiss_s_eff,Lsky=0):
    '''Estimates canopy temperature from the directional LST using 4SAIL parameters.
    
    Parameters
    ----------
    T_R : float
        Directional Radiometric Temperature (K)
    T_S : float
        Soil Temperature (K)
    rdot_star : float
        surface effective reflectance
    emiss_v_eff : float
        canopy effective emissivity
    emiss_s_eff : float
        soil effective emissivity
    Lsky : float
        downwelling atmospheric longwave radiance (W m-2)

    Returns
    -------
    flag : int
        Error flag if inversion not possible (255).
    Ts : float
        Soil temperature (K).'''
    
    import numpy as np    
    Hv=met.CalcStephanBoltzmann(T_C)
    Eo=met.CalcStephanBoltzmann(T_R)
    Hs=np.asarray((Eo - rdot_star*Lsky-emiss_v_eff*Hv)/emiss_s_eff)
    
    flag = np.zeros(T_R.shape)     
    T_S= np.zeros(T_R.shape)
    
    # Succesfull inversion
    T_S[Hs>=0]=(Hs/sb)**0.25        
    
    # Unsuccesfull inversion
    T_S[Hs<0] = 1e-6 #Set this to -999?
    flag[Hs<0] = 255
    
    return np.asarray(flag),np.asarray(T_S)

def CalcT_S_Series(Tr_K,Ta_K,R_a,R_x,R_s,f_theta,H_S,rho,c_p):
    '''Estimates soil temperature from soil sensible heat flux and 
    resistance network in series.
    
    Parameters
    ----------
    Tr_K : float
        Directional Radiometric Temperature (K).
    Ta_K : float
        Air Temperature (K).
    R_a : float
        Aerodynamic resistance to heat transport (s m-1).
    R_x : float
        Bulk aerodynamic resistance to heat transport at the canopy boundary layer (s m-1).
    R_s : float
        Aerodynamic resistance to heat transport at the soil boundary layer (s m-1).
    f_theta : float
        Fraction of vegetation observed.
    H_S : float
        Sensible heat flux of the soil (W m-2).
    rho : float
        Density of air (km m-3).
    c_p : float
        Heat capacity of air at constant pressure (J kg-1 K-1).
        
    Returns
    -------
    T_s: float
        Soil temperature (K).
    T_c : float
        Air temperature at the canopy interface (K).
    
    References
    ----------
    Eqs. A15-A19 from [Norman1995]_'''
    import numpy as np
    #Eq. A.15 Norman 1995
    T_ac_lin=(((Ta_K/R_a)+(Tr_K/(f_theta*R_x))-
        (((1.0-f_theta)/(f_theta*R_x))*H_S*R_s/(rho*c_p))+H_S/(rho*c_p))/
        ((1.0/R_a)+(1.0/R_x)+(1.0-f_theta)/(f_theta*R_x)))    
    #Eq. A.17 Norman 1995
    T_e=T_ac_lin*(1.0+(R_x/R_a))-H_S*R_x/(rho*c_p)-Ta_K*R_x/R_a    
     #Eq. A.16 Norman 1995
    Delta_T_ac=((Tr_K**4-(1.0-f_theta)*(H_S*R_s/(rho*c_p)+T_ac_lin)**4-f_theta*T_e**4)/
        (4*f_theta*T_e**3.0*(1.0+(R_x/R_a))+4.0*(1.0-f_theta)*(H_S*R_s/(rho*c_p)+T_ac_lin)**3))
    #Eq. A.18 Norman 1995
    T_ac=T_ac_lin+Delta_T_ac    
    T_s=T_ac+H_S*R_s/(rho*c_p)
    return np.asarray(T_s),np.asarray(T_ac)
          
    
def _CheckDefaultParameterSize(parameter,input_array):
    import numpy as np
    parameter=np.asarray(parameter)
    if parameter.size==1:
        parameter=np.ones(input_array.shape)*parameter
        return np.asarray(parameter)
    elif parameter.shape!=input_array.shape:
        raise ValueError('dimension mismatch between parameter array and input array with shapes %s and %s'
            %(parameter.shape,input_array.shape))
    else:
        return np.asarray(parameter)
