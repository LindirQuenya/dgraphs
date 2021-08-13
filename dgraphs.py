#!/usr/bin/env python3
#from mpmath import mp

# TODO: Replace all eq:label citations with equation numbers.

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# First one uses script location, second uses cwd. Comment/uncomment to taste.
mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(__file__))
#mpl.rcParams["savefig.directory"] = ""

# The following specifies that titles should be on/off for figures.
TITLES=False

def print_help():
    print("""
Usage: python3 dgraphs.py <STUDY> <COMPARISON> [GRAPHTYPE]
Generate a graph that compares two variables for a given mission concept.

  -h, --help     display this message

  STUDY          one of [habex, luvoir, luvoiralt]
  COMPARISON     integer, [0..3] such that: (x vs y)
                          0  -> N_oplus    vs d
                          1  -> eta_oplus  vs d
                          2  -> d          vs N_oplus
                          3  -> duration   vs d
                          4  -> eta_oplus  vs N_oplus
                          5  -> IWA        vs N_oplus
                          6  -> N_oplus    vs SNR0
                          7  -> N_oplus    vs duration
                          8  -> eta_oplus  vs duration
                          9  -> efficiency vs N_oplus
                          10 -> iwafactor  vs N_oplus
  GRAPHTYPE      one of [linear, logy, logx, loglog]

If an invalid value is supplied for STUDY, it will
default to habex.

If an invalid value is supplied for COMPARISON,
the program will display this message and exit.

If GRAPHTYPE is not supplied or an invalid value
is supplied for GRAPHTYPE, the default will depend
on COMPARISON.
""")

if '-h' in sys.argv or '--help' in sys.argv:
    print_help()
    #exit(0)

HABEX=True
LUVOIR=False
LUVOIRALT=False
# Confusing, got renamed, deal with it.
if len(sys.argv) > 1:
    if 'luvoiralt' == sys.argv[1]:
        HABEX=False
        LUVOIR=True
        LUVOIRALT=False
    elif 'habex' == sys.argv[1]:
        HABEX=True
        LUVOIR=False
        LUVOIRALT=False
    elif 'luvoir' == sys.argv[1]:
        HABEX=False
        LUVOIR=False
        LUVOIRALT=True
    else:
        print_help()
        exit(1)


NOPLUSvD, ETAvD, DvNOPLUS, DvTIME, ETAvNOPLUS, IWAvNOPLUS, SNR0vNOPLUS, NOPLUSvT, ETAvT, IWAFACvNOPLUS, EPSvNOPLUS = [False]*11
if len(sys.argv) > 2:
    if '0' == sys.argv[2]:
        NOPLUSvD=True
    elif '1' == sys.argv[2]:
        ETAvD=True
    elif '2' == sys.argv[2]:
        DvNOPLUS=True
    elif '3' == sys.argv[2]:
        DvTIME=True
    elif '4' == sys.argv[2]:
        ETAvNOPLUS=True
    elif '5' == sys.argv[2]:
        IWAvNOPLUS=True
    elif '6' == sys.argv[2]:
        SNR0vNOPLUS=True
    elif '7' == sys.argv[2]:
        NOPLUSvT = True
    elif '8' == sys.argv[2]:
        ETAvT = True
    elif '9' == sys.argv[2]:
        EPSvNOPLUS=True
    elif '10' == sys.argv[2]:
        IWAFACvNOPLUS=True
    else:
        print_help()
        #exit(1)
else:
    print_help()
    #exit(1)

LINEAR, LOGY, LOGX, LOGLOG, UNSET = [False]*5
if len(sys.argv) > 3:
    if 'linear' == sys.argv[3]:
        LINEAR=True
    elif 'logy' == sys.argv[3]:
        LOGY=True
    elif 'logx' == sys.argv[3]:
        LOGX=True
    elif 'loglog' == sys.argv[3]:
        LOGLOG=True
    else:
        UNSET=True
        print_help()
else:
    UNSET=True


npf64 = np.float64
# Just shortening these, to make expressions more concise.
m = np.multiply; d = np.divide; a = np.add; s = np.subtract; e=np.e; pi=np.pi

h=npf64("6.62607015e-27")      # Planck's constant,            in erg-seconds
c=npf64('299792458')           # Speed of light,               in meters/second
k=npf64("1.380649e-16")        # Boltzmann constant,           in erg-kelvins
l=npf64("5e-7")                # Middle wavelength~=500nm,     in meters
temp=npf64('5778')             # Solar temperature=5778K,      in kelvins
r=npf64('695700e+3')           # Solar radius,                 in meters
avar=npf64('1.496e+11')        # Exo-earth orbital radius=1AU, in meters
etaearth=npf64('0.24')         # Eta-earth (Belikov 2017),     in unitless + (Kopparapu et al. 2018)
rhostar=npf64('4.76514807e-51')# Density of stars,             in meters^{-3} Assumed: 5% stars are
rhostar*=0.05                  # Singling out solar analogs,   in meters^{-3} solar analogues. Fair?
rval=1+npf64('8.75')           # Background counts +/-2.29,    in unitless

# For HabEx mission concept: Objective 1.
if HABEX:
    res=npf64('140')           # Spectral resolution,          in unitless
    dl=d(l,res)                # Wavelength band width,        in meters
    T=npf64('63113851.9')      # Survey duration=2 years,      in seconds
    K=npf64('1e-10')           # Flux contrast ratio,          in unitless
    eps=npf64('0.5')           # Observatory efficiency,       in unitless
    dvar=npf64('4.0')          # HabEx primary diameter,       in meters
    iwafac=npf64('3')          # iwa=3 l/d, iwafac=3,          in unitless
    Noplus=npf64('20')         # Target Exo-Earth yield,       in unitless No-prior
    snr0=npf64('7')            # Target signal-to-noise ratio, in unitless

# For LUVOIR mission concept: Sig. Sci. Case #1.
if LUVOIR:
    res=npf64('140')           # Spectral resolution,          in unitless
    dl=d(l,res)                # Wavelength band width,        in meters
    T=npf64('63113851.9')      # Survey duration=2 years,      in seconds
    K=npf64('1e-10')           # Flux contrast ratio,          in unitless
    eps=npf64('0.5')           # Observatory efficiency,       in unitless
    dvar=npf64('8.0')          # LUVOIR primary diameter,      in meters   Comment: inscribed d=6.7?
    iwafac=npf64('3')          # iwa=3 l/d, iwafac=3,          in unitless Comment: LUVOIR says 4l/d
    Noplus=npf64('28')         # Target HEC yield,             in unitless
    snr0=npf64('5')            # Target signal-to-noise ratio, in unitless

if LUVOIRALT:
    res=npf64('140')           # Spectral resolution,          in unitless
    dl=d(l,res)                # Wavelength band width,        in meters
    T=npf64('63113851.9')      # Survey duration=2 years,      in seconds
    K=npf64('1e-10')           # Flux contrast ratio,          in unitless
    eps=npf64('0.5')           # Observatory efficiency,       in unitless
    dvar=npf64('6.7')          # LUVOIR primary diameter,      in meters   Comment: inscribed d=6.7?
    iwafac=npf64('4')          # iwa=3 l/d, iwafac=3,          in unitless Comment: LUVOIR says 4l/d
    Noplus=npf64('28')         # Target HEC yield,             in unitless
    snr0=npf64('5')            # Target signal-to-noise ratio, in unitless

if not (LUVOIR^HABEX^LUVOIRALT):
    raise ValueError('Exactly one of HABEX, LUVOIR, and LUVOIRALT must be set!')

#Extra factor of 4, because surface area of sun, rather than cross-sectional area.
R=d(m(32*pi*pi, m(r, m(r, m(c, dl)))), m(np.power(l, 4), np.expm1(d(m(h, c), m(l, m(k, temp))))))

# Note that 'a' is appended to the end of all the names, to prevent name conflicts
# with the global vars.
# Eq Label: "eq:dphotonnoprioriwa"
def d_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    n = d(m(iwafaca, la), avara)
    cuberoot = np.cbrt(d(m(3, Noplusa), m(4*pi, m(rhostara, etaeartha))))
    is_low = m(5, m(Ra, m(Ka, m(epsa, m(np.power(n, 2), m(Ta, etaeartha))))))
    is_high = m(96, m(m(rvala,np.power(snr0a, 2)), Noplusa))
    innersqrt = np.sqrt(a(1, np.power(d(is_high, is_low), 2)))
    outersqrt = np.sqrt(a(1, innersqrt))
    return m(d(n, np.sqrt(2)), m(cuberoot, outersqrt))

# Eq Label: "eq:iwaApprox"
def iwa_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    return d(m(iwafaca,la),d_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala))

# Eq Label: "eq:dphotonnoise"
def d_photon_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    cuberoot = np.cbrt(d(3, m(4*pi, m(rhostara, etaeartha))))
    onlysqrt = np.sqrt(d(m(m(3, rvala), np.power(Noplusa, 5/3)), m(5, m(Ra, m(Ka, m(Ta, epsa))))))
    return m(m(4, snr0a), m(cuberoot, onlysqrt))

# Eq Label: "eq:dminiwa"
def d_iwa_all(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    cuberoot = np.cbrt(d(m(3, Noplusa), m(4*pi, m(rhostara, etaeartha))))
    if type(Ta)==np.ndarray:
        return np.array([m(d(m(iwafaca, la), avara), cuberoot) for x in Ta])
    return m(d(m(iwafaca, la), avara), cuberoot)

# Eq Label: "eq:iwaphotonint_noplus"
def intersection_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    return m(m(5/48, np.power(iwafaca, 2)), d(m(Ra, m(Ka, m(Ta, m(epsa, np.power(la, 2))))), m(rvala, np.power(m(snr0a, avara), 2))))

# No eq label, just a combo of the above.
def d_prior_noplus(Noplusarr, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    noplus_intersecta = np.array([intersection_prior([], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)])
    intersect_ind = int(np.where(np.isclose(Noplusarr,noplus_intersecta))[0])
    iwalim = d_iwa_all(Noplusarr[:intersect_ind], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    iwalim2 = d_iwa_all(Noplusarr[intersect_ind:], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    photonlim = d_photon_prior(Noplusarr[intersect_ind:], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    photonlim2 = d_photon_prior(Noplusarr[:intersect_ind], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    return np.concatenate((iwalim, photonlim)), np.concatenate((photonlim2, iwalim2))

# Trivial rearrangement of Eq Label : "eq:iwaphotonint_noplus" TODO: check if this needs an extra iwafac
def T_intersection_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    return d(m(Noplusa,m(m(48/5, rvala), np.power(m(snr0a, avara), 2))), m(np.power(m(iwafaca, la), 2), m(m(Ra, Ka), epsa)))

# No eq label, just a combo of below.
def d_prior_T(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Tarr, iwafaca, rvala):
    T_intersecta = T_intersection_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, [], iwafaca, rvala)
    intersect_ind =  int(np.where(np.isclose(Tarr,T_intersecta))[0])
    iwalim = d_iwa_all(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Tarr[:intersect_ind], iwafaca, rvala)
    iwalim2 = d_iwa_all(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Tarr[intersect_ind:], iwafaca, rvala)
    photonlim = d_photon_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Tarr[intersect_ind:], iwafaca, rvala)
    photonlim2 = d_photon_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Tarr[:intersect_ind], iwafaca, rvala)
    return np.concatenate((iwalim, photonlim)), np.concatenate((photonlim2, iwalim2))

# Eq Label: "eq:Tnoprioriwa_sub"
def T_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    mvar = d(m(16,m(rvala,m(snr0a,snr0a))),m(Ra,m(Ka,epsa)))
    mvar = m(mvar,np.power(d(3,m(4*pi,rhostara)),2/3))
    nvar = d(m(iwafaca,la),avara)
    sigvar = m(3/5,np.power(d(Noplusa,etaeartha),5/3))
    Dlim = np.cbrt(d(m(3,Noplusa),m(4*pi,m(rhostara, etaeartha))))
    #return d(m(mvar,sigvar),m(m(da,da),np.sqrt(np.abs(s(1,np.power(d(m(Dlim,nvar),da),2))))))
    return d(m(mvar,sigvar),m(m(da,da),np.sqrt(s(1,np.power(d(m(Dlim,nvar),da),2)))))

# Eq Label: "eq:Tsumknow_sub"
def T_photon_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    front = d(m(m(16,rvala),m(snr0a,snr0a)),m(m(Ra,Ka),m(m(da,da),epsa)))
    middle = np.power(d(3,m(4*pi,m(rhostara,etaeartha))),2/3)
    end = m(3/5,np.power(Noplusa,5/3))
    return m(front,m(middle,end))

T_iwa_prior=T_photon_prior

# Eq Label: "eq:noplusmaxiwa"
def get_Noplus_sca(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    cubed = np.power(d(m(avara,da),m(iwafaca,la)),3)
    modifier = m(4*pi/3,m(rhostara,etaeartha))
    return m(cubed, modifier)

# Eq Label: "eq:etaminiwa"
def get_etaearth_sca(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    cubed = np.power(d(m(iwafaca,la),m(avara,da)),3)
    modifier = d(m(3,Noplusa),m(4*pi,rhostara))
    return m(cubed, modifier)

# Eq Label: "eq:SNR0noprior_iwa"
def SNR_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, da, Ra, Ka, epsa, Ta, iwafaca, rvala):
    Dlim = np.cbrt(d(m(3,Noplusa),m(4*pi,m(rhostara, etaeartha))))
    innersqrt=np.sqrt(s(1,np.power(d(m(iwafaca,m(Dlim,la)),m(da,avara)),2)))
    top=m(m(m(5/48,Ra),m(Ka,Ta)),m(m(epsa,innersqrt),m(np.power(etaeartha,5/3),np.power(da,2))))
    bottom=m(rvala,np.power(Noplusa,5/3))
    outersqrt=np.sqrt(d(top,bottom))
    cuberoot=np.cbrt(d(m(4*pi,rhostara),3))
    return m(outersqrt,cuberoot)

# Eq Label: "eq:SNR0prior"
def SNR_photon_prior(Noplusa, la, avara, rhostara, etaeartha, da, Ra, Ka, epsa, Ta, iwafaca, rvala):
    top=m(m(5/48,np.power(da,2)),m(m(Ra,Ka),m(Ta,epsa)))
    bottom=m(rvala,np.power(Noplusa,5/3))
    squareroot=np.sqrt(d(top,bottom))
    cuberoot=np.cbrt(d(m(4*pi,m(rhostara,etaeartha)),3))
    return m(squareroot,cuberoot)

#Just a helper function, which uses np.polynomial to find a real positive root.
def get_real_pos_root(coeff0a,coeff2a,coeff10a):
    roots=np.polynomial.polynomial.Polynomial((coeff0a,0,coeff2a,0,0,0,0,0,0,0,coeff10a)).roots()
    realroots = roots[np.isreal(roots)]
    posrealroots = realroots[realroots>0]
    return npf64(posrealroots[0])

# Eq Label: "eq:Nstarpolynomial"
def Noplus_photon_noprior(da, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    mvar=m(d(m(m(16,rvala),m(snr0a,snr0a)),m(Ra,m(Ka,epsa))),np.power(d(3,m(4*pi,rhostara)),2/3))
    nvar=d(m(iwafaca,la),avara)
    alphsq=np.power(d(m(m(5,Ta),m(da,da)),m(3,mvar)),2)
    beta=m(np.power(d(nvar,da),2),np.power(d(3,m(4*pi,rhostara)),2/3))
    coeff0=-alphsq
    coeff2=m(alphsq,beta)
    coeff10=1
    #Check if coeff2 is iterable
    try:
        len(coeff2)
        results=np.array([])
        for i in range(len(coeff2)):
            #We don't know if alphsq or beta was iterable, so try alphsq.
            try:
                len(coeff0)
                c0=coeff0[i]
            except TypeError:
                c0=coeff0
            #coeff2 should always be iterable, and coeff10 is just 1.
            c2=coeff2[i]
            c10=coeff10
            #Now cube, and multiply by etaearth
            results=np.append(results,get_real_pos_root(c0,c2,c10))
    except TypeError:
        results=get_real_pos_root(coeff0,coeff2,coeff10)
    return m(np.power(results,3),etaeartha)

# Eq Label: "eq:noplusprior"
def Noplus_photon_prior(da, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    cubestop=np.power(m(m(Ra,Ka),m(Ta,epsa)),3)
    squarestop = np.power(m(m(pi,np.power(da,3)),m(rhostara,etaeartha)),2)
    bottom=np.power(m(rvala,m(snr0a,snr0a)),3)
    fraction=m(125/62208,d(m(cubestop,squarestop),bottom))
    try:
        len(fraction)
        return np.power(fraction,1/5)
    except:
        return [np.power(fraction,1/5) for i in m(m(la,avara),iwafaca)]

def Noplus_iwa_prior(da, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    cubes = np.power(d(m(avara,da),la),3)
    linear = m(4*pi/81, m(rhostara, etaeartha))
    res = m(linear, cubes)
    try:
        len(res)
        return res
    except:
        return [res for i in m(m(m(snr0a, Ra), m(Ka, epsa)),m(Ta,m(iwafaca,rvala)))]

# Now we get to actually graph all these things.

if NOPLUSvD:
    noplus_intersect = intersection_prior([], l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    NoplusVals = np.concatenate((np.array(range(int(noplus_intersect+1))),
                                 np.array([noplus_intersect]),
                                 np.array(range(int(noplus_intersect+1), int(8*noplus_intersect)))))
    dvals_p, dvals_pa  = d_prior_noplus(NoplusVals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    dvals_np = d_photon_noprior(NoplusVals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)

    ax=plt.figure()
    if LINEAR:
        plt.plot(NoplusVals[1:], dvals_p[1:], color='b', label='Prior Knowledge')
        plt.plot(NoplusVals[1:], dvals_pa[1:], color='b', linestyle='--')
        plt.plot(NoplusVals[1:], dvals_np[1:], color='r', label='No Prior Knowledge')
        plt.plot(Noplus, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Properties', markersize=3)
        plt.figtext(0.135, 0.56, 'IWA-Limited', fontsize=8);plt.figtext(0.54, 0.56, 'Photon-Limited', fontsize=8)
    elif LOGY:
        plt.semilogy(NoplusVals[1:], dvals_p[1:], color='b', label='Prior Knowledge')
        plt.semilogy(NoplusVals[1:], dvals_pa[1:], color='b', linestyle='--')
        plt.semilogy(NoplusVals[1:], dvals_np[1:], color='r', label='No Prior Knowledge')
        plt.semilogy(Noplus, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Properties', markersize=3)
        plt.figtext(0.135, 0.8, 'IWA-Limited', fontsize=8);plt.figtext(0.54, 0.65, 'Photon-Limited', fontsize=8)
    elif LOGX or UNSET:
        plt.semilogx(NoplusVals[1:], dvals_p[1:], color='b', label='Prior Knowledge')
        plt.semilogx(NoplusVals[1:], dvals_pa[1:], color='b', linestyle='--')
        plt.semilogx(NoplusVals[1:], dvals_np[1:], color='r', label='No Prior Knowledge')
        plt.semilogx(Noplus, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Properties', markersize=3)
        if HABEX:
            plt.figtext(0.45, 0.2, 'IWA-Limited', fontsize=8);plt.figtext(0.74, 0.2, 'Photon-Limited', fontsize=8)
        else:
            plt.figtext(0.45, 0.23, 'IWA-Limited', fontsize=8);plt.figtext(0.74, 0.23, 'Photon-Limited', fontsize=8)
    elif LOGLOG:
        plt.loglog(NoplusVals[1:], dvals_p[1:], color='b', label='Prior Knowledge')
        plt.loglog(NoplusVals[1:], dvals_pa[1:], color='b', linestyle='--')
        plt.loglog(NoplusVals[1:], dvals_np[1:], color='r', label='No Prior Knowledge')
        plt.loglog(Noplus, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Properties', markersize=3)
        plt.figtext(0.45, 0.23, 'IWA-Limited', fontsize=8);plt.figtext(0.74, 0.23, 'Photon-Limited', fontsize=8)
    # plt.axvline(x=Noplus, color='g', linestyle='--', label=HABEX*'HabEx' + LUVOIR*'LUVOIR'+' Expected Yield')
    plt.axvline(x=noplus_intersect, color='k', linestyle=':', label='IWA-Photon intersection (prior)')
    #plt.figtext(0.45, 0.23, 'IWA-Limited', fontsize=8);plt.figtext(0.74, 0.23, 'Photon-Limited', fontsize=8)
    # plt.axhline(y=dvar, color='c', linestyle='--', label=HABEX*'HabEx' + LUVOIR*'LUVOIR'+' Telescope Diameter')
    plt.ylabel('Minimum Telescope Diameter (m)'); plt.xlabel('EEC/HEC Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + 
                      ' Minimum Diameter vs. Yield, for $\eta_\oplus='+str(float(etaearth))+'$.'))
    plt.legend(prop={'size':6}); plt.ylim([-10, 200])
    plt.show()

if ETAvD:
    # Skip 1st element of both lists to prevent div0 and doubling of eta-earth.
    etavals = np.concatenate((np.linspace(0, etaearth, 100)[1:], np.linspace(etaearth, 1, 300)[1:]))
    noplus_intersect = intersection_prior([], l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    titlestr = 'Minimum Diameter vs. $\eta_\oplus$, for $N_\oplus='+str(int(Noplus))+'$.'
    if Noplus>=noplus_intersect:
        dvals_p = d_photon_prior(Noplus, l, avar, rhostar, etavals, snr0, R, K, eps, T, iwafac, rval)
        titlestr=' Photon-limited '+titlestr
    else:
        dvals_p = d_iwa_all(Noplus, l, avar, rhostar, etavals, snr0, R, K, eps, T, iwafac, rval)
        titlestr=' IWA-limited '+titlestr
    dvals_np = d_photon_noprior(Noplus, l, avar, rhostar, etavals, snr0, R, K, eps, T, iwafac, rval)
    ax = plt.figure()
    if LINEAR:
        plt.plot(etavals, dvals_p, color='b', label='Prior Knowledge')
        plt.plot(etavals, dvals_np, color='r', label='No Prior Knowledge')
        plt.plot(etaearth, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGY:
        plt.semilogy(etavals, dvals_p, color='b', label='Prior Knowledge')
        plt.semilogy(etavals, dvals_np, color='r', label='No Prior Knowledge')
        plt.semilogy(etaearth, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGX:
        plt.semilogx(etavals, dvals_p, color='b', label='Prior Knowledge')
        plt.semilogx(etavals, dvals_np, color='r', label='No Prior Knowledge')
        plt.semilogx(etaearth, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGLOG or UNSET:
        plt.loglog(etavals, dvals_p, color='b', label='Prior Knowledge')
        plt.loglog(etavals, dvals_np, color='r', label='No Prior Knowledge')
        plt.loglog(etaearth, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Minimum Telescope Diameter (m)'); plt.xlabel('$\eta_\oplus$')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    currlim=plt.ylim(); plt.ylim([currlim[0], 25])
    plt.legend(prop={'size':6}); plt.show()

if DvNOPLUS:
    noplusvals = list(range(1,1001))
    noplus_intersect = intersection_prior([], l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    noplusvals.append(noplus_intersect)
    noplusvals.sort()
    noplusvals = np.array(noplusvals)
    titlestr = ' Yield vs. Telescope Diameter, for $\eta_\oplus='+str(float(etaearth))+'$.'
    dvals_p, dvals_pa = d_prior_noplus(noplusvals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    dvals_np = d_photon_noprior(noplusvals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    ax = plt.figure()
    if LINEAR or UNSET:
        plt.plot(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.plot(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.plot(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.plot(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGY:
        plt.semilogy(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.semilogy(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.semilogy(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.semilogy(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGX:
        plt.semilogx(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.semilogx(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.semilogx(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.semilogx(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGLOG:
        plt.loglog(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.loglog(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.loglog(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.loglog(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)

    plt.xlabel('Telescope Diameter (m)'); plt.ylabel('Exo-Earth Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    # currlim=plt.ylim(); plt.ylim([currlim[0], 25])
    plt.legend(prop={'size':6}); plt.show()

if DvTIME:
    tintersect = T_intersection_prior(Noplus, l, avar, rhostar, etaearth, snr0, R, K, eps, [], iwafac, rval)
    timevals = list(np.linspace(0.1*tintersect, 10*T, 10000))
    tval2 = np.linspace(1, 2*tintersect, 1000)
    tval3 = np.linspace(1, 100000,100000)
    timevals.append(T)
    timevals.append(tintersect)
    timevals.sort()
    timevals = np.array(timevals)
    timevals = np.concatenate([timevals, tval2, tval3])
    timevals = np.sort(timevals)
    timevalsy = d(timevals, 3600*24*365)
    titlestr = ' Telescope Diameter vs. Survey Duration, for $\eta_\oplus='+str(float(etaearth))+'$.'
    dvals_p, dvals_pa = d_prior_T(Noplus, l, avar, rhostar, etaearth, snr0, R, K, eps, timevals, iwafac, rval)
    dvals_np = d_photon_noprior(Noplus, l, avar, rhostar, etaearth, snr0, R, K, eps, timevals, iwafac, rval)
    ax = plt.figure()
    if LINEAR:
        plt.plot(timevalsy, dvals_pa, color='b', label='Prior Knowledge')
        plt.plot(timevalsy, dvals_p, color='b', linestyle='--')
        plt.plot(timevalsy, dvals_np, color='r', label='No Prior Knowledge')
        plt.plot(T/(3600*24*365), dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGY:
        plt.semilogy(timevalsy, dvals_pa, color='b', label='Prior Knowledge')
        plt.semilogy(timevalsy, dvals_p, color='b', linestyle='--')
        plt.semilogy(timevalsy, dvals_np, color='r', label='No Prior Knowledge')
        plt.semilogy(T/(3600*24*365), dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGX:
        plt.semilogx(timevalsy, dvals_pa, color='b', label='Prior Knowledge')
        plt.semilogx(timevalsy, dvals_p, color='b', linestyle='--')
        plt.semilogx(timevalsy, dvals_np, color='r', label='No Prior Knowledge')
        plt.semilogx(T/(3600*24*365), dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGLOG or UNSET:
        plt.loglog(timevalsy, dvals_pa, color='b', label='Prior Knowledge')
        plt.loglog(timevalsy, dvals_p, color='b', linestyle='--')
        plt.loglog(timevalsy, dvals_np, color='r', label='No Prior Knowledge')
        plt.loglog(T/(3600*24*365), dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Telescope Diameter (m)'); plt.xlabel('Survey Duration (years)')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    # currlim=plt.ylim(); plt.ylim([currlim[0], 25])
    plt.legend(prop={'size':6}); plt.show()

if IWAvNOPLUS:
    noplusvals = list(range(1,1001))
    noplus_intersect = intersection_prior([], l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    noplusvals.append(noplus_intersect)
    noplusvals.sort()
    noplusvals = np.array(noplusvals)
    titlestr = ' Yield vs. Inner Working Angle, for $\eta_\oplus='+str(float(etaearth))+'$.'
    dvals_p, dvals_pa = d_prior_noplus(noplusvals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    dvals_np = d_photon_noprior(noplusvals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    ax = plt.figure()
    dvals_p = d(m(3,l),dvals_p);dvals_pa = d(m(3,l),dvals_pa);dvals_np = d(m(3,l),dvals_np);dvar = d(m(3,l),dvar);
    if LINEAR or UNSET:
        plt.plot(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.plot(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.plot(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.plot(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGY:
        plt.semilogy(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.semilogy(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.semilogy(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.semilogy(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGX:
        plt.semilogx(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.semilogx(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.semilogx(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.semilogx(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGLOG:
        plt.loglog(dvals_p, noplusvals, color='b', label='Prior Knowledge')
        plt.loglog(dvals_pa, noplusvals, color='b', linestyle='--')
        plt.loglog(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
        plt.loglog(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.xlabel('IWA (rad)'); plt.ylabel('Exo-Earth Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    # currlim=plt.ylim(); plt.ylim([currlim[0], 25])
    plt.legend(prop={'size':6}); plt.show()

if NOPLUSvT:
    titlestr = ' Yield vs. Survey Duration, for $\eta_\oplus='+str(float(etaearth))+'$.'
    noplus_limit = get_Noplus_sca([], l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    noplusvals = list(range(1,int(noplus_limit)-2))
    #Aggressively sample near the asymptote, to give it better clarity.
    nearval = np.linspace(int(noplus_limit)-2, noplus_limit, 1000)[0:999]
    noplusvals = np.concatenate((noplusvals, nearval))
    tvals_p = T_photon_prior(noplusvals, l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    tvals_np = T_photon_noprior(noplusvals, l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    ax = plt.figure()
    if LINEAR:
        plt.plot(noplusvals, d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.plot(noplusvals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
        plt.axvline(x=noplus_limit, color='y', linestyle='--',label='$s_c=a$ limit on $N_\oplus$')
        plt.plot(Noplus,T/(3600*24*365), 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGY:
        plt.semilogy(noplusvals, d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.semilogy(noplusvals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
        plt.axvline(x=noplus_limit, color='y', linestyle='--',label='$s_c=a$ limit on $N_\oplus$')
        plt.semilogy(Noplus,T/(3600*24*365), 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGX:
        plt.semilogx(noplusvals, d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.semilogx(noplusvals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
        plt.axvline(x=noplus_limit, color='y', linestyle='--',label='$s_c=a$ limit on $N_\oplus$')
        plt.semilogx(Noplus,T/(3600*24*365), 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif LOGLOG or UNSET:
        plt.loglog(noplusvals, d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.loglog(noplusvals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
        plt.axvline(x=noplus_limit, color='y', linestyle='--',label='$s_c=a$ limit on $N_\oplus$')
        plt.loglog(Noplus,T/(3600*24*365), 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Survey duration (years)'); plt.xlabel('Exo-Earth Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()

if ETAvT:
    titlestr = ' Yield vs. $\eta_\oplus$'
    minetaval = get_etaearth_sca(Noplus, l, avar, rhostar, [], snr0, R, K, eps, dvar, iwafac, rval)
    etavals = np.linspace(minetaval,1,1000)[1:999]
    tvals_p = T_photon_prior(Noplus, l, avar, rhostar, etavals, snr0, R, K, eps, dvar, iwafac, rval)
    tvals_np = T_photon_noprior(Noplus, l, avar, rhostar, etavals, snr0, R, K, eps, dvar, iwafac, rval)
    ax = plt.figure()
    plt.xlim([0,1])
    if LINEAR:
        plt.plot(etavals,d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.plot(etavals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
    elif LOGY or UNSET:
        plt.semilogy(etavals,d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.semilogy(etavals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
    elif LOGX:
        plt.semilogx(etavals,d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.semilogx(etavals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
    elif LOGLOG:
        plt.loglog(etavals,d(tvals_p,3600*24*365), color='b', label='Prior Knowledge')
        plt.loglog(etavals, d(tvals_np,3600*24*365), color='r', label='No Prior Knowledge')
    plt.axvline(x=minetaval, color='y', linestyle='--',label='$s_c=a$ limit on $\eta_\oplus$')
    plt.plot(etaearth,T/(3600*24*365), 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Survey duration (years)'); plt.xlabel('$\eta_\oplus$')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()

if SNR0vNOPLUS:
    titlestr = ' $SNR_0$ vs $N_\oplus$'
    noplus_limit = get_Noplus_sca([], l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    noplusvals = np.linspace(1,int(noplus_limit)-2,1000)
    #Aggressively sample near the asymptote, to give it better clarity.
    nearval = np.linspace(int(noplus_limit)-2, noplus_limit, 1000)[0:999]
    noplusvals = np.concatenate((noplusvals, nearval))
    snrvals_p = SNR_photon_prior(noplusvals, l, avar, rhostar, etaearth, dvar, R, K, eps, T, iwafac, rval)
    snrvals_np = SNR_photon_noprior(noplusvals, l, avar, rhostar, etaearth, dvar, R, K, eps, T, iwafac, rval)
    ax = plt.figure()
    if LINEAR:
        plt.plot(noplusvals, snrvals_p, color='b', label='Prior Knowledge')
        plt.plot(noplusvals, snrvals_np, color='r', label='No Prior Knowledge')
    elif LOGY or UNSET:
        plt.semilogy(noplusvals, snrvals_p, color='b', label='Prior Knowledge')
        plt.semilogy(noplusvals, snrvals_np, color='r', label='No Prior Knowledge')
    elif LOGX:
        plt.semilogx(noplusvals, snrvals_p, color='b', label='Prior Knowledge')
        plt.semilogx(noplusvals, snrvals_np, color='r', label='No Prior Knowledge')
    elif LOGLOG:
        plt.loglog(noplusvals, snrvals_p, color='b', label='Prior Knowledge')
        plt.loglog(noplusvals, snrvals_np, color='r', label='No Prior Knowledge')
    plt.plot(Noplus,snr0, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('$SNR_0$'); plt.xlabel('Exo-Earth Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()
    
if ETAvNOPLUS:
    titlestr = ' $\eta_\oplus$ vs $N_\oplus$'
    etavals = np.linspace(0,1,10000)[1:9999]
    nvals_np=Noplus_photon_noprior(dvar, l, avar, rhostar, etavals, snr0, R, K, eps, T, iwafac, rval)
    nvals_pp=Noplus_photon_prior(dvar, l, avar, rhostar, etavals, snr0, R, K, eps, T, iwafac, rval)
    nvals_ip=Noplus_iwa_prior(dvar, l, avar, rhostar, etavals, snr0, R, K, eps, T, iwafac, rval)
    nvals_p = [min(nvals_pp[i], nvals_ip[i]) for i in range(len(nvals_ip))]
    ax = plt.figure()
    plt.xlim([0,1])
    if LINEAR or UNSET:
        plt.plot(etavals,nvals_p, color='b', label='Prior Knowledge')
        plt.plot(etavals, nvals_np, color='r', label='No Prior Knowledge')
        plt.gca().set_ylim(bottom=0)
    elif LOGY:
        plt.semilogy(etavals,nvals_p, color='b', label='Prior Knowledge')
        plt.semilogy(etavals, nvals_np, color='r', label='No Prior Knowledge')
    elif LOGX:
        plt.semilogx(etavals,nvals_p, color='b', label='Prior Knowledge')
        plt.semilogx(etavals, nvals_np, color='r', label='No Prior Knowledge')
    elif LOGLOG:
        plt.loglog(etavals,nvals_p, color='b', label='Prior Knowledge')
        plt.loglog(etavals, nvals_np, color='r', label='No Prior Knowledge')
    
    plt.plot(etaearth,Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Survey Yield'); plt.xlabel('$\eta_\oplus$')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()
if EPSvNOPLUS:
    titlestr=' Survey Efficiency vs $N_\oplus$'
    epsvals=np.linspace(0,1,10000)[1:9999]
    nvals_np=Noplus_photon_noprior(dvar, l, avar, rhostar, etaearth, snr0, R, K, epsvals, T, iwafac, rval)
    nvals_pp=Noplus_photon_prior(dvar, l, avar, rhostar, etaearth, snr0, R, K, epsvals, T, iwafac, rval)
    nvals_ip=Noplus_iwa_prior(dvar, l, avar, rhostar, etaearth, snr0, R, K, epsvals, T, iwafac, rval)
    nvals_p = [min(nvals_pp[i], nvals_ip[i]) for i in range(len(nvals_ip))]
    if LINEAR or UNSET:
        plt.plot(epsvals,nvals_p, color='b', label='Prior Knowledge')
        plt.plot(epsvals, nvals_np, color='r', label='No Prior Knowledge')
        plt.gca().set_ylim(bottom=0)
        plt.xlim([0,1])
    elif LOGY:
        plt.semilogy(epsvals,nvals_p, color='b', label='Prior Knowledge')
        plt.semilogy(epsvals, nvals_np, color='r', label='No Prior Knowledge')
    elif LOGX:
        plt.semilogx(epsvals,nvals_p, color='b', label='Prior Knowledge')
        plt.semilogx(epsvals, nvals_np, color='r', label='No Prior Knowledge')
    elif LOGLOG:
        plt.loglog(epsvals,nvals_p, color='b', label='Prior Knowledge')
        plt.loglog(epsvals, nvals_np, color='r', label='No Prior Knowledge')
    
    plt.plot(eps, Noplus, 'go', label='Assumption', markersize=3)
    plt.ylabel('Survey Yield'); plt.xlabel('Survey Efficiency')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()

if IWAFACvNOPLUS:
    titlestr = ' IWA ($\lambda$/d) vs $N_\oplus$'
    iwavals = np.linspace(1,6,10000)
    nvals_np=Noplus_photon_noprior(dvar, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwavals, rval)
    nvals_pp=Noplus_photon_prior(dvar, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwavals, rval)
    nvals_ip=Noplus_iwa_prior(dvar, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwavals, rval)
    nvals_p = [min(nvals_pp[i], nvals_ip[i]) for i in range(len(nvals_ip))]
    ax = plt.figure()
    if LINEAR or UNSET:
        plt.plot(iwavals,nvals_p, color='b', label='Prior Knowledge')
        plt.plot(iwavals, nvals_np, color='r', label='No Prior Knowledge')
        plt.gca().set_ylim(bottom=0)
    elif LOGY:
        plt.semilogy(iwavals,nvals_p, color='b', label='Prior Knowledge')
        plt.semilogy(iwavals, nvals_np, color='r', label='No Prior Knowledge')
    elif LOGX:
        plt.semilogx(iwavals,nvals_p, color='b', label='Prior Knowledge')
        plt.semilogx(iwavals, nvals_np, color='r', label='No Prior Knowledge')
    elif LOGLOG:
        plt.loglog(iwavals,nvals_p, color='b', label='Prior Knowledge')
        plt.loglog(iwavals, nvals_np, color='r', label='No Prior Knowledge')
    
    plt.plot(iwafac, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Survey Yield'); plt.xlabel('IWA ($\lambda$/d)')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()
# Sources: HabEx
# l:              (HabEx Final Report 3-6) (Standards Team Final Report Table 5)
# eta-earth:      (HabEx Final Report 3-3) (Belikov 2017)
# res:            (
# T:
# eps:
# dvar:
# iwafac:
# Noplus:
# snr0:
# Sources: HabEx
# l:              
# eta-earth:
# res:
# T:
# eps:
# dvar:
# iwafac:
# Noplus:
# snr0:
