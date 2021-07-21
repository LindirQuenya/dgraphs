#!/usr/bin/env python3
#from mpmath import mp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# First one uses script location, second uses cwd. Comment/uncomment to taste.
mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(__file__))
# mpl.rcParams["savefig.directory"] = ""

# The following specifies that titles should be on/off for figures.
TITLES=False

def print_help():
    print("""
Usage: python3 dgraphs.py [STUDY] [COMPARISON]
Generate a graph that compares two variables for a given mission concept.

  -h, --help     display this message

  STUDY          one of [habex, luvoir, luvoiralt]
  COMPARISON     integer, [0..3] such that: (x vs y)
                          0 -> N_oplus   vs d
                          1 -> eta_oplus vs d
                          2 -> d         vs N_oplus
                          3 -> d         vs duration
                          4 -> eta_oplus vs N_oplus  (Not Implemented)
                          5 -> IWA       vs N_oplus
                          6 -> SNR0      vs N_oplus
                          7 -> N_oplus   vs duration
                          8 -> eta_oplus vs duration

If an invalid value is supplied for STUDY, it will
default to habex.
If an invalid value is supplied for COMPARISON,
the program will display this message and exit.""")

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


NOPLUSvD, ETAvD, DvNOPLUS, DvTIME, ETAvNOPLUS, IWAvNOPLUS, SNR0vNOPLUS, NOPLUSvT, ETAvT = [False]*9
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
    else:
        print_help()
        #exit(1)
else:
    print_help()
    #exit(1)

#mp.dps=50
#mp.prec=171

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
def d_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    n = d(m(iwafaca, la), avara)
    cuberoot = np.cbrt(d(m(3, Noplusa), m(4*pi, m(rhostara, etaeartha))))
    is_low = m(5, m(Ra, m(Ka, m(epsa, m(np.power(n, 2), m(Ta, etaeartha))))))
    is_high = m(48, m(m(rvala,np.power(snr0a, 2)), Noplusa))
    innersqrt = np.sqrt(a(1, np.power(d(is_high, is_low), 2)))
    outersqrt = np.sqrt(a(1, innersqrt))
    return m(d(n, np.sqrt(2)), m(cuberoot, outersqrt))

def iwa_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    return d(m(3,la),d_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala))

def d_photon_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    cuberoot = np.cbrt(d(3, m(4*pi, m(rhostara, etaeartha))))
    onlysqrt = np.sqrt(d(m(m(3, rvala), np.power(Noplusa, 5/3)), m(5, m(Ra, m(Ka, m(Ta, epsa))))))
    return m(m(4, snr0a), m(cuberoot, onlysqrt))

def d_iwa_all(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    cuberoot = np.cbrt(d(m(3, Noplusa), m(4*pi, m(rhostara, etaeartha))))
    return m(d(m(iwafaca, la), avara), cuberoot)

def intersection_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    return m(m(5/48, np.power(iwafaca, 2)), d(m(Ra, m(Ka, m(Ta, m(epsa, np.power(la, 2))))), m(rvala, np.power(m(snr0a, avara), 2))))

def d_prior_noplus(Noplusarr, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala):
    noplus_intersecta = intersection_prior([], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    intersect_ind = int(np.where(Noplusarr == noplus_intersecta)[0])
    iwalim = d_iwa_all(Noplusarr[:intersect_ind], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    iwalim2 = d_iwa_all(Noplusarr[intersect_ind:], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    photonlim = d_photon_prior(Noplusarr[intersect_ind:], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    photonlim2 = d_photon_prior(Noplusarr[:intersect_ind], la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, Ta, iwafaca, rvala)
    return np.concatenate((iwalim, photonlim)), np.concatenate((photonlim2, iwalim2))

def T_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    mvar = d(m(16,m(rvala,m(snr0a,snr0a))),m(Ra,m(Ka,epsa)))
    mvar = m(mvar,np.power(d(3,m(4*pi,rhostara)),2/3))
    nvar = d(m(3,la),avara)
    sigvar = m(3/5,np.power(d(Noplusa,etaeartha),5/3))
    Dlim = np.cbrt(d(m(3,Noplusa),m(4*pi,m(rhostara, etaeartha))))
    #return d(m(mvar,sigvar),m(m(da,da),np.sqrt(np.abs(s(1,np.power(d(m(Dlim,nvar),da),2))))))
    return d(m(mvar,sigvar),m(m(da,da),np.sqrt(s(1,np.power(d(m(Dlim,nvar),da),2)))))

def T_photon_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    front = d(m(m(16,rvala),m(snr0a,snr0a)),m(m(Ra,Ka),m(m(da,da),epsa)))
    middle = np.power(d(3,m(4*pi,m(rhostara,etaeartha))),2/3)
    end = m(3/5,np.power(Noplusa,5/3))
    return m(front,m(middle,end))

T_iwa_prior=T_photon_prior

def get_sc(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    Dlim = np.cbrt(d(m(3,Noplusa),m(4*pi,m(rhostara, etaeartha))))
    d(m(3,m(Dlim,la)),da)

def get_Noplus_sca(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    cubed = np.power(d(m(avara,da),m(3,la)),3)
    modifier = m(4*pi/3,m(rhostara,etaeartha))
    return m(cubed, modifier)

def get_etaearth_sca(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
    cubed = np.power(d(m(3,la),m(avara,da)),3)
    modifier = d(m(3,Noplusa),m(4*pi,rhostara))
    return m(cubed, modifier)

def SNR_photon_noprior(Noplusa, la, avara, rhostara, etaeartha, da, Ra, Ka, epsa, Ta, iwafaca, rvala):
    Dlim = np.cbrt(d(m(3,Noplusa),m(4*pi,m(rhostara, etaeartha))))
    nvar = d(m(3,la),avara)
    sq = s(1,np.power(d(m(Dlim,nvar),da),2))
    extramult=d(m(5,m(m(Ta,Ra),m(Ka,epsa))),m(48,r))
    fivethirds = np.power(d(etaeartha,Noplusa),5/3)
    twothirds = np.power(d(m(4*pi,rhostara),3),2/3)
    return np.sqrt(m(m(sq, extramult),m(fivethirds, twothirds)))

def SNR_photon_prior(Noplusa, la, avara, rhostara, etaeartha, da, Ra, Ka, epsa, Ta, iwafaca, rvala):
    top=m(m(5,Ta),m(m(Ra,Ka),m(np.power(da,2),epsa)))
    bottom=m(48,m(r,np.power(Noplusa,5/3)))
    twothirds = np.power(d(m(4*pi,m(rhostara,etaeartha)),3),2/3)
    return m(d(top, bottom),twothirds)

#WRONG!
#def T_iwa_prior(Noplusa, la, avara, rhostara, etaeartha, snr0a, Ra, Ka, epsa, da, iwafaca, rvala):
#    top = m(16,m(m(rvala,m(snr0a,snr0a)),m(Noplusa,m(avara,avara))))
#    bottom = m(15,m(m(Ra,Ka),m(epsa,m(la,la))))
#    return d(top,bottom)

if NOPLUSvD:
    noplus_intersect = intersection_prior([], l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    NoplusVals = np.concatenate((np.array(range(int(noplus_intersect+1))),
                                 np.array([noplus_intersect]),
                                 np.array(range(int(noplus_intersect+1), int(8*noplus_intersect)))))
    dvals_p, dvals_pa  = d_prior_noplus(NoplusVals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    dvals_np = d_photon_noprior(NoplusVals, l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)

    ax=plt.figure()  # ; pltfname=HABEX*'HabEx' + LUVOIR*'LUVOIR'+'plt.pickle'
    if True:
        plt.semilogx(NoplusVals[1:], dvals_p[1:], color='b', label='Prior Knowledge')
        plt.semilogx(NoplusVals[1:], dvals_pa[1:], color='b', linestyle='--')
        plt.semilogx(NoplusVals[1:], dvals_np[1:], color='r', label='No Prior Knowledge')
        plt.semilogx(Noplus, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Properties', markersize=3)
    else:
        plt.loglog(NoplusVals[1:], dvals_p[1:], color='b', label='Prior Knowledge')
        plt.loglog(NoplusVals[1:], dvals_pa[1:], color='b', linestyle='--')
        plt.loglog(NoplusVals[1:], dvals_np[1:], color='r', label='No Prior Knowledge')
        plt.loglog(Noplus, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Properties', markersize=3)
    # plt.axvline(x=Noplus, color='g', linestyle='--', label=HABEX*'HabEx' + LUVOIR*'LUVOIR'+' Expected Yield')
    plt.axvline(x=noplus_intersect, color='k', linestyle=':', label='IWA-Photon intersection (prior)')
    plt.figtext(0.45, 0.23, 'IWA-Limited', fontsize=8);plt.figtext(0.74, 0.23, 'Photon-Limited', fontsize=8)
    # plt.axhline(y=dvar, color='c', linestyle='--', label=HABEX*'HabEx' + LUVOIR*'LUVOIR'+' Telescope Diameter')
    plt.ylabel('Minimum Telescope Diameter (m)'); plt.xlabel('EEC/HEC Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + 
                      ' Minimum Diameter vs. Yield, for $\eta_\oplus='+str(float(etaearth))+'$.'))
    plt.legend(prop={'size':6}); plt.ylim([-10, 200])
    plt.show() # pickle.dump(ax, open(pltfname,'wb'))

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
    if False:
        plt.plot(etavals, dvals_p, color='b', label='Prior Knowledge')
        plt.plot(etavals, dvals_np, color='r', label='No Prior Knowledge')
        plt.plot(etaearth, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    elif False:
        plt.semilogy(etavals, dvals_p, color='b', label='Prior Knowledge')
        plt.semilogy(etavals, dvals_np, color='r', label='No Prior Knowledge')
        plt.semilogy(etaearth, dvar, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    else:
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
    plt.plot(dvals_p, noplusvals, color='b', label='Prior Knowledge')
    plt.plot(dvals_pa, noplusvals, color='b', linestyle='--')
    plt.plot(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
    plt.plot(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.xlabel('Telescope Diameter (m)'); plt.ylabel('Exo-Earth Yield')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    # currlim=plt.ylim(); plt.ylim([currlim[0], 25])
    plt.legend(prop={'size':6}); plt.show()

if DvTIME:
    timevals = list(np.linspace(0.1*T, 10*T, 1000))
    timevals.append(T)
    timevals.sort()
    timevals = np.array(timevals)
    timevalsy = d(timevals, 3600*24*365)
    titlestr = ' Survey Duration vs. Telescope Diameter, for $\eta_\oplus='+str(float(etaearth))+'$.'
    dvals_p, dvals_pa = d_prior_noplus(Noplus, l, avar, rhostar, etaearth, snr0, R, K, eps, timevals, iwafac, rval)
    dvals_np = d_photon_noprior(Noplus, l, avar, rhostar, etaearth, snr0, R, K, eps, timevals, iwafac, rval)
    ax = plt.figure()
    plt.plot(dvals_p, timevals, color='b', label='Prior Knowledge')
    plt.plot(dvals_pa, timevals, color='b', linestyle='--')
    plt.plot(dvals_np, timevals, color='r', label='No Prior Knowledge')
    plt.plot(dvar, T, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.xlabel('Telescope Diameter (m)'); plt.ylabel('Survey Duration (seconds)')
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
    plt.plot(dvals_p, noplusvals, color='b', label='Prior Knowledge')
    plt.plot(dvals_pa, noplusvals, color='b', linestyle='--')
    plt.plot(dvals_np, noplusvals, color='r', label='No Prior Knowledge')
    plt.plot(dvar, Noplus, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
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
    #noplusvals = list(range(1,1001))
    #noplus_intersect = intersection_prior([], l, avar, rhostar, etaearth, snr0, R, K, eps, T, iwafac, rval)
    #NoplusVals = np.concatenate((np.array(range(int(noplus_intersect+1))),
    #                             np.array([noplus_intersect]),
    #                             np.array(range(int(noplus_intersect+1), int(8*noplus_intersect)))))
    #tvals_p_pho = T_photon_prior(np.array(range(int(noplus_intersect+1))), l, avar, rhostar,
    #                             etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    #tvals_p_iwa = T_iwa_prior(np.array(range(int(noplus_intersect+1), int(8*noplus_intersect))),
    #                          l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    #tvals_np = T_photon_noprior(NoplusVals, l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    #noplusvals = np.concatenate((np.array(range(int(noplus_intersect+1))),
    #                             np.array(range(int(noplus_intersect+1), int(8*noplus_intersect)))))
    #tvals_p = np.concatenate((tvals_p_pho,tvals_p_iwa))
    ax = plt.figure()
    plt.loglog(noplusvals,tvals_p, color='b', label='Prior Knowledge')
    plt.loglog(noplusvals, tvals_np, color='r', label='No Prior Knowledge')
    plt.axvline(x=noplus_limit, color='y', linestyle='--',label='$s_c=a$ limit on $N_\oplus$')
    plt.plot(Noplus,T, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Survey duration (seconds)'); plt.xlabel('Exo-Earth Yield')
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
    plt.semilogy(etavals,tvals_p, color='b', label='Prior Knowledge')
    plt.semilogy(etavals, tvals_np, color='r', label='No Prior Knowledge')
    plt.axvline(x=minetaval, color='y', linestyle='--',label='$s_c=a$ limit on $\eta_\oplus$')
    plt.plot(etaearth,T, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('Survey duration (seconds)'); plt.xlabel('$\eta_\oplus$')
    plt.title(TITLES*(HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR' + titlestr))
    plt.legend(prop={'size':6}); plt.show()

if SNR0vNOPLUS:
    titlestr = ' $SNR_0$ vs $N_\oplus$'
    noplus_limit = get_Noplus_sca([], l, avar, rhostar, etaearth, snr0, R, K, eps, dvar, iwafac, rval)
    noplusvals = list(range(1,int(noplus_limit)-2))
    #Aggressively sample near the asymptote, to give it better clarity.
    nearval = np.linspace(int(noplus_limit)-2, noplus_limit, 1000)[0:999]
    noplusvals = np.concatenate((noplusvals, nearval))
    snrvals_p = SNR_photon_prior(noplusvals, l, avar, rhostar, etaearth, dvar, R, K, eps, T, iwafac, rval)
    snrvals_np = SNR_photon_noprior(noplusvals, l, avar, rhostar, etaearth, dvar, R, K, eps, T, iwafac, rval)
    ax = plt.figure()
    plt.semilogy(noplusvals, snrvals_p, color='b', label='Prior Knowledge')
    plt.semilogy(noplusvals, snrvals_np, color='r', label='No Prior Knowledge')
    plt.plot(Noplus,snr0, 'go', label=HABEX*'HabEx' + (LUVOIR or LUVOIRALT)*'LUVOIR'+' Assumption', markersize=3)
    plt.ylabel('$SNR_0$'); plt.xlabel('Exo-Earth Yield')
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
