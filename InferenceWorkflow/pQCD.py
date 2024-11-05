#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adopted from https://zenodo.org/records/7781233
"""

from numpy import log, pi
from numba import njit


GeV3_to_fm3 = 1.0e3/1.9732705**3


@njit
def PNLO(a_s):
	return 1. - 0.637*a_s
@njit
def PNNLO(a_s,X):
	return -a_s**2*(-1.831 + 0.304*log(a_s)) + a_s**2*(-2.706 - 0.912*log(X))
@njit
def PN3LO(a_s):
	return 0.484816*a_s**3
@njit
def alpha_s(mu,X):
	numerator = 4*pi*(1. - (64.*log(log(0.777632*mu**2*X**2)))/(81.*log(0.777632*mu**2*X**2)))
	denominator = (9.*log(0.777632*mu**2*X**2))
	return numerator/denominator
@njit
def das_dmu(mu,X):
	numerator = -2.20644 - 2.79253*log(0.777632*mu**2*X**2) + 4.41288*log(log(0.777632*mu**2*X**2))
	denominator = mu*(log(0.777632*mu**2*X**2))**3
	return numerator/denominator
@njit
def d2as_dmu2(mu,X):
	numerator1 = 22.0644 + 2.79253*(log(0.777632*mu**2*X**2))**2 - 26.4773*log(log(0.777632*mu**2*X**2))
	numerator2 = log(0.777632*mu**2*X**2)*(13.3765 - 4.41288*log(log(0.777632*mu**2*X**2)))
	denominator = mu**2*(log(0.777632*mu**2*X**2))**4
	return (numerator1+numerator2)/denominator
@njit
def dp_das(a_s,X):
	dPNLO = -0.637
	dPNNLO = a_s*(-2.054 - 0.608*log(a_s) - 1.824*log(X))
	dPN3LO = 1.45445*a_s**2 
	return dPNLO + dPNNLO + dPN3LO
@njit
def d2p_das2(a_s,X):
	d2PNNLO = -2.662 - 0.608*log(a_s) - 1.824*log(X)
	d2PN3LO = 2.9089*a_s 
	return d2PNNLO + d2PN3LO
@njit
def pFD(mu):
	return (mu)**4/(108*pi**2)
@njit
def dpFD(mu):
	return mu**3/(27*pi**2)
@njit
def d2pFD(mu):
	return mu**2/(9*pi**2)
@njit
def n_dens_QCD(mu,X): #fm-3

		a_s = alpha_s(mu,X)
		p_as = (PNLO(a_s) + PNNLO(a_s,X) + PN3LO(a_s))

		return (dp_das(a_s,X)*das_dmu(mu,X)*pFD(mu) + p_as*dpFD(mu))*GeV3_to_fm3
@njit
def speed2_QCD(mu,X):

		a_s = alpha_s(mu,X)
		p_as = (PNLO(a_s) + PNNLO(a_s,X) + PN3LO(a_s))

		dn_dmu1 = pFD(mu)*(d2p_das2(a_s,X)*(das_dmu(mu,X))**2 + dp_das(a_s,X)*d2as_dmu2(mu,X))
		dn_dmu2 = 2.*dp_das(a_s,X)*das_dmu(mu,X)*dpFD(mu) + p_as * d2pFD(mu)

		return n_dens_QCD(mu,X)/(mu*GeV3_to_fm3*(dn_dmu1+dn_dmu2))

@njit
def pressure_QCD(mu, X): #MeV/fm-3

		a_s = alpha_s(mu, X)

		return (PNLO(a_s) + PNNLO(a_s, X) + PN3LO(a_s))*pFD(mu)*GeV3_to_fm3*1.e3


@njit
def pressure(mu, X):#GeV/fm3
	a_s = alpha_s(mu, X)
	return (PNLO(a_s) + PNNLO(a_s, X) + PN3LO(a_s))*pFD(mu)*GeV3_to_fm3

@njit
def number_density(mu, X): #fm-3
	a_s = alpha_s(mu, X)
	p_as = (PNLO(a_s) + PNNLO(a_s, X) + PN3LO(a_s))
	return (dp_das(a_s, X)*das_dmu(mu, X)*pFD(mu) + p_as*dpFD(mu))*GeV3_to_fm3

@njit
def edens(mu): #GeV/fm3
	return -pressure(mu)+mu*number_density(mu)

@njit
def speed2(mu, X):
	a_s = alpha_s(mu, X)
	p_as = (PNLO(a_s) + PNNLO(a_s, X) + PN3LO(a_s))

	dn_dmu1 = pFD(mu)*(d2p_das2(a_s, X)*(das_dmu(mu, X))**2 + dp_das(a_s, X)*d2as_dmu2(mu, X))
	dn_dmu2 = 2.*dp_das(a_s, X)*das_dmu(mu, X)*dpFD(mu) + p_as * d2pFD(mu)

	return number_density(mu)/(mu*GeV3_to_fm3*(dn_dmu1+dn_dmu2))

@njit
def constraints(X, e0, p0, n0, muQCD = 2.6, cs2=1): 
	mu0 = (e0 + p0) / n0
	pQCD = pressure(muQCD, X)
	nQCD = number_density(muQCD, X)
	deltaP = pQCD - p0
	pMin = cs2 / (1.0 + cs2) * (muQCD * (muQCD / mu0) ** (1. / cs2) - mu0) * n0
	pMax = cs2 / (1.0 + cs2) * (muQCD - mu0 * (mu0 / muQCD) ** (1. / cs2)) * nQCD
	nMax = nQCD * (mu0 / muQCD) ** (1. / cs2)

	tag = (pMin < deltaP < pMax) and (n0 < nMax)
	return tag
