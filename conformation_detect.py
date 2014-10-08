#!/usr/bin/python

def detect_conf(phi, psi):
	if phi < 0 and psi < 50 and psi > -120:
		return r'$\alpha$'
	elif phi < -120 :
		return r'$\beta$'
	elif phi < -80:
		return r'$\beta/P_{II}$'
	elif phi < -20 :
		return r'$P_{II}$'
	elif psi < 60 and psi > -100:
		return r'$\alpha_L$'
	else:
		return r'$Y$'


def detect_conf_print(phi,psi):
	if phi < 0 and psi < 50 and psi > -120:
		return 'Alpha'
	elif phi < -90 :
		return 'Beta'
	elif phi < -20 :
		return 'PII'
	elif psi < 60 and psi > -100:
		return 'Alpha_L'
	else:
		return 'Y'
