__author__="huziy"
__date__ ="$Jun 15, 2011 3:57:03 PM$"

from math import *

#f([sigma, mu, ksi], x) = 1 - (ksi/sigma) * (x - mu)
def yi(params, value):
    sigma, mu, ksi = params
    return 1 - (ksi / sigma) * (value - mu)

##derivative of the -ln(joint probability density function)
## after Stedinger (2000)
def deriv_of_objective_func_stedinger(params, extremes):
    sigma, mu, ksi = params
    ksi = -ksi
    params1 = [sigma, mu, ksi]

    d1 = 0.0
    d2 = - len(extremes)
    d3 = 0.0
    for extreme in extremes:
        y = yi(params1, extreme)
        similar_term = ( 1.0 - ksi - y ** ksi)
        d1 +=  similar_term / y
        d2 += similar_term / y * (1.0 - y) / ksi
        d3 += ( log(y) + 1 - y ) * similar_term

    d1 *= 1.0 / sigma
    d2 *= 1.0 / sigma
    d3 *= 1.0 / ksi ** 2

    return [d1, d2, d3]

    pass

def test_deriv_stedinger():
    print deriv_of_objective_func_stedinger([1000,1,1], [1])

if __name__ == "__main__":
    test_deriv_stedinger()
    print "Hello World"
