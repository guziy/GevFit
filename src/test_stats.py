__author__ = "huziy"
__date__ = "$16 fevr. 2011 17:24:40$"


from scipy.stats.distributions import genextreme_gen

def test_func():

    data = range(100, 150, 1)
    x = genextreme_gen(name = 'genextreme', longname = 'A generalized extreme value', shapes = 'c',
                extradoc = """ Generalized extreme value (see gumbel_r for c=0)
                genextreme.pdf(x,c) = exp(-(1-c*x)**(1/c))*(1-c*x)**(1/c-1)
                for x <= 1/c, c > 0 """)
    x.fit(data)
    help(x.fit)
    print x.__dict__

test_func()





if __name__ == "__main__":
    print "Hello World"
