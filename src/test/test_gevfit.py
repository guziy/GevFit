
import os
import application_properties
import gevfit
import members
import pickle
import numpy as np

import numpy.random as random

import matplotlib.pyplot as plt


def gev_cdf_inverse(p, sigma, mu, ksi):
    y = (-np.log(p)) ** (-ksi)
    return sigma / ksi * (y - 1.0) + mu

def gev_cdf(x, sigma, mu, ksi):
    y = (1 + ksi * (x - mu) / sigma) ** (-1.0 / ksi)
    return np.exp(-y)


class  GevfitTest():
    def setUp(self):
        application_properties.set_current_directory()
        #the files with dumped parameters using pickle.dump()
        self.gml_result_folder = 'old_rl/rl_14_GML_half_and_half'
        self.lm_result_folder = 'old_rl/rl_7_LM'

        self.gml_parameters = {}
        self.lm_parameters = {}


        types = ['high', 'low']
        prefix = 'gev_params_stationary'

        name_pattern = prefix + '_%s_%s'

        for the_type in types:
            lm = []
            gml = []

            for the_id in members.all_members:
                #read lm and gml params
                path = os.path.join(self.lm_result_folder, name_pattern % (the_id, the_type))
                lm_pars_list = pickle.load(open(path))
                
                path = os.path.join(self.gml_result_folder, name_pattern % (the_id, the_type))
                gml_pars_list = pickle.load(open(path))
                
                for lm_pars, gml_pars in zip(lm_pars_list, gml_pars_list):
                    if lm_pars[0] == None or gml_pars[0] == None:
                        continue
                    lm.append(lm_pars[0:3])
                    gml.append(gml_pars[0:3])


            ##save to the object's fields
            print the_type
            assert len(gml) == len(lm), 'Not None: gml = %d, lm = %d' % (len(gml), len(lm))
            self.gml_parameters[the_type] = np.array(gml)
            self.lm_parameters[the_type] = np.array(lm)

            pass



    def plot_scatter(self, xs, ys, xlabel = 'x', ylabel = 'y'):
        plt.scatter(xs, ys, linewidth = 0)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        vals = []
        vals.extend(xs)
        vals.extend(ys)

        the_min = 0.1 * np.min(vals)
        the_max = 1.1 * np.max(vals)
        plt.xlim(the_min, the_max)
        plt.ylim(the_min, the_max)

        plt.plot([the_min, the_max],[the_min, the_max], color = 'k')
        #plt.savefig('%s_%s_scatter.png' % (xlabel, ylabel))
        pass


    #compare GML and LM fitting approaches
    #read saved parameters obtained using the 2 methods
    #and plot scatter plot
    #data in the dumped files saved as lists:
    # x = [[sigma, mu, ksi, zero fraction], ... ]
    # plot pannel of 6 subplots comparing 3 parameters (sigma, mu, ksi) for high
    # and 3 parameters for the low flow
    def test_compare_GML_and_LM(self):
        plt.figure()

        plt.subplots_adjust(hspace = 0.3)
        the_types = ['high', 'low']
        var_names = ['$\\sigma$', '$\\mu$', '$\\xi$']
        k = 1
        for the_type in the_types:
            for i, var_name in enumerate(var_names):
                plt.subplot(2,3,k)
                plt.title(var_name + ', ' + the_type)
                xs = self.gml_parameters[the_type]
                ys = self.lm_parameters[the_type]

                xlabel = 'GML' if the_type == 'high' else 'ML'
                self.plot_scatter(xs[:,i], ys[:,i], xlabel = xlabel, ylabel = 'LM')
                k += 1

        plt.savefig('lm_vs_gml(rl_14_GML_half_and_half).png')
        pass


    #generate samples using GEV distribution and then using GML
    #try to determine the parameters of the distribution
    def test_ML(self):

        ksi0 = -0.1
        mu0 = 1000.0
        sigma0 = 100.0


        vals = []
        ps_0 = []
        ps_lm = []
        ps_gml = []
        ps_for_vals = []
        for i in range(30):
           p = random.random()
           ps_for_vals.append(p)
           vals.append(gev_cdf_inverse(p, sigma0, mu0, ksi0))

        vals = np.array(vals)

        xs = np.arange(800, 1500, 10)


        
        sigma_lm, mu_lm, ksi_lm, dummy = gevfit.optimize_stationary_for_period(vals, use_lmoments = True)
        sigma_gml, mu_gml, ksi_gml, dummy = gevfit.optimize_stationary_for_period(vals, use_lmoments = False)

        for x in xs:
            ps_0.append(gev_cdf(x, sigma0, mu0, ksi0))
            ps_lm.append(gev_cdf(x, sigma_lm, mu_lm, ksi_lm))
            ps_gml.append(gev_cdf(x, sigma_gml, mu_gml, ksi_gml))


        print 'LM: ', sigma_lm, mu_lm, ksi_lm
        print 'GML: ', sigma_gml, mu_gml, ksi_gml
        print ps_lm
        plt.title(r'$\sigma = %.2f , \mu = %.2f , \xi = %.2f$' % (sigma0, mu0, ksi0))
        plt.plot(xs, ps_0, linewidth = 2, label = 'Theory')
        plt.plot(xs, ps_lm, linewidth = 2, label = 'LM')
        plt.plot(xs, ps_gml, linewidth = 2, label = 'GML')

        plt.plot(vals, ps_for_vals, 'o',  label = 'Sample')

        plt.xlabel('$x$')
        plt.ylabel('CDF')

        plt.legend()
        plt.show()



        pass

    ##Plot as function of parameters
    ##f(x0, ksi, mu, sigma), ksi in [-0.5, 0.5],
    ## mu = 1000, sigma in [10, 10000]
    def test_unimodality_of_GEV(self):
        x0 = 1500

        mu = 1000
        data = np.array([x0])

        ksi = np.arange(-2, 2, 0.01)
        sigma = np.arange(10, 8000, 10)

        n_ksi = len(ksi)
        n_sigma = len(sigma)

        z = np.zeros((n_ksi, n_sigma))

        for i, the_ksi in enumerate(ksi):
            for j, the_sigma in enumerate(sigma):
                z[i, j] = gevfit.objective_function_stationary_high([the_sigma, mu, the_ksi], data)


        sigma, ksi = np.meshgrid(sigma, ksi)
        z = np.ma.masked_where(z == gevfit.BIG_NUM, z)
        z = np.ma.masked_where(z > 9, z)

        plt.figure()
        plt.pcolormesh(ksi, sigma, z)
        plt.colorbar()
        plt.xlabel('$\\xi$')
        plt.ylabel('$\\sigma$')
        plt.title('$\\mu = %.1f, x = %.1f$' % (mu, x0))

        plt.show()


        pass


    def runTests(self):
#        self.setUp()
        self.test_ML()
#        self.test_compare_GML_and_LM()
#        self.test_unimodality_of_GEV()


if __name__ == '__main__':
    t = GevfitTest()
    t.runTests()

