import newton_rule
from models import ReactorModels
import unittest

class OverallCaseTest (unittest.TestCase):
    """ Checking if the entire functions in each module will works """


    def test_numerical_analysis(self):
        """ Can it run numerical analysis for the sample praam data? """
        sample_param = [0,3,5,5,4,2,1,0]
        soln = newton_rule.numerical_rules(sample_param, 5)
        self.assertEqual(soln, 100.0)

    def test_tank_model(self,):
        """ Run the conversion prediction for tank in series model with the values from the sample data"""
        conc = [0,3,5,5,4,2,1,0]  
        conversion = ReactorModels.tank_in_series(14.667, 48.212, 0.3, conc, 1)
        self.assertEqual(conversion, 0.953)

    def test_dispersion_model(self):
        """ Run the conversion predication for dispersion model"""
        conc = [0,3,5,5,4,2,1,0]  
        conversion = ReactorModels.dispersion(14.667, 48.212, 0.3, conc, 1)
        self.assertEqual(conversion, 0.96)

    def test_idealcstr_model(self):
        """ Run the conversion prediction for ideal continuous stirred tank reactor model"""
        conc = [0,3,5,5,4,2,1,0]  
        conversion = ReactorModels.ideal_cstr_model(14.667, 0.3, conc, 1)
        self.assertEqual(conversion, 0.815)

    def test_idealpfr_model(self):
        """ Run the conversion prediction for ideal plug flow reactor model"""
        conc = [0,3,5,5,4,2,1,0]  
        conversion = ReactorModels.ideal_pfr_model(14.667, 0.3, conc, 1)
        self.assertEqual(conversion, 0.988)

if __name__ == '__main__':
    unittest.main()