#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Sobol sensitivity analysis

Steady-state functions (bone remodeling):
	- Bone_xC, Bone_xB, Bone_xT, Bone_xW

Steady-state functions (bone metastasis):
	- Cancer_xC, Cancer_xB, Cancer_xT, Cancer_xW
"""

#import sys
#sys.path.append('../..')


from SALib.analyze import sobol
from SALib.sample import saltelli
from SALib.util import read_param_file

# steady-state func
import Cancer_xW 

# bone_problem.txt or metastasis_problem.txt
problem = read_param_file('metastasis_problem.txt')

# Execute Sobol
NumberOfSamples = 400000
param_values = saltelli.sample(problem, NumberOfSamples, calc_second_order=True)
Y = Cancer_xW.evaluate(param_values)
Si = sobol.analyze(problem, Y, calc_second_order=True, conf_level=0.95, print_to_console=True)
