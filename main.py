# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 17:14:04 2022

@author: Abhishek Kumar
"""

import numpy as np
from CEC2022 import cec2022_func
import sys

#np.set_printoptions(precision=15)

#if number of arguments is less than 2, then use default values

mx = 1 #number of runs
if len(sys.argv) < 3:
    
    nx = 10 #dimensions

    fx_n = 1 #function number
    argNum = 1 #argument number
else:
    nx = int(sys.argv[2])
    fx_n = int(sys.argv[1])
    argNum = int(sys.argv[1])
with open(f'test_data/shift_data_{argNum}.txt', 'r') as file:
    data = file.read()


dataNumbers = [float(i) for i in data.split()]
dataNumbers = np.array(dataNumbers)


test = dataNumbers.reshape(-1, 1)
#and keep there only the first nx values
test = test[:nx]

CEC = cec2022_func(func_num = fx_n)

x = test
F = CEC.values(x)
for num in F.ObjFunc:
    print(f"{num:.9f}")
    with open(f"test_data/current_result_{argNum}.txt", 'w') as file:
        file.write(f"{num:.9f}")