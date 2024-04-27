# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 17:14:04 2022

@author: Abhishek Kumar
"""

import numpy as np
from myCEC2022 import cec22_test_func
import os
import re
from myCEC2022 import ini_flag

def main(func_num, n):
    m = 2
    x = np.zeros(m * n)
    f = np.zeros(m)

    for i in range(func_num, func_num + 1):
        file_name = f"test_data/shift_data_{func_num}.txt"
        if not os.path.exists(file_name):
            print("\n Error: Cannot open input file for reading \n")
            return

        with open(file_name, 'r')as fpt:
            data = fpt.read()

        numbers = re.findall(r'(-?\d+\.\d+e[+-]\d+)', data)
        for j in range(n):
            x[j] = float(numbers[j])


        x[n:] = 0.0

        for k in range(1):
            cec22_test_func(x, f, n, m, func_num)
            for j in range(1):
                print(f" f{func_num}(x[{j+1}]) = {f[j]:.10f},")

                result_file = f"test_data/current_result_{func_num}.txt"
                with open(result_file, 'w+') as rst:
                    rst.write(f" f{func_num}(x[{j+1}]) = {f[j]:.10f},")
                print()

main(1, 10)

