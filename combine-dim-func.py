import runpy
import re
import sys

skipped = {6,7,8,12}
for i in range(6, 13):
    if i in skipped:
        continue
    # CEC 2022 has 12 functions
    for j in [10, 20]:
        for k in range(1, j+1):
            args = {
                'arg1': str(i),
                'arg2': str(j),
                'arg3': str(k)
            }
            
            runpy.run_path('./run-tests.py', init_globals=args)


