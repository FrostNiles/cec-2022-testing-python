import subprocess
import sys

func_num = sys.argv[1]
dimension = sys.argv[2]
# Check if dimension is valid
if dimension not in ['2', '10','20', '30', '50', '100']:
    print("Invalid dimension. Please choose from 2, 10, 20, 30, 50, or 100.")
    sys.exit(1)

# run the main.py file with the specified arguments
args = ['python', 'main.py', func_num, dimension]

# Run the compiled file with the specified arguments
subprocess.run(args, check=True)