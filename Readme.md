### Testing CEC benchmark 2022 Python version of code 

Testing CEC benchmark 2022 from https://github.com/P-N-Suganthan/2022-SO-BO

The codes could be found at https://github.com/P-N-Suganthan/2022-SO-BO/blob/main/Python-CEC2022.zip

The raw codes were taken there and modified for testing purposes.

### Running the tests

To run tests check combine-dim-func.py where you set number of functions i, dimensions j, and number of element k.

Test run is paralel for all the functions.

For manual testing check run-tests.py. It is possilbe to run this script with arguments number of function, dimension and number of element.

### Making graphs and tables

For graphs and comparisons check C version of the code at https://github.com/FrostNiles/cec-2022-testing-C

### Results

Results are find in test_data/result

### portability

It is possilbe to run all the tests in different years, just check the years documentation, change mains so it writes results in appropriate files.

Change the precission - it is now to 10e-08 in every document, it is possilbe to change it globally one day not now. (seven digits now in scripts.py and 9 floating point numbers in mains)

First try run-tests.py for single element few times.

Lastly change the combine-dim-func.py so you can run it automatically for the years functions.