
# Windows

Used git bash for doing all of this.

## Get Python installed, and make/activate a virtual environment

1. Install Python, if you need to. See https://www.python.org/downloads/. note that the global planner code is currently tested with Python 3.5.4.
2. Install virtual environment if needed, `pip install virtualenv`
3. Make a virtual environment for installing all the right packages to run the code. Create this virtual environment in a convenient directory. E.g. from that directory: `virtualenv --python=/c/Users/STARLab/AppData/Local/Programs/Python/Python35/python.exe venv` ( where the path is to the Python 3.5 executable)
4. Activate the virtual environment: `source venv/Scripts/activate`

## Clone the repo

In git bash, In the desired location to locate the repo, run `git clone --recurse-submodules git@github.mit.edu:star-lab/GlobalPlanner.git`. The recursive is needed to get the circinus_tools submodule as well.

## Install the requirements

In GlobalPlanner/python_runner, run `pip install -r requirements.txt`

## Install python matlab package.

Follow the instructions here: http://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html?refresh=true
1. Run git bash as an administrator (Find the program, right-click, run as administrator)
2. You'll cd to somewhere like:  "/c/Program Files/MATLAB/R2017a/extern/engines/python"
3. In the desired python (virtual) environment, run `python setup.py install`
4. Now the "matlab" package should be available when running Python in this environment

## Install optimization solver backend

The wonderful folks over at Gurobi make academic license licenses easy to get, so I recommend installing Gurobi for use as the optimization backend. It's a state-of-the-art commercial software, and on par with CPLEX for performance in running the global planner in my experience. In general, the commercial solvers are much faster than open source solvers. However, if you need to use an open source solver, the GNU Linear Programming Kit is a good place to start (https://www.gnu.org/software/glpk/).

### Gurobi

I was able to download the software and license for Gurobi here: http://www.gurobi.com/academia/for-universities. The process is straightforward, and doesn't require any additional explanation here. After installation, you end up with a commandline tool gurobi_cl.exe which you can run yourself, and also is used by the Pyomo optimization framework used in the global planner.

The code has been extensively tested with version 7.5.2 of Gurobi, and a bit with version 8.0.0

### CPLEX

Note that CPLEX restricts academic installation to only one machine. Nonetheless it's a great tool as well.

The code has been tested with IBM ILOG CPLEX v12.8 (https://www.ibm.com/products/ilog-cplex-optimization-studio). 

It can be downloaded free for students (this link worked for me: https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=733c3d21-0ce1-e711-80fa-000d3af41938&pmv=00000000-0000-0000-0000-000000000000)

