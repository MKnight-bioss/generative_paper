# generative_paper
Code for model described in "Tractable generative network models provide insights into the effect of trade and partnership dynamics on endemic livestock disease".

The code contained in main.cpp allows for users to run stochastic simulations of the model described in the above article. Specifically we use the parameterisation for the homogeneous system described in the Supplementary Material that can be run without the use of data.
Code is written in C++ and we use the standard Gillespie Stochastic Simulation Algorithm (Gillespie SSA).

Code requires installation of GSL random number generator. This can be obtained from https://www.gnu.org/software/gsl/ or code can be edited to use user desired random number generator. Users can replace any calls to the gsl rng with their desired rng by removing code where necessary.

Animal movement data (CTS dataset) used in the article is not publicly available as it contains confidential information of individuals, e.g. names, addresses.
