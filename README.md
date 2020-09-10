# generative_paper
Code for model described in "Tractable generative network models provide insights into the effect of trade and partnership dynamics on endemic livestock disease".

The code contained within allows for users to run stochastic simulations of the model described in the above article.
Code is written in C++ and we use the standard Gillespie Stochastic Simulation Algorithm (Gillespie SSA).

Code requires installation of GSL random number generator. This can be obtained from https://www.gnu.org/software/gsl/ or code can be edited to use user desired random number generator.

Animal movement data (CTS dataset) used in the article is not included as this contains personal information of individuals and is not available for open access.
