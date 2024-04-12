Run driver.m.
Function driver.m will call the optimizer LevenbergMarquardt.m.
The objective function F(r) = 0.5||r||^2 is constructed as the mean squared error. 
To generate it, LevenbergMarquardt.m calls res_and_Jac.m that returns the vector-function r and its Jacobian matrix J. 
The vector r is computed in res.m as in equation (10) in NeuralNetworks4PDEs.pdf. 
Function res.m calls function setup.m where the functions involved into the solution model and their derivatives are defined.
The computed and the exact solution are evaluated at training points in function evaluateNNsolution.m

In order to solve another PDE, you need to modify res.m, setup.m, and evaluateNNsolution.m.