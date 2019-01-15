# PDMCSampler

This repository deals with 3 different non-reversible continuous Markov Chain Monte Carlo based on piecewise deterministic Markov processes (see *Davis, Markov Models and Optimization*). The source folder contains:
* a subfolder *Models* containing the target distributions we want to sample from (currently only Gaussian)
* *PDSamplers.jl* and *kernel.jl* contains the implementation of the ZigZag Sampler, Bouncy Particle Sampler and Coordinate Sampler. Kernels implement the change of the velocity coordinates when the inhomogeneous Poisson event occours.
* *Simulations.jl* calls and plots the different algorithms for comparisons  
* *Faber1.jl* implements the sparse matrix we find when we want to run the ZigZag for sampling a linear SDE (See *Notes.pdf*) for explanations.
* Finally *Experiment.jl* calls the ZigZag inside the coordinates of the basis expansion of the Bridge.  

## Some Plots



## Usage
Once you import the folder Source and install the package *Ciesielski*, you can run both *Experiment.jl* and *Simulations.jl*. You can play with the parameters of the Multivariate Gaussian RV and with the coefficients of the Linear sde
$$dX_t = (a + bX_t)dt + dW_t$$

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[MIT](https://choosealicense.com/licenses/mit/)
