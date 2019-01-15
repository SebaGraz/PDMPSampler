# PDMCSampler

This repository deals with 3 different non-reversible continuous Markov Chain Monte Carlos based on piecewise deterministic Markov processes (see *Davis, Markov Models and Optimization*). The source folder contains:
* a subfolder *Models* containing the target distributions we want to sample from (currently only Gaussian)
* *PDSamplers.jl* and *kernel.jl* contain the implementation of the ZigZag Sampler, Bouncy Particle Sampler and Coordinate Sampler. *Kernels.jl* implements the change of the velocity coordinates when the inhomogeneous Poisson event occours.
* *Simulations.jl* calls and plots the different algorithms for comparisons  
* *Faber1.jl* implements the sparse matrix we need to create when we want to run the ZigZag for sampling a linear SDE (See *Notes.pdf*) for explanations.
* Finally *Experiment.jl* calls the ZigZag inside the coordinates of the basis expansion of the Bridge.  

## Some Plots
Plots coming from the Script *Simulations.jl*: Comparison of the three algorithms: respectively Bouncy, ZigZag and Coordinate Sampler with the same internal clock fixed to 1000 units. The target distribution is 2d-Gaussian with zero mean and correlation equal to -0.5.
![Bouncy](/images/ciao.png)
![ZigZag](/images/ZigZag.png)
![Coordinate](/images/Coordinate.png)
Plots coming from the Script *Experiment.jl*: Final Bridge after 21 seconds of internal clock of ZigZag algorithm and samples every 0.5 seconds from time 1 to time 21 of internal clock (mixing of the algorithm). The target measure is given by the sde $$dX_t = (a + bX_t)dt + dW_t$$ with parameter a = 10.0 ; b = -1.0 ; the number of level before truncation in the Faber Schauder expansion is L = 10 implying 2047 dimensions.  
![Coordinate](/images/Bridge.png)
![Coordinate](/images/Mixing.png)
## Usage
Once you import the folder Source and install the package *Ciesielski*, you can run both *Experiment.jl* and *Simulations.jl*. You can play with the parameters of the Multivariate Gaussian RV and with the coefficients of the Linear sde
$$dX_t = (a + bX_t)dt + dW_t$$. It s reccommended

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.


## License
[MIT](https://choosealicense.com/licenses/mit/)
