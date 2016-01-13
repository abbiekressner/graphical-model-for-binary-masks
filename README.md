# Graphical model for binary masks (GMBM)
This toolbox provides a set of functions that implement the graphical model for binary masks presented in the following paper:

> Kressner, A. A. and Rozell, C. J. (2015). “Structure in time-frequency binary
> masking errors and its impact on speech intelligibility”, The Journal of the
> Acoustical Society of America 137, 2025–2035.

The graphical model is a useful tool for 

- *Sampling*: randomly generating binary masks that contain specific amounts of errors and structure
- *Parameter estimation*: learning about the errors that estimation algorithms tend to make

This toolbox relies on the [Undirected Graphical Models toolbox](http://www.cs.ubc.ca/~schmidtm/Software/UGM.html) by Mark Schmidt. To ensure cohesion between versions of the toolbox, the UGM toolbox is wrapped into this one.

## Installation
1. Download the toolbox using Git, SVN, or a zip file.
2. Open Matlab, and run `addpath(genpath('path/to/repo/'));`, where you replace `'path/to/repo/'` with the actual path to the toolbox on your own computer (e.g., `addpath(genpath('/Users/abbie/graphical-model-for-binary-masks/'));`).

## Getting started
### Sampling
To do sampling, use the function `generate_mask` as follows to create a model-generated mask:
```
load('example_masks','mask_ideal');  
params.alpha = 0.2;  
params.beta = 0.1;  
params.gamma = 2.0;  
mask = generate_mask(mask_ideal{1},params);  
```
This minimal working example uses the default for the `verbose`, `accuracy`, and `max_iter` settings. It also loads the default `example_lookup_table` for converting the `alpha` and `beta` to the weights `A` and `B` (for an explanation of these variables, see Kressner and Rozell, JASA, 2015, pg 2027). The mapping between `alpha` and `A` and `beta` and `B` is, in general, dependent on the speech corpus, noise type, signal-to-noise ratio (SNR), and clustering parameter `gamma`. Thus, to generate masks most efficiently, it is best to build a `lookup_table` specifically for each set of corpus, noise, SNR, local criterion threshold, filterbank, etc. and `gamma` combination using the included `build_lookup_table` function:
```
eg_ideal = load(); % TODO load an IBM that you create with your own code  
g = 1.0:0.5:2.5; % specify the gamma you are interested in using  
lookup_table_using_my_settings = build_lookup_table(eg_ideal,g);  
save('lookup_table_using_my_settings.mat','lookup_table_using_my_settings');  
```
Then call `generate_mask` with this new `lookup_table`:
```
mask = generate_mask(eg_ideal,params,true,0.01,50,lookup_table_using_my_settings);
```

### Parameter estimation
To do parameter estimation, use the function `learn_parameters`.
```
load('example_masks');  
gamma = learn_parameters(mask_ideal,mask_gmm);  
```
This minimal working example uses the default for the `verbose` and `optimization` settings. 

## License
Written by Abigail Kressner in 2014.

This file is part of GMBM.

GMBM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GMBM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GMBM.  If not, see <http://www.gnu.org/licenses/>.

## Citation
If you use this toolbox, please cite the following paper:

> Kressner, A. A. and Rozell, C. J. (2015). “Structure in time-frequency binary
> masking errors and its impact on speech intelligibility”, The Journal of the
> Acoustical Society of America 137, 2025–2035.
