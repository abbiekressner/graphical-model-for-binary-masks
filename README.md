# Graphical model for binary masks (GMBM)
This toolbox provides a set of functions that implement the graphical model for binary masks presented in the following paper:

```
Kressner, A. A. and Rozell, C. J. (2015). “Structure in time- frequency binary masking errors and its impact on speech intelligibility”, The Journal of the Acoustical Society of America 137, 2025–2035.
```

The graphical model is a useful tool for 

- *Sampling*: randomly generating binary masks that contain specific amounts of errors and structure
- *Parameter estimation*: learning about the errors that estimation algorithms tend to make

This toolbox relies on the [Undirected Graphical Models toolbox](http://www.cs.ubc.ca/~schmidtm/Software/UGM.html) by Mark Schmidt. To ensure cohesion between versions of the toolbox, the UGM toolbox is wrapped into this one.

## Installation
Download (using Git, SVN, or a zip file), open Matlab, and run `addpath(genpath('path/to/repo/'));`.

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

```
Kressner, A. A. and Rozell, C. J. (2015). “Structure in time- frequency binary masking errors and its impact on speech intelligibility”, The Journal of the Acoustical Society of America 137, 2025–2035.
```
