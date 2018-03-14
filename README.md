# Diffeomorphic shape modelling

This software implements a statistical shape model based on diffeomorphic transforms.

Learning shape parameters &ndash; _i.e._, a mean shape (or template) and a principal space of deformation &ndash; from a series of 2D or 3D images is cast as a maximum _a posterori_ inferrence problem in a graphical model.

This repository contains core functions and executable scripts written in Matlab.

## Content of the repository

- [<tt>backward-compatibility</tt>](backward-compatibility): A few basic Matlab functions that only appeared in recent versions of Matlab and that were reimplemented for backward compatibility.
- [<tt>core</tt>](core): Core functions. Their signature is as straightforward as possible, so they can be used in other contexts as some sort of library. They should in general work on both numerical arrays and SPM's <tt>file_array</tt> (_i.e._, memory mapped arrays).
  * [<tt>core/register</tt>](core/register) contains functions related to the registration of single subjects.
  * [<tt>core/shape</tt>](core/shape) contains functions related to population parameters learning (template, principal subspace...).
- [<tt>scripts</tt>](scripts): Executable functions that implement the complete model.
  * [<tt>pgva_model.m</tt>](scripts/pgva_model.m) implements a purely variational model based on a Riemannian PCA applied to fitted velocity fields. It can also work with known (observed) velocity fields,
  * [<tt>pgra_model.m</tt>](scripts/pgra_model.m) implements a _direct fit_ version, where all shape parameters (latent coordinates, principal subspace, residual field) are obtained by maximising the data term using Gauss-Newton optimisation.

  Both models use a set of subfunctions that have the same organisation:
  * <tt>*_input.m</tt> deals with input files (individual images, eventually some parameters of the model...),
  * <tt>*_default.m</tt> sets all default parameters,
  * <tt>*_data.m</tt> sets all working data structures,
  * <tt>*_init.m</tt> initialises the model.
  
- [<tt>utility-functions</tt>](utility-functions): Various very basic utility functions. Some of them might get moved to the [<tt>auxiliary-functions</tt> toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions) in the future.

## Dependencies

This project has strong dependencies to SPM12 and its <tt>Shoot</tt> toolbox. Both of them should be added to Matlab's path. SPM can be downloaded at [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/).

Core functions also depend on our [<tt>auxiliary-functions</tt> toolbox](https://github.com/WTCN-computational-anatomy-group/auxiliary-functions), which gathers lots of low-level functions.

Furthermore, executable scripts depend on our [<tt>distributed-computing</tt> toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing), which helps parallelising parts of the code either on the local workstation (using Matlab's parallel processing toolbox) or on a computing cluster (see the toolbox help file for use cases and limitations).

Note that if these toolboxes are all located in the same folder, _i.e._:
* <tt>./shape-toolbox</tt>
* <tt>./auxiliary-functions</tt>
* <tt>./distributed-computing</tt>

a call to [<tt>setpath.m</tt>](setpath.m) adds all necessary folders to Matlab's path.
