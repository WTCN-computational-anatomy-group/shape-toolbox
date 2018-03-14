# Diffeomorphic shape modelling

This software implements a statistical shape model based on diffeomorphic transforms.

Learning shape parameters &ndash; _i.e._, a mean shape (or template) and a principal space of deformation &ndash; from a series of 2D or 3D images is cast as a maximum _a posterori_ inferrence problem in a graphical model.

This repository contains core functions and executable scripts written in Matlab.

## Standalone release

## The model

This work relies on a generative model of shape in which individual images (of brains, in particular) are assumed to be generated from a mean shape &ndash; commonly named template &ndash; deformed according to a transformation of the coordinate space. Here, these transformations are diffeomorphisms, _i.e._, one-to-one invertible mappings that allow for very large deformations. By using the _geodesic shooting_ framework, we parameterise these transformations by their _initial velocity_, which can be seen as an infinitesimal (very small) deformation. The _a posteriori_ covariance structure of these velocity fields is infered by making use of a technique related to the well-known _principal component analysis_, adapted to the particular structure of the space on which lie velocity fields, called a Riemannian manifold. Our model also includes a rigid-body transform, whose role is to factor out all deformations induces by brains misalignment.

All considered, the following variables are infered:
- <tt>W = [w1 .. wK]</tt>: the principal subspace of deformation, made of K _principal geodesics_ ;
- <tt>z</tt>: transformation coordinates in the principal subspace, which is a low-dimensional representation of each subject in terms of deformation of the template ;
- <tt>A</tt>: precision matrix (_i.e._, inverse covariance) of the latent coordinates. At the optimum, it should be a diagonal matrix that contains the variance along each principal component, or in other words, their scale ;
- <tt>v</tt>: the velocity field of each subject. It is only an explicit random variable in the PGVA model, in which case the residual field, <tt>r = v - Wz</tt>, can be recovered by substracting the principal representation ;
- <tt>r</tt>: alternatively, the residual field can be explicitely infered, as is the case in the PGRA model. Then, the initial velocity is reconstructed according to <tt>v = Wz + r</tt> ;
- <tt>lam</tt>: precision of the residual field, also named _anatomical noise_ ;
- <tt>q</tt>: parameters of the rigid-body transform. Note that there are options to use different kind of affine transforms instead, however it is not advised, as differences in size should be captured by the shape model.

The following parameters are manually set and impact the model's behaviour:
- <tt>A0</tt> and <tt>n0</tt>: prior expected value of the latent precision matrix and its degrees of freedom, which should be seen as the virtual number of subjects that weight this prior belief ;
- <tt>l0</tt> and <tt>n0</tt>: prior expected value of the residual precision matrix and its degrees of freedom, which should be seen as the virtual number of subjects that weight this prior belief ;

### PGRA vs PGVA


## User documentation

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

## Contributors

This software was developed under the [_Human Brain Project_](https://www.humanbrainproject.eu) (SP2) flagship by John Ashburner's [Computational Anatomy Group](http://www.fil.ion.ucl.ac.uk/Ashburner/) at the [Wellcome Centre for Human Neuroimaging](http://www.fil.ion.ucl.ac.uk/) in UCL.
- The shape toolbox was mainly developed by Yaël Balbastre with invaluable help from John Ashburner
- The <tt>auxiliary-functions</tt> and <tt>distributed-computing</tt> toolboxes were developed by Mikael Brudfors and Yaël Balbastre

## License

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.
