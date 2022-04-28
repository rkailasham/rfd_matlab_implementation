# rfd_matlab_implementation
MATLAB scripts to implement the random finite difference (RFD) method for the calculation of divergence

"iv_compare_rfd_diffmat_divergence.m" compares the analytically known divergence of diffusion tensor for freely-draining chains against that calculated numerically using the RFD method.

For Rouse chains with hydrodynamic interactions (HI) only (no IV), the divergence of the diffusion tensor is zero. 

"hi_compare_rfd_diffmat_divergence.m" compares the numerically calculated divergence for chains with HI against this analytical result.
