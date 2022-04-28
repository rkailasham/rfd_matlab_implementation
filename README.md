# rfd_matlab_implementation
MATLAB scripts to implement the random finite difference (RFD) method [1,2] for the calculation of divergence.

28-APR-2022:

"iv_compare_rfd_diffmat_divergence.m" compares the analytically known divergence of diffusion tensor for freely-draining chains against that calculated numerically using the RFD method. The analytical calculation of the divergence is based on recursive relations, as explained in the "Supplementary Material" section of the paper titled "Rouse model with fluctuating internal friction" [J. Rheol. 65, 903 (2021)]. The MATLAB codes needed for this purpose are all included in this directory, and also in the other repository "div_calc_continued_fractions".

For Rouse chains with hydrodynamic interactions (HI) only (no IV), the divergence of the diffusion tensor is zero. This result follows from the incompressibility of the solvent, and the manner in which the Rotne-Prager-Yamakawa (RPY) tensor is constructed.

"hi_compare_rfd_diffmat_divergence.m" compares the numerically calculated divergence for chains with HI against this analytical result.

The RFD method is used in the evaluation of the stress tensor for bead-spring-dashpot chains with fluctuating internal friction and hydrodynamic interactions, as described in detail in the preprint "Shear viscosity for finitely extensible chains with fluctuating internal friction and hydrodynamic interactions" [3], a copy of which is saved in this repository as "fene_iv_hi_arxiv_paper.pdf".


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

References

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[1] B. Sprinkle, F. Balboa Usabiaga, N. A. Patankar, and A. Donev, J. Chem. Phys. 147, 244103 (2017).

[2] B. Sprinkle, A. Donev, A. P. S. Bhalla, and N. Patankar, J. Chem. Phys. 150 (2019).

[3] R. Kailasham, Rajarshi Chakrabarti, and J. Ravi Prakash, arXiv:2204.10656
