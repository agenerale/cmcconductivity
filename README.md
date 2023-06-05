**Dataset description:**  
This is the repository for the Materials and Design paper titled ["Reduced-order Models for Microstructure-Sensitive Effective Thermal Conductivity of Woven Ceramic Matrix Composites with Residual Porosity"](https://www.sciencedirect.com/science/article/pii/S0263822321008618)

**Files Included:**
* k11.mat, k22.mat, k33.mat - ABAQUS effective thermal conductivity values (k11.mat, k22.mat, k33.mat)
* pcs.mat - The low-dimensional model input of principal component scores for all microstructures. Statistical quantification of microstructures can be accomplished through the PyMKS github page: https://github.com/materialsinnovation/pymks. 
* avg_hflux.py - post-processing script for extracting effective orthotropic thermal conductivity
* analyticalK1.m, analyticalK3.m - implementation of the Hierarchical Two Layer model from ["Modeling the Transverse Thermal Conductivity of 2-D SiC/SiC Composites Made with Woven Fabric"](https://www.tandfonline.com/doi/abs/10.13182/FST04-A533).
* main_gpr_linear.m - main script for building of predictive Gaussian processes for each orthotropic component.

Further information regarding the GPR model implementation can be found in the Mathwords page: https://www.mathworks.com/help/stats/fitrgp.html
