Response of a compressible viscoleastic shell to the tidal potential

Author: Vít Beran, Charles University Prague

Compilation: gfortran -O3 main.f90 mCoef.f90 mEqs.f90 mFunc.f90 mPar.f90 sPolint.f90 sBan* 

Spherical harmonics implementation of the viscoelastic deformation with Maxwell rheology

Load: Tidal periodic / Heaviside load supported

Use cases: Compute tidal Love numbers, tidal heating due to oblique tides,... 

For the compressible models, the density profile must satisfy the Adams Williamson equation (results in divergence sooner or later otherwise) 
