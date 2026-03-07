# API Reference

```@contents
Pages = ["api_reference.md"]
Depth = 2:2
```

## Index
```@index
Pages = ["api_reference.md"]
```

## General functions

```@docs
wavenumber
phase_constant
attenuation_factor
propagation_constant
cutoff_frequency
mode_wavelength
mode_impedance
```

## Rectangular waveguide

```@docs
kc_rwg
te_rwg_fields
tm_rwg_fields
te_normalization_rwg
tm_normalization_rwg
```


## Circular waveguide

```@docs
kc_cwg
te_cwg_fields
tm_cwg_fields 
te_normalization_cwg
tm_normalization_cwg
metric_and_unit_cylindrical
```


## Coaxial waveguide

```@docs
kc_coax
te_coax_fields
tm_coax_fields
tem_coax_fields
te_normalization_coax
tm_normalization_coax
characteristic_coax_equation_te
characteristic_coax_equation_tm
```


## Radial waveguide

```@docs
kc_radial
phase_constant_radial
cutoff_frequency_radial
te_radial_fields
tm_radial_fields
te_zw
tm_zw
te_normalization_radial
tm_normalization_radial
```


## Wedge waveguide

```@docs
kc_wedge
te_wedge_fields
tm_wedge_fields
te_normalization_wedge
tm_normalization_wedge
```


## Elliptical waveguide

```@docs
kc_ewg
te_ewg_fields
tm_ewg_fields 
cart2elliptic
metric_and_unit_elliptic
ce_m
se_m
Ce_m
Se_m
```

## Spherical modes

```@docs
kc_sph
te_sph_fields_lmax
tm_sph_fields_lmax
te_sph_fields
tm_sph_fields
te_normalization_sph
tm_normalization_sph
first_n_modes_sph
spherical_to_cartesian_fields
to_svector
mn_vectors_sph
mn_sph_vectors_lmax
m_sph_vectors_lmax
n_sph_vectors_lmax
m_normalization_sph
n_normalization_sph
te_from_mn_sph
tm_from_mn_sph
```

## Spheroidal modes

```@docs
kc_spheroidal
spheroidal_families
mn_spheroidal_vector
m_spheroidal_vector
n_spheroidal_vector
mn_spheroidal_vectors
m_spheroidal_vectors
n_spheroidal_vectors
mn_spheroidal_vectors_mnmax
ProlateSpheroidalBasis
OblateSpheroidalBasis
SpheroidalB
obl2cart
pro2cart
cart2pro
cart2obl
spheroidal_parameter
scale_factors_prolate
scale_factors_oblate
prolate_vector_to_cartesian
oblate_vector_to_cartesian
```

## Sort Modes
```@docs
first_n_modes_rwg
first_n_modes_cwg
first_n_modes_coax
first_n_modes_radial
first_n_modes_ewg
```
