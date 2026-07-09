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

## Triangular waveguides

```@docs
EquilateralTriangle
RightIsoscelesTriangle
HalfEquilateralTriangle
ReflectionMode
kc_equilateral
kc_right_isosceles
kc_half_equilateral
first_n_modes_equilateral
first_n_modes_right_isosceles
first_n_modes_half_equilateral
te_equilateral_fields
tm_equilateral_fields
te_right_isosceles_fields
tm_right_isosceles_fields
te_half_equilateral_fields
tm_half_equilateral_fields
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
AnalyticEMModes.kc_sph
AnalyticEMModes.te_sph_fields_lmax
AnalyticEMModes.tm_sph_fields_lmax
AnalyticEMModes.te_sph_fields
AnalyticEMModes.tm_sph_fields
AnalyticEMModes.te_normalization_sph
AnalyticEMModes.tm_normalization_sph
AnalyticEMModes.first_n_modes_sph
AnalyticEMModes.spherical_to_cartesian_fields
AnalyticEMModes.to_svector
AnalyticEMModes.mn_vectors_sph
AnalyticEMModes.mn_sph_vectors_lmax
AnalyticEMModes.m_sph_vectors_lmax
AnalyticEMModes.n_sph_vectors_lmax
AnalyticEMModes.m_normalization_sph
AnalyticEMModes.n_normalization_sph
AnalyticEMModes.te_from_mn_sph
AnalyticEMModes.tm_from_mn_sph
```

## Spheroidal modes

```@docs
AnalyticEMModes.kc_spheroidal
AnalyticEMModes.spheroidal_families
AnalyticEMModes.mn_spheroidal_vector
AnalyticEMModes.m_spheroidal_vector
AnalyticEMModes.n_spheroidal_vector
AnalyticEMModes.mn_spheroidal_vectors
AnalyticEMModes.m_spheroidal_vectors
AnalyticEMModes.n_spheroidal_vectors
AnalyticEMModes.mn_spheroidal_vectors_mnmax
AnalyticEMModes.ProlateSpheroidalBasis
AnalyticEMModes.OblateSpheroidalBasis
AnalyticEMModes.SpheroidalB
AnalyticEMModes.obl2cart
AnalyticEMModes.pro2cart
AnalyticEMModes.cart2pro
AnalyticEMModes.cart2obl
AnalyticEMModes.spheroidal_parameter
AnalyticEMModes.scale_factors_prolate
AnalyticEMModes.scale_factors_oblate
AnalyticEMModes.prolate_vector_to_cartesian
AnalyticEMModes.oblate_vector_to_cartesian
```

## Sort Modes
```@docs
first_n_modes_rwg
first_n_modes_cwg
first_n_modes_coax
first_n_modes_radial
first_n_modes_ewg
```
