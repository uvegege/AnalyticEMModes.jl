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

## Sort Modes
```@docs
first_n_modes_rwg
first_n_modes_cwg
first_n_modes_coax
first_n_modes_radial
first_n_modes_ewg
```
