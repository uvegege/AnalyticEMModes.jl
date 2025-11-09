
#using Pkg
#cd("C:\\MisProyecto\\Upload/AnalyticEMModes.jl//")
#Pkg.activate(".")


module AnalyticEMModes

    using Bessels
    using MathieuF

    export wavenumber, phase_constant, attenuation_factor, propagation_constant
    export cutoff_frequency, mode_wavelength, mode_impedance
    
    export kc_rwg, te_rwg_fields, tm_rwg_fields, te_normalization_rwg, tm_normalization_rwg

    export kc_cwg, te_cwg_fields, tm_cwg_fields 
    export te_normalization_cwg, tm_normalization_cwg, metric_and_unit_cylindrical
    
    export kc_coax, te_coax_fields, tm_coax_fields, tem_coax_fields, te_normalization_coax
    export tm_normalization_coax, characteristic_coax_equation_te, characteristic_coax_equation_tm
    
    export kc_radial, phase_constant_radial, cutoff_frequency_radial
    export te_radial_fields, tm_radial_fields, te_zw, tm_zw, te_normalization_radial, tm_normalization_radial
    
    export kc_wedge, te_wedge_fields, tm_wedge_fields, te_normalization_wedge, tm_normalization_wedge
    
    export kc_ewg, te_ewg_fields, tm_ewg_fields 
    export cart2elliptic, metric_and_unit_elliptic, ce_m, se_m, Ce_m, Se_m
    
    export first_n_modes_rwg, first_n_modes_cwg, first_n_modes_coax, first_n_modes_radial, first_n_modes_ewg

    include("./common.jl")
    
    # Rectangular Coordinates
    include("./rectangularwg.jl") # Rectangular Waveguide

    # Cylindrical Coordinates:
    include("./circularwg.jl") # Circular Waveguide
    include("./coaxialwg.jl") # Coaxial Waveguide
    include("./radialwg.jl") # Radial Waveguide
    include("./wedgewg.jl") # Wedge Waveguide

    # Elliptic Coordinates
    include("./ellipticalwg.jl") 
    
    # Spherical Coordinates
    #include("./Spherical.jl") # TODO

    # Sort modes
    include("./sortmodes.jl")

end

using .AnalyticEMModes

