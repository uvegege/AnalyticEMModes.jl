"""
    ModeSort{T1, T2, T3}

Structure to represent and sort electromagnetic waveguide modes by cutoff frequency.

# Fields
- `kc::T1`: Cutoff wavenumber
- `kind::T2`: Mode type symbol
- `m::T3`: First mode index (typically azimuthal)
- `n::T3`: Second mode index (typically radial or vertical)

# Mode type symbols
## Rectangular waveguides
- `:TE`, `:TM`: Standard TE and TM modes

## Circular and coaxial waveguides
- `:TE`, `:TM`: Azimuthally symmetric modes (m=0)
- `:TEc`, `:TEs`: TE modes with cosine/sine angular variation (cos(mθ)/sin(mθ))
- `:TMc`, `:TMs`: TM modes with cosine/sine angular variation (cos(mθ)/sin(mθ))
- `:TEM`: TEM mode (coaxial only, kc=0)

**Note**: For m > 0, the `c` (cosine) and `s` (sine) modes are degenerate polarizations
with **identical cutoff frequencies**.

## Elliptic waveguides
- `:TEe`, `:TMe`: Even modes (using even Mathieu functions, cosine-elliptic)
- `:TEo`, `:TMo`: Odd modes (using odd Mathieu functions, sine-elliptic)

**Note**: For elliptic waveguides, even and odd modes have **different cutoff frequencies**
(non-degenerate).

## Radial and wedge waveguides
- `:TEM`: TEM mode (radial, kc=0)
- `:TE`, `:TM`: Higher-order modes
"""
struct ModeSort{T1, T2, T3}
    kc::T1
    kind::T2  # :TM, :TE, :TEM, :TMe, :TMo, :TEc, :TEs, :TMc, :TMs, :TEe, :TEo
    m::T3
    n::T3
end

Base.isless(a::ModeSort, b::ModeSort) = a.kc < b.kc

function push_heap!(heap::Vector{T}, visited, val::T) where T
    key = (val.kind, val.m, val.n)
    #if key in visited
    #    return 
    #end
    push!(heap, val)
    push!(visited, key)
    i = length(heap)
    while i > 1
        parent = div(i, 2)
        if heap[i] < heap[parent]
            heap[i], heap[parent] = heap[parent], heap[i]
            i = parent
        else
            break
        end
    end
end

function push_heap!(heap::Vector{T}, val::T) where T
    push!(heap, val)
    i = length(heap)
    while i > 1
        parent = div(i, 2)
        if heap[i] < heap[parent]
            heap[i], heap[parent] = heap[parent], heap[i]
            i = parent
        else
            break
        end
    end
end

function pop_heap!(heap::Vector{T}) where T
    if isempty(heap)
        error("Heap is empty")
    end
    minval = heap[1]
    heap[1] = heap[end]
    pop!(heap)
    i = 1
    while true
        left = 2*i
        right = 2*i + 1
        smallest = i
        if left <= length(heap) && heap[left] < heap[smallest]
            smallest = left
        end
        if right <= length(heap) && heap[right] < heap[smallest]
            smallest = right
        end
        if smallest != i
            heap[i], heap[smallest] = heap[smallest], heap[i]
            i = smallest
        else
            break
        end
    end
    return minval
end


function addmode_rectangular!(heap, visited ,dims, m, n, T; m_max = nothing)
    T == :TEM && return nothing
    m == 0 && n == 0 && return nothing

    key = (m, n)
    push!(visited, key)

    kc = kc_rwg(dims[1], dims[2], m, n)
    push_heap!(heap, ModeSort(kc, :TE, m, n))
    if m >= 1 && n >= 1
        push_heap!(heap, ModeSort(kc, :TM, m, n))
    end

    return nothing
end

function addmode_circular!(heap, visited ,dims, m, n, T; m_max = nothing)
    T == :TEM && return nothing
    n <= 0 && return nothing

    key = (m, n)
    push!(visited, key)

    if m == 0
        kc = kc_cwg(dims[1], m, n, :TE)
        push_heap!(heap, ModeSort(kc, :TE, m, n))
        kc = kc_cwg(dims[1], m, n, :TM)
        push_heap!(heap, ModeSort(kc, :TM, m, n))
    else
        kc = kc_cwg(dims[1], m, n, :TE)
        push_heap!(heap, ModeSort(kc, :TEc, m, n))
        push_heap!(heap, ModeSort(kc, :TEs, m, n))
        kc = kc_cwg(dims[1], m, n, :TM)
        push_heap!(heap, ModeSort(kc, :TMc, m, n))
        push_heap!(heap, ModeSort(kc, :TMs, m, n))
    end

    return nothing
end

function addmode_coaxial!(heap, visited, dims, m, n, T; m_max = nothing)
    
    n <= 0 && return nothing
    key = (m, n)
    push!(visited, key)

    if T == :TEM 
        push_heap!(heap, ModeSort(0.0, T, 0, 0))
        return nothing
    end
    n <= 0 && return nothing
    
    if m == 0
        kc = kc_coax(dims[1], dims[2], m, n, :TE)
        push_heap!(heap, ModeSort(kc, :TE, m, n))
        kc = kc_coax(dims[1], dims[2], m, n, :TM)
        push_heap!(heap, ModeSort(kc, :TM, m, n))
    else
        kc = kc_coax(dims[1], dims[2], m, n, :TE)
        push_heap!(heap, ModeSort(kc, :TEc, m, n))
        push_heap!(heap, ModeSort(kc, :TEs, m, n))
        kc = kc_coax(dims[1], dims[2], m, n, :TM)
        push_heap!(heap, ModeSort(kc, :TMc, m, n))
        push_heap!(heap, ModeSort(kc, :TMs, m, n))
    end

    return nothing
end

function addmode_radial!(heap, visited, dims, m, n, T; m_max = nothing)
    
    if !isnothing(m_max)
        m > m_max && return
    end

    key = (m, n)
    push!(visited, key)

    if T == :TEM || (m == 0 && n == 0)
        push_heap!(heap, ModeSort(0.0, :TEM, 0, 0))
        return nothing
    end

    if n >= 1
        kc = phase_constant_radial(dims[1], n)
        push_heap!(heap, ModeSort(kc, :TE, m, n))
        push_heap!(heap, ModeSort(kc, :TM, m, n))
    end

    return nothing
end

addmode_wedge!(heap, visited ,dims, m, n, T; m_max = nothing) = addmode_radial!(heap, visited ,dims, m, n, T; m_max = m_max)

function addmode_elliptic!(heap, visited, dims, m, n, T; m_max = nothing)
    T == :TEM && return nothing
    n <= 0 && return nothing
    m >= 100 && return nothing

    key = (m, n)
    push!(visited, key)

    kc = kc_ewg(dims[1], dims[2], m, n, true, :TE)
    push_heap!(heap, ModeSort(kc, :TEe, m, n))
    kc = kc_ewg(dims[1], dims[2], m, n, true, :TM)
    push_heap!(heap, ModeSort(kc, :TMe, m, n))

    if m > 0
        kc = kc_ewg(dims[1], dims[2], m, n, false, :TE)
        push_heap!(heap, ModeSort(kc, :TEo, m, n))
        kc = kc_ewg(dims[1], dims[2], m, n, false, :TM)
        push_heap!(heap, ModeSort(kc, :TMo, m, n))
    end
    
    return nothing
end



"""
    first_N_modes(addmode!::F, N, dims; m_max = nothing, initial_m = 3) where F

Generic function to find the first `N` modes for any waveguide geometry using a heap-based algorithm.

# Arguments
- `addmode!`: Function that adds mode(s) to the heap for a given (m,n) index pair
- `N`: Number of desired modes
- `dims`: Tuple of waveguide dimensions (varies by geometry type)
- `m_max`: Optional maximum azimuthal index (useful for radial/wedge waveguides)
- `initial_m`: Initial search range for m index (default: 3)

# Waveguide-specific functions
Different waveguide types use different `addmode!` functions and dimensions:
- Rectangular: `addmode_rectangular!`, `dims = (a, b)` → `kc_rwg(a, b, m, n)`
- Circular: `addmode_circular!`, `dims = (r,)` → `kc_cwg(r, m, n, T)`
- Coaxial: `addmode_coaxial!`, `dims = (r1, r2)` → `kc_coax(r1, r2, m, n, T)`
- Radial: `addmode_radial!`, `dims = (h,)` → `phase_constant_radial(h, n)`
- Wedge: `addmode_wedge!`, `dims = (h,)` → `phase_constant_radial(h, n)`
- Elliptic: `addmode_elliptic!`, `dims = (a, b)` → `kc_ewg(a, b, m, n, even, T)`

# Implementation details
The algorithm uses a min-heap to efficiently find modes in order of increasing cutoff frequency,
avoiding the need to pre-compute all possible modes. The result is sorted and the first N modes
are returned.
"""
function first_N_modes(addmode!::F, N, dims; m_max = nothing, initial_m = 3) where F

    visited = Set{Tuple{Int64, Int64}}()

    temp_length = round(Int, 1.05*N)
    heap = Vector{ModeSort{Float64, Symbol, Int64}}()
    sizehint!(heap, 2*N)
    result = Vector{ModeSort{Float64, Symbol, Int64}}()
    sizehint!(result, temp_length)

    # Initialize 
    addmode!(heap, visited, dims, 0, 0, :TEM; m_max = m_max)
    for m in 0:initial_m
        for n in 0:1
            addmode!(heap, visited, dims, m, n, :TE; m_max = m_max)
        end 
    end

    while length(result) < temp_length && !isempty(heap)
        e = pop_heap!(heap)

        m, n, kind = e.m, e.n, e.kind
        push!(result, e)
        kind == :TEM && continue

        next_m = (m+1, n)
        if !(next_m in visited)
            addmode!(heap, visited, dims, m+1, n, e.kind; m_max = m_max)
        end

        next_n = (m, n+1)
        if !(next_n in visited)
            if m == 0 || m == 1
                addmode!(heap, visited, dims, m, n+1, e.kind; m_max = m_max)
            end
        end

    end
    if length(result) < N
        @error "Not enough modes found."
    end
    return sort!(result)[1:N]
end


"""
    first_n_modes_rwg(N, a, b)

Returns the first `N` modes for a rectangular waveguide with width `a` and height `b`.

# Returns
A vector of tuples `(kind, m, n, kc)` where:
- `kind`: Mode type (`:TE` or `:TM`)
- `m`, `n`: Mode indices (number of half-wavelength variations in x and y directions)
- `kc`: Cutoff wavenumber

# Mode restrictions
- TE modes: `m + n > 0` (at least one index must be non-zero)
- TM modes: `m ≥ 1` and `n ≥ 1` (both indices must be positive)

The dominant mode is TE₁₀ (or TE₀₁ if b > a).

# Performance
This implementation uses a heap-based algorithm and is much faster than brute-force sorting
over all possible (m,n) combinations, especially for large N.
"""
function first_n_modes_rwg(N, a, b)
    if N <= 0
        return Tuple{Symbol, Int, Int, Float64}[]
    end
    k = a/b

    heap = Tuple{Float64, Int, Int}[]
    visited = Set{Tuple{Int, Int}}()

    push_heap!(heap, (1.0, 1, 0))
    push_heap!(heap, (k * k, 0, 1))
    push!(visited, (1, 0))
    push!(visited, (0, 1))
    result = Tuple{Symbol, Int, Int, Float64}[]

    while length(result) < N && !isempty(heap)
        e = pop_heap!(heap)
        val, m, n = e
        fc = kc_rwg(a, b, m, n)
        push!(result, (:TE, m, n, fc))
        if m >= 1 && n >= 1 && length(result) < N
            push!(result, (:TM, m, n, fc))
        end

        next_m = (m + 1, n)
        if !(next_m in visited)
            push!(visited, next_m)
            new_val = (m + 1)^2 + (k * n)^2
            push_heap!(heap, (new_val, m + 1, n))
        end

        if m == 0
            next_n = (0, n + 1)
            if !(next_n in visited)
                push!(visited, next_n)
                new_val = (k * (n + 1))^2
                push_heap!(heap, (new_val, 0, n + 1))
            end
        end
    end

    return result[1:N]
end


"""
    first_n_modes_cwg(N, r)

Returns the first `N` modes for a circular waveguide of radius `r`.

# Returns
A vector of `(kind, m, n, kc)` where:
- `kind`: Mode type (`:TE`, `:TM`, `:TEc`, `:TEs`, `:TMc`, `:TMs`)
- `m`, `n`: Mode indices
- `kc`: Cutoff wavenumber

# Mode notation
For `m > 0`, there are two degenerate polarizations with the same cutoff frequency:
- `:TEc`, `:TMc`: Cosine modes (cos(mθ) variation)
- `:TEs`, `:TMs`: Sine modes (sin(mθ) variation)

For `m = 0`, only `:TE` and `:TM` are returned (no polarization degeneracy).
"""
function first_n_modes_cwg(N, r)
    return map(x->(x.kind, x.m, x.n, x.kc), first_N_modes(addmode_circular!, N, (r,)))
end


"""
    first_n_modes_coax(N, r1, r2)

Returns the first `N` modes for a coaxial waveguide with inner radius `r2` and outer radius `r1`.

# Returns
A vector of`(kind, m, n, kc)` where:
- `kind`: Mode type (`:TEM`, `:TE`, `:TM`, `:TEc`, `:TEs`, `:TMc`, `:TMs`)
- `m`, `n`: Mode indices
- `kc`: Cutoff wavenumber

# Mode notation
The TEM mode (`:TEM`) has zero cutoff frequency and is the dominant mode.

For `m > 0`, there are two degenerate polarizations with the same cutoff frequency:
- `:TEc`, `:TMc`: Cosine modes (cos(mθ) variation)
- `:TEs`, `:TMs`: Sine modes (sin(mθ) variation)

For `m = 0`, only `:TE` and `:TM` are returned (no polarization degeneracy).
"""
function first_n_modes_coax(N, r1, r2)
    return map(x->(x.kind, x.m, x.n, x.kc), first_N_modes(addmode_coaxial!, N, (r1,r2)))
end


"""
    first_n_modes_radial(N, h)

Returns the first `N` modes for a radial waveguide with height `h` between parallel plates.

# Returns
A vector of `(kind, m, n, kc)` where:
- `kind`: Mode type (`:TEM`, `:TE`, `:TM`)
- `m`: Azimuthal mode index
- `n`: Vertical mode index
- `kc`: Cutoff wavenumber (radial propagation constant)

# Note
The TEM mode (m=0, n=0) has zero cutoff frequency and corresponds to radial propagation
between parallel plates.
"""
function first_n_modes_radial(N, h)
    return map(x->(x.kind, x.m, x.n, x.kc),  first_N_modes(addmode_radial!, N, (h,)))
end


"""
    first_n_modes_ewg(N, a, b)

Returns the first `N` modes for an elliptic waveguide with semi-major axis `a` and semi-minor axis `b`.

# Returns
A vector of `(kind, m, n, kc)` where:
- `kind`: Mode type (`:TEe`, `:TEo`, `:TMe`, `:TMo`)
- `m`, `n`: Mode indices
- `kc`: Cutoff wavenumber

# Mode notation
Unlike circular or coaxial waveguides, elliptic waveguides have even and odd modes that are
**not degenerate** (they have different cutoff frequencies):
- `:TEe`, `:TMe`: Even modes (using even Mathieu functions, cosine-elliptic)
- `:TEo`, `:TMo`: Odd modes (using odd Mathieu functions, sine-elliptic)

For `m = 0`, only even modes exist (`:TEe`, `:TMe`).
For `m > 0`, both even and odd modes exist with **different** cutoff frequencies.
"""
function first_n_modes_ewg(N, a, b)
    return map(x->(x.kind, x.m, x.n, x.kc), first_N_modes(addmode_elliptic!, N, (a, b)))
end

