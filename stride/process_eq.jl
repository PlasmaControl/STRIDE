


""" Replacing Bfield Class and functions """
struct Bfield
    psirz
    psio
    F_func
    P_func
end
# TODO: I'm unsure if all the dr=0 and dz=0 are correct, should maybe be 1
# TODO: same with dpsi=0
function psi(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr, dy=dz, grid=grid)
end
function F(b::Bfield, r, z; dpsi=0, grid=false)
    return b.F_func(1-psi(b, r, z, grid=grid) / b.psio, nu=dpsi)
end
function F_prime(b::Bfield, r, z; dpsi=0, grid=false)
    return b.F_func(1-psi(b, r, z, grid=grid) / b.psio, nu=dpsi+1)
end
function P(b::Bfield, r, z; dpsi=0, grid=false)
    return b.P_func(1-psi(b, r, z, grid=grid) / b.psio, nu=dpsi)
end
function P_prime(b::Bfield, r, z; dpsi=0, grid=false)
    return b.P_func(1-psi(b, r, z, grid=grid) / b.psio, nu=dpsi+1)
end
function psi_r(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr+1, dy=dz, grid=grid)
end
function psi_z(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr, dy=dz+1, grid=grid)
end
function Br(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return psi_z(b, r, z, dr=dr, dz=dz, grid=grid) ./ r
end
function Bz(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return -psi_r(b, r, z, dr=dr, dz=dz, grid=grid) ./ r
end
function psi_rr(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr+2, dy=dz, grid=grid)
end
function psi_rz(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr+1, dy=dz+1, grid=grid)
end
function psi_zz(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr, dy=dz+2, grid=grid)
end
function Br_r(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return (psi_rz(b, r, z, dr=dr, dz=dz, grid=grid) - Br(b, r, z, dr=dr, dz=dz, grid=grid)) ./ r
end
function Br_z(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return -psi_zz(b, r, z, dr=dr, dz=dz, grid=grid) ./ r
end
function Bz_r(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return -(psi_rz(b, r, z, dr=dr, dz=dz, grid=grid) + Bz(b, r, z, dr=dr, dz=dz, grid=grid)) ./ r
end
function Bz_z(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return -psirz(b, r, z, dr=dr, dz=dz, grid=grid) ./ r
end
function Bp(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return sqrt(Br(b, r, z, dr=dr, dz=dz, grid=grid) .^ 2 + Bz(b, r, z, dr=dr, dz=dz, grid=grid) .^ 2)
end
function Bt(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return F(b, r, z, dpsi=0, grid=grid) ./ r
end
function B(b::Bfield,r,z;grid=false)
    return sqrt(Bp(b,r,z,grid=grid) .^2 + Bt(b,r,z,grid=grid) .^2)
end

""" Other equilibrium functions """

function find_0_point(bf; guess=None)
    # TODO
end

function find_X_points(bf, ro, zo)
    # TODO
end

function convert_efit_equilibrium(g; mpsi=129, mtheta=129, psilow=0.01, psihigh=0.98, return_arrs=false)
    """
    Convert an EFIT equilibrium to grid and profiles for STRIDE
    """
    
    mu0 = pi*4e-7
    R_grid = LinRange(g["rlefft"], g["rleft"] + g["rdim"], g["nw"])
    Z_grid = LinRange(-g["zdim"]/2, g["zdim"] / 2, g["nh"])
    psio = g["boundary_flux"] - g["axis_flux"]
    psigrid_n = LinRange(0, 1, g["nw"])
    fpol = abs(g["fpol"])
    pres = max(g["pres"]*mu0, 0)
    qpsi = g["qpsi"]
    psirz_arr = g["boundary_flux"] - g["psirz"]

    if psio < 0
        psio = -psio
        psirz_arr = -psirz_arr
    end

    direct_out = process_direct_equilibrium(R_grid, Z_grid, psirz_arr, psigrid_in, pres, fpol, 
                                            psilow, psihigh, mpsi, mtheta, psio, return_arrs=return_arrs)

    
    psi_grid = direct_out[1]
    theta_grid = direct_out[2]
    straight_field_line_coords = direct_out[3]
    profiles = direct_out[4]
    ro = direct_out[5]
    zo = direct_out[6]
    if return_arrs
        straight_field_line_coords_arrs = direct_out[7]
        profiles_arrs = direct_out[8]
        bf = direct_out[9]
        temp_data = direct_out[10]
    end

    if return_arrs
        return psi_grid, theta_grid, straight_field_line_coords, profiles, ro, zo, psio, 
                         straight_field_line_coords_arrs, profiles_arrs, bf, temp_data
    else
        return psi_grid, theta_grid, straight_field_line_coords, profiles, ro, zo, psio
    end
end

function process_direct_equilibrium(R_grid, Z_grid, psirz_arr, psigrid_in, pres, fpol, psilow, psihigh, 
                                    mpsi, mtheta, psio; return_arrs=false)
    """
    TODO
    """

    # TODO: fill in rest

    if return_arrs
        return psi_grid, theta_grid, straight_field_line_coords, profiles, ro, zo, 
                         straight_field_line_coords_arrs, profiles_arrs, bf, temp_data
    else
        return psi_grid, theta_grid, straight_field_line_coords, profiles, ro, zo
    end
end