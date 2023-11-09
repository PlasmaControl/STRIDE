


""" Replacing Bfield Class and functions """
struct Bfield
    psirz
    psio
    F_func
    P_func
end
# TODO: I'm unsure if all the dr-9 and dz=0 are correct, should maybe be 1
# TODO: same with dpsi=0
function psi(b::Bfield, r, z; dr=0, dz=0, grid=false)
    return b.psirz(r, z, dx=dr, dy=dz, grid=grid)
end
function F(b::Bfield, r, z; dpsi=0, grid=false)
    return b.F_func(1-psi(b, r, z, grid=grid) ./ b.psio, nu=dpsi)
end
function F_prime(b::Bfield, r, z; dpsi=0, grid=false)
    return b.F_func(1-psi(b, r, z, grid=grid) ./ b.psio, nu=dpsi+1)
end
function P(b::BField, r, z; dpsi=0, grid=false)
    return b.P_func(1-psi(b, r, z, grid=grid) ./ b.psio, nu=dpsi)
end
function P_prime(b::BField, r, z; dpsi=0, grid=false)
    return b.P_func(1-psi(b, r, z, grid=grid) ./ b.psio, nu=dpsi+1)
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

