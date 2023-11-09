using FFTW

function make_metric(psi_grid, theta_grid, m_grid, ro, zo, rzphi)
    """

    """
    metric_names = ["g11", "g22", "g33", "g23", "g31", "g12", "jac", "jac_prime"]
    fourfit_metric_names = ["G11","G22","G33","G23","G31","G12","Jmat","Jmat_prime"]
    metric = Dict{String, Matrix{Float64}}()
    for name in metric_names
        metric[name] = zeros(length(psi_grid), length(theta_grid))
    end

    r = sqrt(rzphi["r_squared"](psi_grid, theta_grid, grid=true))
    eta = 2*pi*(theta_grid + rzphi["delta_eta"](psi_grid, theta_grid, grid=true))
    R = ro + r*cos(eta)
    jac = rzphi["jac"](psi_grid, theta_grid, grid=true)
    jac1 = rzphi["jac"](psi_grid, theta_grid, grid=true)

    v = zeros(3,3, length(psi_grid), length(theta_grid))

    """Copmute contravariant basis vectors"""
    v[1,1,:,:] .= rzphi["r_squared"](psi_grid, theta_grid, grid=true, dx=1) / ( 2 *r .* jac)        # 0, 5
    v[1,2,:,:] .= rzphi["delta_eta"](psi_grid, theta_grid, grid=true, dx=1) * 2 * pi .*r ./ jac     # 0, 5
    v[1,3,:,:] .= rzphi["delta_phi"](psi_grid, theta_grid, grid=true, dx=1) .* R ./ jac             # 0, 5, 4
    v[2,1,:,:] .= rzphi["r_squared"](psi_grid, theta_grid, grid=true, dy=1) / (2 .* r .* jac)       # 1. 5
    v[2,2,:,:] .= (1+rzphi["delta_eta"](psi_grid, theta_grid, grid=true, dy=1)) * 2*pi .*r ./ jac   # 1, 5
    v[2,3,:,:] .= rzphi["delta_phi"](psi_grid, theta_grid, grid=true, dy=1) .*R ./jac               # 1, 3, 5
    v[3,3,:,:] .= 2*pi .* R ./ jac                                                                  # 2, 3, 4

    """Compute metric tensor components"""
    metric["g11"] .= sum(v[1,:,:,:] .^ 2, dims=1) .* jac            # g11
    metric["g22"] .= sum(v[2,:,:,:] .^ 2, dims=1) .* jac            # g22
    metric["g33"] .= v[3,3,:,:] .^ 2 .* jac                         # g33
    metric["g23"] .= v[2,3,:,:] .* v[3,3,:,:] .* jac                # g23
    metric["g31"] .= v[3,3,:,:] .* v[1,3,:,:] .* jac                # g31
    metric["g12"] .= sum(v[1,:,:,:] .* v[2,:,:,:], dims=1) .* jac   # g12
    metric["jac"] .= jac
    metric["jac_prime"] .= jac1

    delta_m = (m_grid' .- m_grid)' # TODO: check this, I think m_grid might be a boolean and Julia might not like this
    fourfit_metric = Dict{String, Matrix{Complex{Float64}}}()
    for (name1, name2) in zip(fourfit_metric_names, metric_names) # TODO: This loop is a little funky, check it
        metric_subset = metric[name2][:, 1:end] / (length(theta_grid) -1)
        fourfit_metric[name1] = fft(metric_subset, 1)[:, delta_m]
    end

    return metric, fourfit_metric
end