

function plot_gfile(g)
    """
    
    """

    if typeof(g) == String
        g = read_gfile(g)
    end

    

end

function write_gfile(filename, date, shot, time, efit, nw, nh, rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, axis_flux, boundary_flux, bcentr, current, fpol, pres, ffprime, pprime, psirz, qpsi)
    """

    """


end

function read_gfile(filename)
    """

    """

    g = Dict()

    return g
end