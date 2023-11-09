

function ode_itime(coeffs, pt1, sing_loc; pt2=nothing)
    """
    Computes estimated cost to integrate from pt1 to pt2 using adaptive method 
    (ie, number of steps to take over the interval)
    """
    a_axis = coeffs[1]
    b_axis = coeffs[2]
    a_sing = coeffs[3]
    b_sing = coeffs[4]
    a_edge = coeffs[5]
    b_edge = coeffs[6]

    itime = similar(pt1, zero(eltype(pt1)))

    if pt2 !== nothing
        # estimate integration time of the whole interval
        itime += (a_axis/b_axis) * abs(log(1+b_axis*abs(pt2-sing_loc[1])) - log(1+b_axis*abs(pt1-sing_loc[1])))
        itime += (a_edge/b_edge) * abs(log(1+b_edge*abs(pt2-sing_loc[end])) - log(1+b_edge*abs(pt1-sing_loc[end])))
        for s in sin_loc[2:end-1]
            itime += (a_sing/b_sing) * abs(log(1+b_sing*abs(pt2-s)) - log(1+b_sing*abs(pt1-s)))
        end
    else
        # esimate instananeous time of current location
        itime += a_axis/(1+b_axis*abs(pt1-sing_loc[1]))
        itime += a_edge/(1+b_axis*abs(pt1-sing_loc[end]))
        for s in sing_loc[2:end-1]
            itime += a_sing/(1+b_sing*abs(pt1-s))
        end
    end
    return itime
end

function set_intervals(sing, itime_coeffs, start_ind, end_ind,nInters; method="naive")
    """
    Divides domains into subintervals for integration in parallel
    """ 

    # Make sure there are enough intervals so that each will touch at most 1 singularity
    @assert nInters >= length(sing) * 3 + 2
    # each row defines one interval: [start, stop, sing_left, sing_right, direction, q]
    inters = fill(NaN, (length(sing) * 2 + 2, 6))
    sing_loc = [0; mean(sing[:,1:2], dims=2); 1]

    # TODO: finish function

end