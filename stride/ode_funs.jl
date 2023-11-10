

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
        # Estimate integration time of the whole interval
        itime += (a_axis/b_axis) * abs(log(1+b_axis*abs(pt2-sing_loc[1])) - log(1+b_axis*abs(pt1-sing_loc[1])))
        itime += (a_edge/b_edge) * abs(log(1+b_edge*abs(pt2-sing_loc[end])) - log(1+b_edge*abs(pt1-sing_loc[end])))
        for s in sin_loc[2:end-1]
            itime += (a_sing/b_sing) * abs(log(1+b_sing*abs(pt2-s)) - log(1+b_sing*abs(pt1-s)))
        end
    else
        # Estimate instananeous time of current location
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
    # Each row defines one interval: [start, stop, sing_left, sing_right, direction, q]
    inters = fill(NaN, (length(sing) * 2 + 2, 6))
    sing_loc = [0; mean(sing[:,1:2], dims=2); 1]

    # First set minimal intervals defined by axis and singular surfaces. 
    # TODO: I don't trust this
    temp = vcat(start_ind, reshape(sing[:, 1:2], (length(sing), 2))', end_ind)'
    for (i, elem) in enumerate(eachrow(temp))
        inters[2*i-1,1] = elem[1]
        inters[2*i-1,2] = elem[1] + (elem[2] - elem[1])/2
        inters[2*i-1,3:end] = [1,0,1,0]
        inters[2*i,1] = elem[1] + (elem[2] - elem[1])/2
        inters[2*i,2] = elem[2]
        inters[2*i,3:end] = [0,1,-1,0]
    end

    # Evenly subdivide intervals further
    while length(inters) < nInters-length(sing)
        if method == "naive"
            # Find the largest interval
            interval_lengths = inters[:,2] - inters[:,1]
            max_ind = argmax(interval_lengths)
            split_pt = inters[max_ind,1] + interval_lengths[max_ind]/2
        elseif method == "sing"
            # Find the inverval that will take the most steps
            interval_times = ode_itime(itime_coeffs, inters[:,1] sing_loc, inters[:,2])
            max_ind = argmax(interval_times)
            # Find where to split to make each half take the same time
            s1 = odeitime(itime_coeffs, inters[max_ind,1], sing_loc)
            s2 = odeitime(itime_coeffs, inters[max_ind,2], sing_loc)
            alpha = (2*s1 - sqrt(2* s1^2 + 2* s2^2)) / (2*(s1-s2))
            split_pt = alpha*inters[max_ind,2] + (1-alpha)*inters[max_ind,1]
        end

        # Divide that interval, keeping track of where singularities are
        if inters[max_ind][3] == 1 # singularity to the left of interval, split to the sing_right
            temp = [split_pt, inters[max_ind,2], 0, 0, inters[max_ind, end-1], 0]
            inters[max_ind,1] = split_pt
            inters = vcat(inters[:max_ind+1], temp', inters[max_ind+1:end]) # TODO: make sure vcats are done correctly
        elseif inters[max_ind][4] == 1 # singularity to the right of interval, split to the eft
            temp = [inters[max_ind, 1], split_pt, 0, 0, inters[max_ind, end-1], 0]
            inters[max_ind,0] = split_pt
            inters = vcat(inters[:max_ind], temp', inters[max_ind:end])
        else # no singularity on either side, doesn't matter which way you split
            temp = [split_pt, inters[max_ind,2], 0, 0, inters[max_ind, end-1], 0]
            inters[max_ind,1] = split_pt
            inters = vcat(inters[1:max_ind, :], temp', inters[max_ind+1:end,:])
        end
    end

    # Insert singular intervals
    for s in sing
        idx = searchsortedfirst(inters[:,1], s[1])
        row_to_insert = [s[1], s[2], 1, 1, 0, s[3]]
        inters = vcat(inters[1:idx, :], row_to_insert', inters[idx+1:end, :])
    end

    return inters
end

function solve(inters, L, mpert, method)
    # TODO
end

function fixup(uT)
    # TODO
end

function prpagate(all_soln, inters, psio)
    # TODO
end

function wrapper(nInters, ode_method, interval_method, start_ind, end_ind, sing, mpert, L, psio, itime_coeffs)
    """
    Wrapper function that computes intervals, integrates ODE, combines solutions
    """
    inters = set_intervals(sing, itime_coeffs, start_ind, end_ind, nInters, method=interval_method)
    all_soln, status, nfev, njev, nlu, nmm, nsteps = solve(inters, L, mpert, ode_method)
    Wp, soln = propagate(all_soln, inters, psio)
    stats = {"nfev": nfev,
             "njev": njev,
             "nlu": nlu,
             "nmm": nmm,
             "nsteps": nsteps}
             "flops": 38/3*nfev + 8/3*nlu + 2*(nmm+1) # 1 extra matmul per interval for propagation
end