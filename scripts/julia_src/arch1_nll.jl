function arch1_nll(x::Array{Float64}, theta::Array{Float64, 2})
    n = length(x)
    xbar = 1 / n * sum(x)

    diff2 = (x .- xbar).^2

    uv = 1 / n * sum(diff2)

     function nll(thet)
        a0 = thet[1]
        a1 = thet[2]

        cv = a0 + a1 * uv

        ll = -1 / 2 * (log(cv) + diff2[1] / cv)

        for i in 2:n
            ll += -1 / 2 * (log(cv) + diff2[i] / cv)
        end

        return -ll
    end

    res = optimize(nll, theta)
    # return Optim.minimizer(res)
    return res
end
