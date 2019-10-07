arch1_nll = function(x, theta)
    n = length(x)
    xbar = 1 / n * sum(x)

    diff2 = (x .- xbar).^2

    uv = 1 / n * sum(diff2)

    nll = function(thet)
        a0 = thet[1]
        a1 = thet[2]

        cv = a0 + a1 * uv

        ll = -1 / 2 * (log(cv) + diff2[1] / cv)

        for i in 2:n
            ll += -1 / 2 * (log(cv) + diff2[i] / cv)
        end

        return -ll
    end

    return optimize(nll, theta)
end
