function arch(x::Float64, q::Int64, θ0::Float64)
    #=
        function to estimate an ARCH(q) model on x, where X_t = σ_t * ϵ_t,
        with ϵ_t ~ N(0,1)
        Args:
            x: data vector
            q: lag order of the ARCH model
            θ0: initial parameters for optimisation procedure
        Returns:
            Optim output
    =#

    function nll(θ::Float64)::Float64
        #=
            function to calculate the negative log-likelihood value of an ARCH(q)
            model at θ
            N, x, σ2, q of the ARCH(q) are required in the parent environment
            σ2 in the parent environment is mutated by this function
            Args:
                θ: vector of parameters [ω, α_1, ..., α_q] at which to evaluate
            Returns:
                the numerical value of the negative logL
        =#

        ω = θ[1]
        α = θ[collect(2:(q + 1))]

        # initialise logL sum
        ll = 0

        for t in (q + 1):N
            # collect indices for correct temporal subsetting
            t_i = collect((t - 1):-1:(t - q))
            # estimate conditional variance at t
            σ2[t] = ω + α' * x2[t_i]
            # add value of logL at time t to the sum
            ll += -1 / 2 * (log(σ2[t]) + x2[t] / σ2[t])
        end

        # return negative of log-likelihood for minimisation
        return -ll
    end

    function ∇nll!(G::Float64, θ::Float64)
        #=
            function to evaluate the gradient of the negative log-likelihood of
            an ARCH(q) model
            N, x, σ2, q of the ARCH(q) are required in the parent environment
            Args:
                G: fixed-size storage array that is mutated by the function; to
                    avoid creating a new array with every function call
                θ: vector of parameters [ω, α_1, ..., α_q] at which to evaluate
            Returns:
                a mutation of G with the numerical values of the gradient
                evaluated at θ
        =#

        ω = θ[1]
        α = θ[collect(2:(q + 1))]

        # initialise σ^2 with unconditional variance
        σ2 = fill(1 / N * sum(x2), N)

        # calculate conditional variance
        for t in (q + 1):N
            t_i = collect((t - 1):-1:(t - q))
            σ2[t] = ω + α' * x2[t_i]
        end

        # derivative of -logL w.r.t. ω
        G[1] = 1 / 2 * sum((1 ./ σ2) .- x2 ./ σ2.^2)

        # derivates w.r.t. α_i
        function ∇_helper(i)
            #=
                function to evaluate the derivative of -logL w.r.t. α_i
                Args:
                    i: the α index for which to evaluate the derivate
                Returns:
                    the numerical value of the derivate
            =#

            z = 1 / 2 .* ((1 ./ σ2) .- x2 ./ σ2.^2)
            return x2[collect((N - i):-1:1)]' * z[collect(N:-1:(i + 1))]
        end

        # mutate G with numerical values of the partial derivates
        for i in 2:(q + 1)
            G[i] = ∇_helper(i - 1)
        end
    end

    # estimate model
    N = length(x)

    # square the zero mean series x
    x2 = x.^2

    # initialise with unconditional variance to not have 0 for the first q elements
    σ2 = fill(1 / N * sum(x2), N)

    # define parameter bounds, 0 < θ < Inf
    lower = zeros(q + 1)
    upper = fill(Inf, q + 1)
    inner_optimizer = LBFGS()
    res = optimize(nll, ∇nll!, lower, upper, θ0, Fminbox(inner_optimizer))

    #return (res, σ2)
    return res

end
