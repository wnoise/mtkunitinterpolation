using ModelingToolkit
using DataInterpolations
using Unitful
using Random
using DifferentialEquations

@variables t, [unit=u"s"]
@variables T(t), [unit=u"K"]

function generate_data()
    rng = MersenneTwister(54321)

    time = 0.0u"s":10.0u"minute":24.0u"hr"
    N = length(time)
    A = 280.0u"K" .+ rand(rng, N)*u"K"
    B = rand(rng, N)*u"V"

    data = (; time, A, B)
end


function interpolation_wrapper(vals, val_unit, args, arg_unit)
    interpolator = LinearInterpolation(vals, args)
    function wrapper(x)
        interpolator(x)
    end
    # this is tricky but appears to work; can only be handle at top-level
    @eval @register_symbolic wrapper(x)

    # Should probably check args is 1-tuple (arg_unit,) (i.e. u"s")
    # ERROR: Can only be evaluated at top-level.  @eval loses val_unit binding
    function ModelingToolkit.get_unit(wrapper, args)
        return val_unit
    end
    return wrapper
end

@connector function Pin(;name)
    sts = @variables begin
        V(t), [unit=u"V"]
        I(t), [unit=u"A", connect=Flow]
    end
    ODESystem(Equation[], t, sts, []; name=name)
end

function Source(data; name)
    Aitp = interpolation_wrapper(ustrip.(data.A), u"K", ustrip.(data.time), u"s")
    Bitp = interpolation_wrapper(ustrip.(data.B), u"V", ustrip.(data.time), u"s")
    @variables T(t), [unit=u"K"]
    T = GlobalScope(T)
    @named out = Pin()
    @named gnd = Pin()
    eqs = [
        out.V - gnd.V ~ Bitp(t)
        0 ~ out.I + gnd.I
        T ~ Aitp(t)
    ]
    compose(ODESystem(eqs, t, [T], []; name=name), out, gnd)
end

function Ground(;name)
    @named a=Pin()
    eqs = [
        a.V ~ 0
    ]
    compose(ODESystem(eqs, t, [], []; name=name), a)
end

function Resistor(R₀;name)
    @named a=Pin()
    @named b=Pin()
    @parameters begin
        R₀=R₀, [unit=u"Ω"]
        γ=0.1, [unit=u"Ω/K"]
        T₀=300.0, [unit=u"K"]
    end
    @variables R(t), [unit=u"Ω"]
    @variables T(t), [unit=u"K"]
    T = GlobalScope(T)
    eqs = [
        R ~ R₀ + (T-T₀)*γ
        a.V - b.V ~ R*a.I
        0 ~ a.I + b.I
    ]
    compose(ODESystem(eqs, t, [T], [T₀, R₀, γ]; name=name), a, b)
end

function Capacitor(C;name)
    @named a=Pin()
    @named b=Pin()
    @parameters begin
        C₀=C, [unit=u"F"]
        γ=0.1, [unit=u"F/K"]
        T₀=300.0, [unit=u"K"]
    end
    @variables V(t), [unit=u"V"]
    @variables T(t), [unit=u"K"]
    T = GlobalScope(T)
    D = Differential(t)
    C = C₀ + (T - T₀)*γ
    eqs = [
       V ~ a.V - b.V
       a.I + b.I ~ 0
       a.I ~ C*D(V)
    ]
    compose(ODESystem(eqs, t, [V, T], [T₀, C₀, γ]; name=name), a, b)
end

function example()
    data = generate_data()
    @named src = Source(data)
    @named R = Resistor(1.0)
    @named C = Capacitor(1.0)
    @named gnd = Ground()
    @variables T(t), [unit=u"K"]
    T = GlobalScope(T)
    eqs = [
        connect(src.out, R.a)
        connect(R.b, C.a)
        connect(C.b, gnd.a, src.gnd)
    ]
    sys = compose(ODESystem(eqs, t, [], []; name=:example), src, R, C, gnd)
    sys = structural_simplify(sys)
    u0 = [ C.V => 0.0, T => ustrip(data.A[begin]) ]
    prob = ODAEProblem(sys, u0, extrema(ustrip.(data.time)))
    soln = solve(prob)
end
