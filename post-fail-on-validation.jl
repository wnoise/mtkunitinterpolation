using ModelingToolkit
using DataInterpolations
using Unitful
using Random
using DifferentialEquations

@variables t, [unit=u"s"]

function generate_data()
    rng = MersenneTwister(54321)

    time = 0.0u"s":10.0u"minute":24.0u"hr"
    N = length(time)
    A = 280.0u"K" .+ rand(rng, N)*u"K"
    B = rand(rng, N)*u"V"

    data = (; time, A, B)
end

@connector function Pin(;name)
    sts = @variables begin
        V(t), [unit=u"V"]
        I(t), [unit=u"A", connect=Flow]
    end
    ODESystem(Equation[], t, sts, []; name=name)
end

function Source(data; name)
    Aitp = LinearInterpolation(ustrip.(data.A), ustrip.(data.time))
    Bitp = LinearInterpolation(ustrip.(data.B), ustrip.(data.time))
    @parameters T=300.0, [unit=u"K"]
    T = GlobalScope(T)
    @named out = Pin()
    @named gnd = Pin()
    evts = [data.time => [T ~ Aitp(t)]]
    eqs = [
        out.V - gnd.V ~ Bitp(t)
        0 ~ out.I + gnd.I
    ]
    compose(ODESystem(eqs, t, [], [T]; name=name), out, gnd)
end

function Ground(;name)
    @named a=Pin()
    eqs = [
        a.V ~ 0
    ]
    compose(ODESystem(eqs, t, [], []; name=name), a)
end

function Resistor(R;name)
    @named a=Pin()
    @named b=Pin()
    @parameters begin
        R₀=R, [unit=u"Ω"]
        γ=0.1, [unit=u"Ω/K"]
        T₀=300.0, [unit=u"K"]
        T, [unit=u"K"]
    end
    T = GlobalScope(T)
    R = R₀ + (T-T₀)*γ
    eqs = [
        a.V - b.V ~ R*a.I
        0 ~ a.I + b.I
    ]
    compose(ODESystem(eqs, t, [], [T₀, T, R₀, γ]; name=name), a, b)
end

function Capacitor(C;name)
    @named a=Pin()
    @named b=Pin()
    @parameters begin
        C₀=C, [unit=u"F"]
        γ=0.1, [unit=u"F/K"]
        T₀=300.0, [unit=u"K"]
        T, [unit=u"K"]
    end
    @variables V(t), [unit=u"V"]
    D = Differential(t)
    T = GlobalScope(T)
    C = C₀ + (T - T₀)*γ
    eqs = [
       V ~ a.V - b.V
       a.I + b.I ~ 0
       a.I ~ C*D(V)
    ]
    compose(ODESystem(eqs, t, [V], [T₀, T, C₀, γ]; name=name), a, b)
end

function example()
    data = generate_data()
    @named src = Source(data)
    @named R = Resistor(1.0)
    @named C = Capacitor(1.0)
    @named gnd = Ground()
    eqs = [
        connect(src.out, R.a)
        connect(R.b, C.a)
        connect(C.b, gnd.a, src.gnd)
    ]
    sys = compose(ODESystem(eqs, t, [], []; name=:example), src, R, C, gnd)
    sys = structural_simplify(sys)
    u0 = [ C.V => 0.0 ]
    prob = ODAEProblem(sys, u0, extrema(data.time))
    soln = solve(prob)
end
