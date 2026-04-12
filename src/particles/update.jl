xss_update(::LCS.RKStage1, props::Properties) = props.xss2
xss0_update(::LCS.RKStage1, props::Properties) = props.xss

xss_update(::LCS.RKStage2, props::Properties) = props.xss
xss0_update(::LCS.RKStage2, props::Properties) = props.xss

uss_update(::LCS.RKStage1, props::Properties) = props.uss2
uss0_update(::LCS.RKStage1, props::Properties) = props.uss

uss_update(::LCS.RKStage2, props::Properties) = props.uss
uss0_update(::LCS.RKStage2, props::Properties) = props.uss

"""
    update!(buf, state, config, topo, stage, iprofile)

Update particle positions and velocities using explicit Runge-Kutta time integration.
"""
function update!(
    buf::LCS.Buffer, state::LCS.State, config::LCS.Config, ::Topologies.Topology, stage::LCS.RKStage, iprofile::Integer
)
    update!(;
        xss=xss_update(stage, buf.particles[iprofile].props),
        xss0=xss0_update(stage, buf.particles[iprofile].props),
        uss=uss_update(stage, buf.particles[iprofile].props),
        uss0=uss0_update(stage, buf.particles[iprofile].props),
        buf.particles[iprofile].props,
        dt=LCS.dtrk(stage, state),
        state.particles[iprofile].nvalid,
        config.backend,
    )
end
function update!(;
    xss::Tuple3{Property{<:Real}},
    xss0::Tuple3{Property{<:Real}},
    uss::Tuple3{Property{<:Real}},
    uss0::Tuple3{Property{<:Real}},
    props::Properties,
    dt::Real,
    nvalid::Integer,
    backend::KA.Backend,
)
    xs, ys, zs = xss
    xs0, ys0, zs0 = xss0
    us, vs, ws = uss
    us0, vs0, ws0 = uss0
    dudts, dvdts, dwdts = props.dudtss

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            x0, y0, z0 = xs0[i], ys0[i], zs0[i]
            u0, v0, w0 = us0[i], vs0[i], ws0[i]
            dudt, dvdt, dwdt = dudts[i], dvdts[i], dwdts[i]

            xs[i] = x0 + dt * u0
            ys[i] = y0 + dt * v0
            zs[i] = z0 + dt * w0

            us[i] = u0 + dt * dudt
            vs[i] = v0 + dt * dvdt
            ws[i] = w0 + dt * dwdt
        end
    end

    LCS.EMPTY_LOG
end
