Us_motion(::LCS.RKStage1, flowbuf::Flows.Fields) = flowbuf.Us
Us_motion(::LCS.RKStage2, flowbuf::Flows.Fields) = flowbuf.Us2

xss_motion(::LCS.RKStage1, props::Properties) = props.xss
xss_motion(::LCS.RKStage2, props::Properties) = props.xss2

uss_motion(::LCS.RKStage1, props::Properties) = props.uss
uss_motion(::LCS.RKStage2, props::Properties) = props.uss2

"""
    motion!(buf, state, config, topo, stage, iprofile)

Compute particle equations of motion.
"""
function motion!(
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology,
    stage::LCS.RKStage,
    iprofile::Integer,
)
    motion!(;
        Us=Us_motion(stage, buf.flow.fields),
        xss=xss_motion(stage, buf.particles[iprofile].props),
        uss=uss_motion(stage, buf.particles[iprofile].props),
        buf.particles[iprofile].props,
        state.particles[iprofile].nvalid,
        config.grid,
        fparams=config.flow.params,
        pparams=config.particles[iprofile].params,
        config.particles[iprofile].drag,
        config.backend,
        topo,
    )
end

function motion!(;
    Us::Tuple3{Field},
    xss::Tuple3{Property{<:Real}},
    uss::Tuple3{Property{<:Real}},
    props::Properties,
    nvalid::Integer,
    grid::LCS.Grid,
    fparams::Flows.FlowParams,
    pparams::Particles.ParticleParams,
    drag::DragModel,
    backend::KA.Backend,
    topo::Topologies.Topology,
)
    U, V, W = Us
    xs, ys, zs = xss
    us, vs, ws = uss
    dudts, dvdts, dwdts = props.dudtss
    diams = props.diams

    dx, dy, dz = LCS.spacings(grid)
    Re = Flows.Re(fparams)
    g = pparams.gravity
    ox, oy, oz = LCS.origins_l(topo)

    Parallel.foraxes(backend, (1:nvalid,)) do i
        @inbounds begin
            xp, yp, zp = xs[i], ys[i], zs[i]
            up, vp, wp = us[i], vs[i], ws[i]
            diam = diams[i]

            uf = interp(U, grid, xmapping(xp - ox, yp - oy, zp - oz, dx, dy, dz))
            vf = interp(V, grid, ymapping(xp - ox, yp - oy, zp - oz, dx, dy, dz))
            wf = interp(W, grid, zmapping(xp - ox, yp - oy, zp - oz, dx, dy, dz))

            vrel = Utils.norm((uf - up, vf - vp, wf - wp))
            Rep = diam * vrel * Re
            f = coeff(drag, Rep)

            taui = 1 / Particles.taup(pparams, fparams, diam)
            dudts[i] = f * (uf - up) * taui - g
            dvdts[i] = f * (vf - vp) * taui
            dwdts[i] = f * (wf - wp) * taui
        end
    end

    LCS.EMPTY_LOG
end
