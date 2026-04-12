"""
    evolve!(buf, state, config, topo) -> NamedTuple

Advance simulation by one time step using 2-stage Runge-Kutta.
"""
function evolve!(buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    evolve!(config.mode, buf, state, config, topo)
end

function evolve!(state::State)
    state.step += 1
    state.t += state.dt
end

function evolve!(::FlowMode, buf::LCS.Buffer, state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    backend = config.backend
    stage = LCS.RKStage1()

    @profile backend "flow/synchalo/1" Topologies.synchalo!(
        Topologies.FullHalo(config.grid),
        (buf.flow.fields.Us..., buf.flow.fields.P),
        buf.halo,
        config.grid,
        backend,
        topo,
    )

    @profile backend "flow/rhs/1" f_rhs1 = Flows.rhs!(buf, state, config, topo, stage)
    @profile backend "flow/update/1" f_update1 = Flows.update!(buf, state, config, topo, stage)
    @profile backend "flow/couple/1" f_couple1 = Flows.couple!(buf, state, config, topo, stage)

    stage = LCS.RKStage2()

    @profile backend "flow/synchalo/2" Topologies.synchalo!(
        Topologies.FullHalo(config.grid),
        (buf.flow.fields.Us2..., buf.flow.fields.P),
        buf.halo,
        config.grid,
        backend,
        topo,
    )

    @profile backend "flow/rhs/2" f_rhs2 = Flows.rhs!(buf, state, config, topo, stage)
    @profile backend "flow/update/2" f_update2 = Flows.update!(buf, state, config, topo, stage)
    @profile backend "flow/forcing" f_forcing = Flows.forcing!(buf, state, config, topo, stage)
    @profile backend "flow/couple/2" f_couple2 = Flows.couple!(buf, state, config, topo, stage)

    evolve!(state)

    (;
        flow=(;
            rhs1=f_rhs1,
            update1=f_update1,
            couple1=f_couple1,
            rhs2=f_rhs2,
            update2=f_update2,
            forcing=f_forcing,
            couple2=f_couple2,
        )
    )
end

function evolve!(
    ::FlowParticleMode,
    buf::LCS.Buffer,
    state::LCS.State,
    config::LCS.Config,
    topo::Topologies.Topology;
    cb_between_rk_stages=nothing,
)
    backend = config.backend
    stage = LCS.RKStage1()

    @profile backend "flow/synchalo/1" Topologies.synchalo!(
        Topologies.FullHalo(config.grid),
        (buf.flow.fields.Us..., buf.flow.fields.P),
        buf.halo,
        config.grid,
        backend,
        topo,
    )

    @profile backend "particles/boundary/1" p_boundary1 = ntuple(
        i -> Particles.boundary!(buf, state, config, topo, stage, i), length(buf.particles)
    )
    @profile backend "particles/motion/1" p_motion1 = ntuple(
        i -> Particles.motion!(buf, state, config, topo, stage, i), length(buf.particles)
    )
    @profile backend "particles/update/1" p_update1 = ntuple(
        i -> Particles.update!(buf, state, config, topo, stage, i), length(buf.particles)
    )

    @profile backend "flow/rhs/1" f_rhs1 = Flows.rhs!(buf, state, config, topo, stage)
    @profile backend "flow/update/1" f_update1 = Flows.update!(buf, state, config, topo, stage)
    @profile backend "flow/couple/1" f_couple1 = Flows.couple!(buf, state, config, topo, stage)

    !isnothing(cb_between_rk_stages) && cb_between_rk_stages()

    stage = LCS.RKStage2()

    @profile backend "flow/synchalo/2" Topologies.synchalo!(
        Topologies.FullHalo(config.grid),
        (buf.flow.fields.Us2..., buf.flow.fields.P),
        buf.halo,
        config.grid,
        backend,
        topo,
    )

    @profile backend "particles/boundary/2" p_boundary2 = ntuple(
        i -> Particles.boundary!(buf, state, config, topo, stage, i), length(buf.particles)
    )
    @profile backend "particles/motion/2" p_motion2 = ntuple(
        i -> Particles.motion!(buf, state, config, topo, stage, i), length(buf.particles)
    )
    @profile backend "particles/update/2" p_update2 = ntuple(
        i -> Particles.update!(buf, state, config, topo, stage, i), length(buf.particles)
    )

    @profile backend "flow/rhs/2" f_rhs2 = Flows.rhs!(buf, state, config, topo, stage)
    @profile backend "flow/update/2" f_update2 = Flows.update!(buf, state, config, topo, stage)
    @profile backend "flow/forcing" f_forcing = Flows.forcing!(buf, state, config, topo, stage)
    @profile backend "flow/couple/2" f_couple2 = Flows.couple!(buf, state, config, topo, stage)

    evolve!(state)

    (;
        particles=(;
            boundary1=p_boundary1,
            motion1=p_motion1,
            update1=p_update1,
            boundary2=p_boundary2,
            motion2=p_motion2,
            update2=p_update2,
        ),
        flow=(;
            rhs1=f_rhs1,
            update1=f_update1,
            couple1=f_couple1,
            rhs2=f_rhs2,
            update2=f_update2,
            forcing=f_forcing,
            couple2=f_couple2,
        ),
    )
end
