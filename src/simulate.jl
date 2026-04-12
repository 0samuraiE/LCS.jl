struct Callback
    interval :: Int
    func     :: Function
end

"""
    simulate(file)
    simulate(config)
    simulate(buffer, state, config, topo)

Run complete LCS simulation from configuration or existing objects.
"""
function simulate(
    file::AbstractString;
    callbacks::Vector{Callback}=Callback[],
    statfilter::Base.Callable=identity,
    logfilter::Base.Callable=identity,
)
    config = LCSIO.readconfig(file)
    simulate(config; callbacks, statfilter, logfilter)
end

function simulate(
    config::AbstractConfig;
    callbacks::Vector{Callback}=Callback[],
    statfilter::Base.Callable=identity,
    logfilter::Base.Callable=identity,
)
    backend = config.backend
    @log backend "setup/topology" topo = Topologies.Topology(config)
    @log backend "setup/device" Topologies.device!(config)
    @log backend "setup/buffer" buffer = LCS.Buffer(config, topo)
    @log backend "setup/state" state = LCS.State(config, topo)

    simulate(buffer, state, config, topo; callbacks, statfilter, logfilter)
end

function simulate(
    buffer::AbstractBuffer,
    state::AbstractState,
    config::AbstractConfig,
    topo::Topologies.Topology;
    callbacks::Vector{Callback}=Callback[],
    statfilter::Base.Callable=identity,
    logfilter::Base.Callable=identity,
)
    config = LCSIO.init(config, topo)
    backend = config.backend
    @profile backend "init" LCS.init!(buffer, state, config, topo)

    _should_print = get_log_level() != LCS_LOG_QUIET && Topologies.isroot(topo)

    GC.enable(false)
    for _ in (state.step):(config.simulate.last_step)
        @profile backend "timestep" LCS.timestep!(buffer, state, config, topo)
        @log backend "evolve" log = LCS.evolve!(buffer, state, config, topo)

        _should_print && LCS.printstate(state)
        _should_print && LCS.printlog(log; filter=logfilter)

        LCSIO.oninterval(state, config.simulate.interval_stat) do
            @profile backend "stat" stat = LCS.stat(buffer, state, config, topo)
            _should_print && LCS.printstat(stat; filter=statfilter)
            @profile backend "io/savestat" LCSIO.savestat(config, state, topo; stat, log)
        end

        LCSIO.oninterval(state, config.simulate.interval_restart) do
            @profile backend "io/saverestart" LCSIO.saverestart(buffer, state, config, topo)
            LCSIO.cleanrestart(state, config, topo)
        end

        for cb in callbacks
            LCSIO.oninterval(state, cb.interval) do
                cb.func(buffer, state, config, topo; log)
            end
        end

        LCSIO.oninterval(state, config.simulate.interval_gc) do
            Topologies.barrier(topo)
            GC.enable(true)
            GC.gc()
            GC.gc()
            GC.gc()
            GC.enable(false)
        end
    end
    GC.enable(true)

    buffer, state, config, topo
end
