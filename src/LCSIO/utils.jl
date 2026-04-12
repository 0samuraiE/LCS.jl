function require_outdir(config::LCS.Config)
    outdir = config.outdir
    isdir(outdir) || throw(ArgumentError("output directory $outdir does not exist"))
    outdir
end

function oninterval(f, state::LCS.State, interval::Integer)
    if rem(state.step, interval) == 0
        f()
    end
end

function cleanrestart(state::LCS.State, config::LCS.Config, topo::Topologies.Topology)
    (; interval_restart, num_restart_files) = config.simulate
    outdir = require_outdir(config)

    iold = state.step - interval_restart * num_restart_files
    old_restart_file = joinpath(outdir, Printf.format(FILE_RESTART_FMT, iold))
    if Topologies.isroot(topo) && isfile(old_restart_file)
        rm(old_restart_file)
    end
end
