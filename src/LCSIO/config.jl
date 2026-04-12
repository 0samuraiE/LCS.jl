"""
    readconfig(file; sync=MPI.Initialized(), patch=identity)

Read and parse simulation configuration from YAML.
"""
function readconfig(file::AbstractString; sync::Bool=MPI.Initialized(), patch::Base.Callable=identity)
    if sync
        comm = Topologies.comm()
        config = if MPI.Comm_rank(comm) == 0
            _readconfig(file, patch)
        else
            nothing
        end
        MPI.bcast(config, 0, comm)
    else
        _readconfig(file, patch)
    end
end

function _readconfig(file::AbstractString, patch::Base.Callable)
    dict = YAML.load_file(file; dicttype=Dict{String,Any})
    "resume" in keys(dict) && throw(
        ArgumentError(
            "resume must not be set in configuration, resume mode is controlled by LCS_RESUME environment variable"
        ),
    )

    dict["resume"] = false
    config = PolySerde.deserialize(LCS.Config, dict)
    patch(config)
end

function isresumable(config::LCS.Config)
    isdir(config.outdir) || return false
    restart_files = find_restart_files(config.outdir)
    !isempty(restart_files) && isfile(joinpath(config.outdir, FILE_CONFIG))
end

"""
    init(config, topo)

Resolve resume state and initialize output directory.
"""
function init(config::LCS.Config, topo::Topologies.Topology)
    config = if Topologies.isroot(topo)
        if isdir(config.outdir) && !Base.get_bool_env(LCS.ENV_LCS_RESUME, false)
            throw(ArgumentError("output directory $(config.outdir) already exists, set LCS_RESUME=1 to resume"))
        end

        if isdir(config.outdir) && Base.get_bool_env(LCS.ENV_LCS_RESUME, false)
            if isresumable(config)
                config = @set config.resume = true
                config = patch_resume(config.mode, config)
            else
                throw(
                    ArgumentError(
                        "output directory $(config.outdir) exists but is not resumable, missing restart files or config.lcs-yaml",
                    ),
                )
            end
        end

        config
    end
    config = MPI.bcast(config, 0, Topologies.comm())

    if !config.resume
        if Topologies.isroot(topo)
            mkpath(config.outdir)
            dict = PolySerde.serialize(LCS.Config, config; dicttype=OrderedDict{String,Any})
            delete!(dict, "resume")
            YAML.write_file(joinpath(config.outdir, FILE_CONFIG), dict)
        end
        Topologies.barrier(topo)
    end
    config
end

function patch_resume(::LCS.FlowMode, config::LCS.Config)
    resume_file = find_latest_restart_file(config.outdir)
    config = @set config.flow.init = LCS.Restart(resume_file)
    config
end

function patch_resume(::LCS.FlowParticleMode, config::LCS.Config)
    resume_file = find_latest_restart_file(config.outdir)
    config = @set config.flow.init = LCS.Restart(resume_file)

    particles = if isnothing(config.particles)
        nothing
    else
        ntuple(length(config.particles)) do iprofile
            p = config.particles[iprofile]
            p = @set p.init.id = Particles.ParticleRestart(resume_file, iprofile)
            p = @set p.init.position = Particles.ParticleRestart(resume_file, iprofile)
            p = @set p.init.velocity = Particles.ParticleRestart(resume_file, iprofile)
            p = @set p.init.size = Particles.ParticleRestart(resume_file, iprofile)
            p
        end
    end

    @set config.particles = particles
end

function find_restart_files(dir::AbstractString)
    files = readdir(dir; join=true)
    filter(f -> occursin(FILE_RESTART_REGEX, f), files)
end

function find_latest_restart_file(dir::AbstractString)
    restart_files = find_restart_files(dir)
    !isempty(restart_files) || throw(ArgumentError("no restart files found in directory $dir"))
    last(sort(restart_files; by=f -> parse(Int, extract_step(f))))
end

function extract_step(restart_file::AbstractString)
    m = match(FILE_RESTART_REGEX, restart_file)
    !isnothing(m) || throw(ArgumentError("invalid restart file name $restart_file"))
    m.captures[1]
end
