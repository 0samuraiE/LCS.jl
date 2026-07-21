"""
    @profile backend msg ex

Execute `ex`, printing `msg` with elapsed time when log level is `PROFILE`.
"""
macro profile(backend, msg, ex)
    quote
        if get_log_level() == LCS_LOG_PROFILE
            mpiprintln($(esc(msg)), ": starting...")
            t0 = time()
            local x = $(esc(ex))
            KA.synchronize($(esc(backend)))
            elapsed = time() - t0
            mpiprintln($(esc(msg)), ": done. (", _fmt_elapsed(elapsed), ")")
            x
        else
            $(esc(ex))
        end
    end
end

"""
    @log backend msg ex

Execute `ex`, printing `msg` with elapsed time unless log level is `QUIET`.
"""
macro log(backend, msg, ex)
    quote
        if get_log_level() == LCS_LOG_QUIET
            $(esc(ex))
        else
            mpiprintln($(esc(msg)), ": starting...")
            KA.synchronize($(esc(backend)))
            t0 = time()
            local x = $(esc(ex))
            KA.synchronize($(esc(backend)))
            elapsed = time() - t0
            mpiprintln($(esc(msg)), ": done. (", _fmt_elapsed(elapsed), ")")
            x
        end
    end
end

function _fmt_elapsed(s::Real)
    s < 1E-3 && return "< 1ms"
    s < 1 && return @sprintf("%.0fms", s * 1000)
    s < 60 && return @sprintf("%.2fs", s)
    m = floor(Int, s / 60)
    @sprintf("%dm%ds", m, floor(Int, s) % 60)
end
