"""
    RKStage

Abstract type for Runge-Kutta stages.
"""
abstract type RKStage end

"""
    RKStage1 <: RKStage

First stage of Runge-Kutta method.
"""
struct RKStage1 <: RKStage end

"""
    RKStage2 <: RKStage

Second stage of Runge-Kutta method.
"""
struct RKStage2 <: RKStage end
