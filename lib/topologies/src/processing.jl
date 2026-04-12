"""
    Processing

Abstract type for processing modes in distributed computing.
"""
abstract type Processing end

"""
    SingleProcessing <: Processing

Single-process processing mode for non-distributed execution.
"""
struct SingleProcessing <: Processing end

"""
    MultiProcessing <: Processing

Multi-process processing mode for distributed execution.
"""
struct MultiProcessing <: Processing end

"""
    processing(::Topology)

Determine processing mode based on topology dimensions.
"""
processing(::Topology{1,1,1}) = SingleProcessing()
processing(::Topology) = MultiProcessing()

"""
    is_multi_processing(topo::Topology)

Check if topology uses multi-processing mode.
"""
is_multi_processing(topo::Topology) = is_multi_processing(processing(topo))
is_multi_processing(::SingleProcessing) = false
is_multi_processing(::MultiProcessing) = true
