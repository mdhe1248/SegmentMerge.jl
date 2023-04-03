module SegmentMerge
using Images, ImageView
using Distances, Statistics
using CSV, DataFrames
using PyPlot

# Write your package code here.
export Punctum, Segment, SegmentMerged, Centroid, get_size, find_centroid, get_amplitude, get_Segment, find_dist, findmin_paireddistance, mergeSegments, _mergeSegments, get_fieldval, get_intensity, Punctum, dist_nearest_lbl, find_nearest_lbl, find_nearest_lbls, find_cellID, get_puncta, count_puncta, merge_segimgs

include("merge_segments.jl")
include("merge_blobs.jl")

end
