struct Punctum
  location::Tuple{Int, Int}
  amplitude::Float64
  cellID::Int
end

function dist_nearest_lbl(segimg, blobpos::Vector{CartesianIndex{2}})
  blobpostuple = getfield.(blobpos, :I)
  pxpos, dist = dist_nearest_lbl(segimg, blobpostuple)
  return(pxpos, dist)
end

function dist_nearest_lbl(segimg, blobpos::Vector{Tuple{Int, Int}})
  lbls = unique(segimg)
  lbls = lbls[lbls .!= 0]
  pxpos = findall(segimg .!= 0)
  dist = zeros(length(pxpos), length(blobpos))
  for i in eachindex(pxpos), j in eachindex(blobpos)
	  dist[i, j] = evaluate(Euclidean(), pxpos[i].I, blobpos[j]) 
  end
  return(pxpos, dist)
end

function find_nearest_lbl(segimg1::AbstractMatrix, blobpos::AbstractVector)
  lbl_pxpos1, lbl_dist1 = dist_nearest_lbl(segimg1, blobpos)
  mn_dist1, mn_dist1_idx = vec.(findmin(lbl_dist1, dims = 1))
  mnidx1 = getfield.(mn_dist1_idx, :I) #map(x -> x.I, mn_dist1_idx) #pixel pos and blobID
  #mnidx1 = map(x -> x.I, mn_dist1_idx) #pixel pos and blobID
  coord1 = [lbl_pxpos1[x[1]] for x in mnidx1] #map(x -> lbl_pxpos1[x[1]], mnidx1)
  nearest_lbls1 = segimg1[coord1]
  return(nearest_lbls1, mn_dist1)
end

function find_nearest_lbls(segimgs::Tuple, blobpos::AbstractVector)
  nearest_lbls = Vector{Tuple{AbstractVector, AbstractVector}}()
  for segimg in segimgs
    lbls, dists = find_nearest_lbl(segimg, blobpos)
    push!(nearest_lbls, (lbls, dists)) 
  end
  dists = hcat([x[2] for x in nearest_lbls]...) #map(x -> x[2], nearest_lbls)...)
  mndists, i = findmin(dists, dims = 2)  
  lbls = hcat([x[1] for x in nearest_lbls]...) #map(x -> x[1], nearest_lbls)...)
  channel = [x[2] for x in i] #map(x -> x[2], i)
  mnidx = [CartesianIndex(i, channel[i]) for i in eachindex(blobpos)]
  nearest_lbls = lbls[mnidx]
  return(mndists,  nearest_lbls, channel)
end

function find_cellID(mergedsegments, segmentID, channel)
  for mergedsegment in mergedsegments
    segments = mergedsegment.segments
    for segment in segments
      if segment.channel == channel && segment.segmentID == segmentID
        return(mergedsegment.cellID)
      end
    end
  end
end

function get_puncta(mergedsegments, segimgs, blobs, dist_thresh)
  blobpos = getfield.(blobs, :location) #x, y changed here
  blobamp = getfield.(blobs, :amplitude) #x, y changed here
  mn_dists, nearest_segmentID, channels = find_nearest_lbls(segimgs, blobpos)
  cellids = zeros(Int, length(blobs))
  puncta = Vector{Punctum}(undef, length(blobs));
  for i in eachindex(blobs)
    if mn_dists[i] < dist_thresh
      cellids[i] = find_cellID(mergedsegments, nearest_segmentID[i], channels[i])
    end
  end 
  for i in eachindex(blobs)
    puncta[i] = Punctum(blobpos[i], blobamp[i], cellids[i]) 
  end
 return(puncta)
end

## DataFrame intensity, puncta counts
function count_puncta(mergedsegments, puncta)
  cellids = getfield.(mergedsegments, :cellID)
  puncta_cellids = getfield.(puncta, :cellID)
  counts = zeros(Int, length(cellids))
  for idx in puncta_cellids
    countidx = findfirst(cellids .== idx)
    if !isnothing(countidx)
      counts[countidx] += 1
    end
  end
  return(counts)
end

"""
merge segmented images based on `mergedsegments`
segmentIDs will be the same as cellID
"""
function merge_segimgs(segimgs, mergedsegments)
  segimg = zeros(eltype(segimgs[1]), size(segimgs[1]))
  for ch in eachindex(segimgs)
    segids = []
    for mergedsegment in mergedsegments
      for segments in mergedsegment.segments
        segid = get_fieldval(segments, :segmentID, ch; otherwise_return_val = nothing)
        if !isnothing(segid)
          push!(segids, segid => mergedsegment.cellID)
        end
      end
    end
    [segimg[segimgs[ch] .== segid.first] .= segid.second for segid in segids]
    end
  return(segimg)
end
