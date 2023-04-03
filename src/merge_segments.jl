struct Segment
  size::Int #npixel
  centroid::Tuple{Int, Int}
  amplitude::Float64
  segmentID::Int
  channel::Int
end

struct SegmentMerged
  segments::Vector{Segment}
  cellID::Int
end

struct Centroid
  segmentID::Int
  centroid::Tuple{Int, Int}
end

function get_size(segimg, lbl)
  return(sum(segimg .== lbl))
end

function get_size(segimg, lbls::Vector)
  output = Vector{Tuple{Int, Int}}(undef,length(lbls))
  for (i, lbl) in enumerate(lbls)
    output[i] = (lbl, get_size(segimg, lbl))
  end
  return(output)
end

function find_centroid(segimg, lbl)
  coords = findall(segimg .== lbl)
  xi, yi = 0, 0
  for coord in coords
    xi = xi+coord[1]
    yi = yi+coord[2]
  end
  centroid_xi = round(Int, xi/length(coords))
  centroid_yi = round(Int, yi/length(coords))
  return(Centroid(lbl, (centroid_xi, centroid_yi)))
end

function find_centroid(segimg, lbls::Vector)
  output = Vector{Centroid}(undef,length(lbls))
  for (i, lbl) in enumerate(lbls)
    output[i] = find_centroid(segimg, lbl)
  end
  return(output)
end

function find_centroid(img, segimg, lbl)
  lblimg = copy(img)
  coords = findall(segimg .== lbl)
  xi, yi = 0, 0
  for coord in coords
    xi = xi+coord[1]*img[coord]
    yi = yi+coord[2]*img[coord]
  end
  centroid_xi = round(Int, xi/length(coords))
  centroid_yi = round(Int, yi/length(coords))
  return(Centroid(lbl, (centroid_xi, centroid_yi)))
end

function find_centroid(img, segimg, lbls::Vector)
  output = Vector{Centroid}(undef,length(lbls))
  for (i, lbl) in enumerate(lbls)
    output[i] = find_centroid(img, segimg, lbl)
  end
  return(output)
end

"""
`thresh` is either true or specific value for thresholding before mean calculation.
"""
function get_amplitude(img2d, segimg, lbl; thresh = true)
  pxinten = img2d[segimg .== lbl]
  if thresh 
    thresh = mean(pxinten)
  end
  amp = mean(pxinten[pxinten .> thresh])
  return(amp)
end
function get_amplitude(img2d, segimg, lbls::Vector; thresh = true)   
  output = Vector{Tuple{Int, Float64}}(undef, length(lbls))
  for (i, lbl) in enumerate(lbls)
    output[i] = (lbl, get_amplitude(img2d, segimg, lbl; thresh = thresh))
  end
  return(output)
end

function get_Segment(img2d, segimg, lbl, channel)
  Segment(
    get_size(segimg, lbl),
    find_centroid(segimg, lbl).centroid,
    get_amplitude(img2d, segimg, lbl),
    lbl,
    channel)
end

function find_dist(pos1::Vector{Centroid}, pos2::Vector{Centroid})
  pos1 = getfield.(pos1, :centroid)
  pos2 = getfield.(pos2, :centroid)
  dist = find_dist(pos1, pos2)
end

function find_dist(pos1, pos2)
  dist = zeros(length(pos1), length(pos2))
  for (i, p1) in enumerate(pos1), (j, p2) in enumerate(pos2)
    dist[i, j] = evaluate(Euclidean(), p1, p2)
  end
  dist
end

function findmin_paireddistance(dist, dist_thresh)
  d1, p1 = vec.(findmin(dist, dims = 2))
  d2, p2 = vec.(findmin(dist, dims = 1))
  coords = Vector{CartesianIndex}();
  for pp1 in p1, pp2 in p2
    if pp1 == pp2
      push!(coords, pp1)
    end
  end 
  d = [dist[x] for x in coords] #map(x -> dist[x], coords) #distance
  keep = d .< dist_thresh #thresholding
  return(coords[keep], d[keep])
end

function findmin_paireddistance(centoirds1::Vector{Centroid}, centroids2::Vector{Centroid}, dist_thresh::Number)
  dist = find_dist(centroids1, centroids2)
  d1, p1 = vec.(findmin(dist, dims = 2))
  d2, p2 = vec.(findmin(dist, dims = 1))
  centroid_pairs = Vector{Pair{Centroid, Centroid}}();
  d = Vector{Float64}();
  for (i, pp1) in enumerate(p1), (j, pp2) in enumerate(p2)
    if pp1 == pp2
      push!(centroid_pairs, Pair(centroids1[i], centroids2[j]))
      push!(d, d1[i])
    end
  end
  keep = d .< dist_thresh #thresholding
  return(centroid_pairs[keep], d[keep])
end

function mergeSegments(s1::Vector{Segment}, s2::Vector{Segment}, dist_thresh::Number)
  centroids1 = [Centroid(s.segmentID, s.centroid) for s in s1]
  centroids2 = [Centroid(s.segmentID, s.centroid) for s in s2]
  centroidpairs, centroidpair_dists = findmin_paireddistance(centroids1, centroids2, dist_thresh)
  mergedsegments = Vector{SegmentMerged}();
  idx = 1
  c1 = 0
  c2 = 0
  nrounds = 0
  # push non-merged seg1
  for seg1 in s1
    c1 = Centroid(seg1.segmentID, seg1.centroid)
    if !(c1 ∈  getfield.(centroidpairs, :first)) #map(x -> x.first, centroidpairs))
      push!(mergedsegments, _mergeSegments(seg1, idx))
      idx += 1
    end
  end
  # push non-merged seg2
  for seg2 in s2
    c2 = Centroid(seg2.segmentID, seg2.centroid) 
    if !(c2 ∈  getfield.(centroidpairs, :second)) #map(x -> x.second, centroidpairs))
      push!(mergedsegments, _mergeSegments(seg2, idx))
      idx += 1
    end
  end
  # push merged segments
  for seg2 in s2
    for seg1 in s1
      c1 = Centroid(seg1.segmentID, seg1.centroid)
      c2 = Centroid(seg2.segmentID, seg2.centroid) 
      if Pair(c1, c2) ∈  centroidpairs
        push!(mergedsegments, _mergeSegments(seg1, seg2, idx))  #push merged segments
        idx += 1
      end
    end
  end
  return(mergedsegments)
end

function _mergeSegments(s1::Segment, s2::Segment, idx)
  SegmentMerged([s1, s2], idx)
end

function _mergeSegments(s1::Segment, idx)
  SegmentMerged([s1], idx)
end

function get_fieldval(segment::Segment, fieldname::Symbol, channel; otherwise_return_val = 0.0)
  if segment.channel == channel 
    return(getproperty(segment, fieldname))
  else
    return(otherwise_return_val)
  end
end 

function get_intensity(mergedsegments::Vector{SegmentMerged}, channel)
  inten = zeros(length(mergedsegments))
  for (i, mergedsegment) in enumerate(mergedsegments)
    inten[i] = sum(get_fieldval.(mergedsegment.segments, :amplitude, channel))
  end
  inten
end

