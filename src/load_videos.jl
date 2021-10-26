using HDF5
using NPZ
#%% md
# Load all the matrices. The matrices are exported from the recorded videos and saved to a folder.
#Each matrix has two spatial dimensions (the screen) and a temporal dimension, the first two reflect the pixels of the image, the third corresponds to the video frames.
# Our video has 250x420 pixels and 15 f/s
##

data_path = "/run/media/cocconat/data/ownCloud_MPI/granulari"
(root,dirs,files) =  first(walkdir(data_path))
speed_folders = filter(dir->startswith(dir,"V"),dirs)
for dir in speed_folders
    file = joinpath(root,dir, "image_matrix.npy")
    speed = parse(Float64,dir[3:end])
    try
        _data = npzread(file)
        target = joinpath(root,dir*".h5")
        isfile(target) && rm(target)
        h5open(target, "w") do fid
            fid["matrix"]=_data
            fid["speed"] = speed
        end
        println("loaded:", speed)
    catch
        println("failed:", speed)
    end
end
