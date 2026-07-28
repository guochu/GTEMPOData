using DelimitedFiles
using Printf

function diff(itr)
    A = readdlm(@sprintf("Guu10_iter%d.dat", itr))
    B = readdlm(@sprintf("Guu10_iter%d_ctqmc.dat", itr))
    C = A[:,2] - B[1:40:size(B)[1],2]
    writedlm(@sprintf("Guu10_iter%d_diff.dat", itr), [A[:,1] C])
end

diff(1)
diff(3)
diff(8)
diff(9)
diff(10)

function diff(itr)
    A = readdlm(@sprintf("Gud10_iter%d.dat", itr))
    B = readdlm(@sprintf("Gud10_iter%d_ctqmc.dat", itr))
    C = A[:,2] - B[1:40:size(B)[1],2]
    writedlm(@sprintf("Gud10_iter%d_diff.dat", itr), [A[:,1] C])
end
diff(1)
diff(3)
diff(8)
diff(9)
diff(10)

function diff(itr)
    A = readdlm(@sprintf("Guu100_iter%d.dat", itr))
    B = readdlm(@sprintf("Guu100_iter%d_ctqmc.dat", itr))
    C = A[:,2] - B[1:20:size(B)[1],2]
    writedlm(@sprintf("Guu100_iter%d_diff.dat", itr), [A[:,1] C])
end

diff(1)
diff(10)

function diff(itr)
    A = readdlm(@sprintf("Gud100_iter%d.dat", itr))
    B = readdlm(@sprintf("Gud100_iter%d_ctqmc.dat", itr))
    C = A[:,2] - B[1:20:size(B)[1],2]
    writedlm(@sprintf("Gud100_iter%d_diff.dat", itr), [A[:,1] C])
end

diff(1)
diff(10)
