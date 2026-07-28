using DelimitedFiles
using Printf

function diff(chi)
    A = readdlm(@sprintf("vacuum_Nt100_chi%d.dat", chi))
    B = readdlm(@sprintf("vacuum_Nt200_chi%d.dat", chi))
    C = A[:,2] - B[1:2:200,2]
    writedlm(@sprintf("vacuum_diff_chi%d.dat", chi), [A[:,1] C])

    A = readdlm(@sprintf("nup_Nt100_chi%d.dat", chi))
    B = readdlm(@sprintf("nup_Nt200_chi%d.dat", chi))
    C = A[:,2] - B[1:2:200,2]
    writedlm(@sprintf("nup_diff_chi%d.dat", chi), [A[:,1] C])

    A = readdlm(@sprintf("ndown_Nt100_chi%d.dat", chi))
    B = readdlm(@sprintf("ndown_Nt200_chi%d.dat", chi))
    C = A[:,2] - B[1:2:200,2]
    writedlm(@sprintf("ndown_diff_chi%d.dat", chi), [A[:,1] C])

    A = readdlm(@sprintf("nn_Nt100_chi%d.dat", chi))
    B = readdlm(@sprintf("nn_Nt200_chi%d.dat", chi))
    C = A[:,2] - B[1:2:200,2]
    writedlm(@sprintf("nn_diff_chi%d.dat", chi), [A[:,1] C])
end

diff(200)
diff(250)
diff(300)
diff(400)
