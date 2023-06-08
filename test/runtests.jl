using Test
using NextGenSeqUtils
using RobustAmpliconDenoising

@testset "RAD" begin
    nl43env_2_chars = collect(nl43env)
    mut_positions = [729, 1729]
    for mut_pos in mut_positions
        nl43env_2_chars[mut_pos] = rand([x for x in "ACGT" if x != nl43env_2_chars[mut_pos]])
    end
    nl43env_2 = join(nl43env_2_chars)

    seqs, _ = denoise(vcat(env_pb_seq_sim(nl43env, 1000), env_pb_seq_sim(nl43env_2, 1000)))
    @test Set(seqs) == Set([nl43env, nl43env_2])
end
