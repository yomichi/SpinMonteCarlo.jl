using Test
using SpinMonteCarlo

@testset "unify!" begin
    uf = UnionFind(0)
    addnode!(uf)  # node 1
    addnode!(uf)  # node 2
    addnode!(uf)  # node 3
    addnode!(uf)  # node 4

    # merge two pairs
    r12 = unify!(uf, 1, 2)
    r34 = unify!(uf, 3, 4)
    @test uf.weights[r12] == 2
    @test uf.weights[r34] == 2
    @test uf.nclusters == 2

    # merge the resulting sets
    r = unify!(uf, 1, 3)
    @test r == SpinMonteCarlo.root!(uf, 2) == SpinMonteCarlo.root!(uf, 4)
    @test uf.weights[r] == 4
    @test uf.nclusters == 1
end
