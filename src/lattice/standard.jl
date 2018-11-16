const stdbravais = Dict{String, LatticeParameter}()
const stdunitcells = Dict{String, LatticeParameter}()
const stdlattices = Dict{String, LatticeParameter}()

## 1D bravais
stdbravais["chain"] = P(:name =>"chain",
                        :dimension => 1,
                        :parameters => [(:a,1.0)],
                        :basis => :(reshape([a],1,1)))

## 2D bravais
stdbravais["orthorhombic2d"] = P(:name => "orthorhombic2d",
                                 :dimension => 2,
                                 :parameters => [(:a,1.0), (:b,1.0)],
                                 :basis => :([a 0 
                                              0 b]))
stdbravais["hexagonal2d"] = P(:name => "hexagonal2d",
                              :dimension => 2,
                              :parameters => [(:a,1.0)],
                              :basis => :([a  -a/2
                                           0  a*sqrt(3.0)/2]))
stdbravais["monoclinic2d"] = P(:name => "monoclinic2d",
                               :dimension => 2,
                               :parameters => [(:a,1.0), (:b,1.0), (:theta, 90.0)],
                               :basis => :([a b*cosd(theta)
                                            0 b*sind(theta)]))

## 3D bravais
stdbravais["orthorhombic3d"] = P(:name => "orthorhombic3d",
                                 :dimension => 3,
                                 :parameters => [(:a,1.0),(:b,1.0),(:c,1.0)],
                                 :basis => :([a 0 0
                                              0 b 0
                                              0 0 c]))

## 1D unitcells
stdunitcells["simple1d"] = P(:name => "simple1d",
                             :dimension => 1,
                             :sites => [P(:id=>1,
                                          :sitetype=>1,
                                          :coord=>[0.0]
                                         )],
                             :bonds => [P(:bondtype=>1,
                                          :source=>P(:id=>1, :offset=>[0]),
                                          :target=>P(:id=>1, :offset=>[1]),
                                         )])

stdunitcells["bond-alternating simple1d"] = P(:name => "bond-alternating simple1d",
                                              :dimension => 1,
                                              :sites => [P(:id=>1,
                                                           :sitetype=>1,
                                                           :coord=>[0.0],
                                                          ),
                                                         P(:id=>2,
                                                           :sitetype=>2,
                                                           :coord=>[0.5],
                                                          ),
                                                        ],
                                              :bonds => [P(:bondtype=>1,
                                                           :source=>P(:id=>1, :offset=>[0]),
                                                           :target=>P(:id=>2, :offset=>[0]),
                                                          ),
                                                         P(:bondtype=>2,
                                                           :source=>P(:id=>2, :offset=>[0]),
                                                           :target=>P(:id=>1, :offset=>[1]),
                                                           )])

## 2D unitcells
stdunitcells["simple2d"] = P(:name => "simple2d",
                             :dimension => 2,
                             :sites => [P(:id=>1,
                                          :sitetype=>1,
                                          :coord=>[0.0, 0.0],
                                         )],
                             :bonds => [P(:bondtype=>1,
                                          :source=>P(:id=>1, :offset=>[0,0]),
                                          :target=>P(:id=>1, :offset=>[1,0]),),
                                        P(:bondtype=>1,
                                          :source=>P(:id=>1, :offset=>[0,0]),
                                          :target=>P(:id=>1, :offset=>[0,1]),),
                                       ])
stdunitcells["triangular cell"] = P(:name => "triangular cell",
                                    :dimension => 2,
                                    :sites => [P(:id=>1,
                                                 :sitetype=>1,
                                                 :coord=>[0.0,0.0],
                                                )],
                                    :bonds => [P(:bondtype=>1,
                                                 :source=>P(:id=>1, :offset=>[0,0]),
                                                 :target=>P(:id=>1, :offset=>[1,0]),),
                                               P(:bondtype=>1,
                                                 :source=>P(:id=>1, :offset=>[0,0]),
                                                 :target=>P(:id=>1, :offset=>[0,1]),),
                                               P(:bondtype=>1,
                                                 :source=>P(:id=>1, :offset=>[0,0]),
                                                 :target=>P(:id=>1, :offset=>[1,1]),),
                                              ])
stdunitcells["honeycomb cell"] = P(:name => "honeycomb cell",
                                   :dimension => 2,
                                   :sites => [P(:id=>1,
                                                :sitetype=>1,
                                                :coord=>[0.0, 0.0],
                                               ),
                                              P(:id=>2,
                                                :sitetype=>2,
                                                :coord=> (2//3).*[0.5, 1.0],
                                               )],
                                   :bonds => [P(:bondtype=>1,
                                                :source=>P(:id=>1, :offset=>[0,0]),
                                                :target=>P(:id=>2, :offset=>[0,0]),),
                                              P(:bondtype=>2,
                                                :source=>P(:id=>2, :offset=>[0,0]),
                                                :target=>P(:id=>1, :offset=>[1,1]),),
                                              P(:bondtype=>3,
                                                :source=>P(:id=>2, :offset=>[0,0]),
                                                :target=>P(:id=>1, :offset=>[0,1]),),
                                             ])



## 3D unitcells
stdunitcells["simple3d"] = P(:name => "simple3d",
                             :dimension => 3,
                             :sites => [P(:id=>1,
                                          :sitetype=>1,
                                          :coord=>[0.0, 0.0, 0.0],
                                         )],
                             :bonds => [P(:bondtype=>1,
                                          :source=>P(:id=>1, :offset=>[0,0,0]),
                                          :target=>P(:id=>1, :offset=>[1,0,0]),),
                                        P(:bondtype=>1,
                                          :source=>P(:id=>1, :offset=>[0,0,0]),
                                          :target=>P(:id=>1, :offset=>[0,1,0]),),
                                        P(:bondtype=>1,
                                          :source=>P(:id=>1, :offset=>[0,0,0]),
                                          :target=>P(:id=>1, :offset=>[0,0,1]),),
                                       ])

## 1D lattices
stdlattices["chain lattice"] = P(:name => "chain lattice",
                                 :dimension => 1,
                                 :bravais => "chain",
                                 :unitcell => "simple1d",
                                 :parameters => [],
                                 :periodic => [true],
                                )

stdlattices["bond-alternating chain lattice"] = P(:name => "bond-alternating chain lattice",
                                 :dimension => 1,
                                 :bravais => "chain",
                                 :unitcell => "bond-alternating simple1d",
                                 :parameters => [(:a, 2.0)],
                                 :periodic => [true],
                                )


## 2D lattices
stdlattices["square lattice"] = P(:name => "square lattice",
                                  :dimension => 2,
                                  :bravais => "orthorhombic2d",
                                  :unitcell => "simple2d",
                                  :parameters => [],
                                  :periodic => [true, true],
                                 )
stdlattices["triangular lattice"] = P(:name => "triangular lattice",
                                      :dimension => 2,
                                      :bravais => "hexagonal2d",
                                      :unitcell => "triangular cell",
                                      :parameters => [],
                                      :periodic => [true, true],
                                     )
stdlattices["honeycomb lattice"] = P(:name => "honeycomb lattice",
                                     :dimension => 2,
                                     :bravais => "hexagonal2d",
                                     :unitcell => "honeycomb cell",
                                     :parameters => [(:a, sqrt(3.0))],
                                     :periodic => [true, true],
                                    )
stdlattices["ladder"] = P(:name => "ladder",
                          :dimension => 2,
                          :bravais => "orthorhombic2d",
                          :unitcell => "simple2d",
                          :parameters => [],
                          :periodic => [true, false],
                         )

## 3D lattices
stdlattices["cubic lattice"] = P(:name => "cubic lattice",
                                 :dimension => 3,
                                 :bravais => "orthorhombic3d",
                                 :unitcell => "simple3d",
                                 :periodic => [true, true, true],
                                )

## special lattices
stdlattices["fully connected graph"] = P(:name => "fully connected graph",
                                           :dim => 1,
                                           :generator => generate_fully_connected_graph,
                                          )
