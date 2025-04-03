include "matrix.mc"
include "ext/matrix-ext.mc"
include "ext/dist-ext.mc"
include "phylo.mc"

con Node : {age: Float, seq: [Int], left: Tree, right: Tree} -> Tree
let getAge = lam n.
    match n with Node r then r.age else
        match n with Leaf r then r.age else
            never 