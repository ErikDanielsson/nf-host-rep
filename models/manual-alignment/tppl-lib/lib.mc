include "seq.mc"
include "common.mc"
include "float.mc"

let cons = cons 
let absf = absf
let mulf = mulf
let isNaN = isNaN

let nestList : all a. [Int] -> Int -> Int -> [[Int]] = lam l. lam r. lam c.
  let row = lam i.
    let rs = addi (muli i c) 1 in
    let re = addi (muli (addi i 1) c) 1 in
    slice l rs re in
  map row (range 0 r 1)