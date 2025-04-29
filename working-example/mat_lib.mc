include "ext/mat-ext.mc"
include "ext/dist-ext.mc"
include "seq.mc"
include "math.mc"
include "common.mc"
include "bool.mc"
include "int.mc"

-- This mirrors the tensorFold function 
let mtxFold
  : all a. all b. (b -> a -> b) -> b -> Mat a -> b =
    lam f. lam acc. lam m. foldl f acc (extArrToSeq m.arr)
  
let seqKroneckerDelta : Int -> Int -> [Float] = lam i. lam n.
  let k = subi i 1 in
  let delta = compose bool2real (eqi k) in
  map delta (range 0 n 1)

-- This creates a row vector with a single one
let mtxColKroneckerDelta : Int -> Int -> Mat Float = lam i. lam n. 
  mtxCreate n 1 (seqKroneckerDelta i n)

let mtxRowKroneckerDelta : Int -> Int -> Mat Float = lam i. lam n. 
  mtxCreate 1 n (seqKroneckerDelta i n)

let mtxOnes : Int -> Int -> Mat Float = lam n. lam m. 
  mtxCreate n m (make (muli n m) 1. )
  
let minf : Float -> Float -> Float = lam a. lam b. minf a b
let maxf : Float -> Float -> Float = lam a. lam b. maxf a b
let maxi : Int -> Int -> Int = lam a. lam b. maxi a b
let mini : Int -> Int -> Int = lam a. lam b. mini a b

let pow = lam x. lam y. pow x y

let floor = floorfi