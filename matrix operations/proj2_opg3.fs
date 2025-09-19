module proj2_opg3
open ProjectA
open System
open LinAlgDat.Core

// Helper functions for matrix arithmetic operations

// Matrix subtraction: A - B
let matrixSubtract (A: Matrix) (B: Matrix) =
    if A.M_Rows <> B.M_Rows || A.N_Cols <> B.N_Cols then
        failwith "Matrix dimensions must match for subtraction"
    let result = Matrix(A.M_Rows, A.N_Cols)
    for i in 0..A.M_Rows-1 do
        for j in 0..A.N_Cols-1 do
            result.[i, j] <- A.[i, j] - B.[i, j]
    result

// Matrix negation: -A
let matrixNegate (A: Matrix) =
    let result = Matrix(A.M_Rows, A.N_Cols)
    for i in 0..A.M_Rows-1 do
        for j in 0..A.N_Cols-1 do
            result.[i, j] <- -A.[i, j]
    result

// Matrix scaling: c * A
let matrixScale (c: float) (A: Matrix) =
    let result = Matrix(A.M_Rows, A.N_Cols)
    for i in 0..A.M_Rows-1 do
        for j in 0..A.N_Cols-1 do
            result.[i, j] <- c * A.[i, j]
    result

// Helper function to create a block matrix since there's no built-in Block method
let createBlockMatrix (blocks: Matrix list list) =
    let rowCount = blocks.Length
    let colCount = if rowCount > 0 then blocks.[0].Length else 0
    
    let totalRows = blocks |> List.sumBy (fun row -> 
        if row.Length > 0 then row.[0].M_Rows else 0)
    let totalCols = 
        if blocks.Length > 0 then
            blocks.[0] |> List.sumBy (fun m -> m.N_Cols)
        else 0
    
    let result = Matrix(totalRows, totalCols)
    
    // Fill in the blocks
    let mutable rowOffset = 0
    for i in 0..rowCount-1 do
        let mutable colOffset = 0
        for j in 0..colCount-1 do
            let block = blocks.[i].[j]
            for r in 0..block.M_Rows-1 do
                for c in 0..block.N_Cols-1 do
                    result.[rowOffset + r, colOffset + c] <- block.[r, c]
            colOffset <- colOffset + blocks.[i].[j].N_Cols
        rowOffset <- rowOffset + (if blocks.[i].Length > 0 then blocks.[i].[0].M_Rows else 0)
    
    result

// 3) Hjælpefunktioner til at bygge 4×4-matricer ud fra 2×2-rotation og forward
let buildRotation4 θ =
    let c = cos θ
    let s = sin θ

    let R2 = Matrix(array2D [[c; -s]; [s; c]])

    let I2 = MatrixFactory.Identity 2
    let Z2 = Matrix(2, 2)

    let Lθ = createBlockMatrix [[I2; Z2]; [matrixSubtract I2 R2; R2]]
    let Rθ = createBlockMatrix [[I2; Z2]; [matrixSubtract I2 (matrixScale (-1.0) R2); R2]]
    (Lθ, Rθ)

// 4) Build forward‐matrixen F
let buildForward4 () =
    let Z2 = Matrix(2, 2)
    let I2 = MatrixFactory.Identity 2
    let A  = Z2
    let B  = I2
    let C  = matrixNegate I2
    let D  = matrixScale 2.0 I2
    createBlockMatrix [[A; B]; [C; D]]

let θ      = System.Math.PI * 20.0 / 180.0
let Lθ, Rθ = buildRotation4 θ
let F4     = buildForward4 ()

let x0 = Vector([|0.0; 0.0; 0.0; 1.0|])

let keys = [ F4; Rθ; Rθ; Rθ; Lθ; Lθ; F4 ]


let Mtotal =
  keys
  |> List.reduce (fun acc M -> BasicOps.MatrixProduct M acc)

let xSlut = BasicOps.MatVecProduct Mtotal x0


printfn "Slut-vektor:"
printfn "["
for i in 0..xSlut.Size-1 do
    printfn "  %.6f" xSlut.[i]
printfn "]"

printfn "C = (%.6f, %.6f), S = (%.6f, %.6f)" xSlut.[0] xSlut.[1] xSlut.[2] xSlut.[3]
