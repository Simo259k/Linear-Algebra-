module ProjectB
// LinalgDat23
// Authors: Francois Lauze

open System
open LinAlgDat.Core

type GaussOps = class
    static member AugmentRight (A : Matrix) (v : Vector) : Matrix =
        let m_rows = A.M_Rows
        let n_cols = A.N_Cols

        let retval = Array2D.zeroCreate m_rows (n_cols + 1)

        for i in 0..m_rows-1 do
            for j in 0..n_cols-1 do
                retval.[i,j] <- A.[i,j]
            retval.[i,n_cols] <- v.[i]
        Matrix retval

    static member ElementaryRowReplacement (A : Matrix) (i : int) (m : float) (j : int) : Matrix =
        let n_rows = A.M_Rows
        let n_cols = A.N_Cols
        let retval = Matrix(A)
        for k in 0..n_cols-1 do
            retval.[i,k] <- retval.[i,k] + m * retval.[j,k]
        retval

    static member ElementaryRowInterchange (A : Matrix) (i : int) (j : int) : Matrix =
        let n_rows = A.M_Rows
        let n_cols = A.N_Cols
        let retval = Matrix(A)
        for k in 0..n_cols-1 do
            let temp = retval.[i,k]
            retval.[i,k] <- retval.[j,k]
            retval.[j,k] <- temp
        retval

    static member ElementaryRowScaling (A : Matrix) (i : int) (c : float) : Matrix =
        let n_rows = A.M_Rows
        let n_cols = A.N_Cols
        let retval = Matrix(A)
        for k in 0..n_cols-1 do
            retval.[i,k] <- c * retval.[i,k]
        retval

    static member ForwardReduction (M : Matrix) : Matrix =
        let n_rows = M.M_Rows
        let n_cols = M.N_Cols
        let mutable M = Matrix(M)
        let tol = 0.00000001
        for i in 0 .. min (n_rows-2) (n_cols-1) do
            let mutable j = i
            while j < n_cols && abs M.[i,j] < tol do
                j <- j + 1

            if j = n_cols then
                ()
            else
                match [i..n_rows-1] |> List.tryFind (fun r -> abs M.[r,j] > tol) with
                | None -> ()
                | Some k ->
                    if k <> i then M <- GaussOps.ElementaryRowInterchange M i k
                    for r in i+1 .. n_rows-1 do
                        let pivot  = M.[i,j]
                        let factor = - M.[r,j] / pivot
                        M <- GaussOps.ElementaryRowReplacement M r factor i
        M

    static member BackwardReduction (M0 : Matrix) : Matrix =
        let m = M0.M_Rows
        let n = M0.N_Cols

        let mutable M = Matrix M0
        let tol = 1e-7


        let pivots =
            [ for i in 0 .. m-1 do
                match [0 .. n-1] |> List.tryFind (fun j -> abs M.[i,j] > tol) with
                | Some j -> yield (i,j)
                | None   -> () ]

        pivots
        |> List.rev
        |> List.iter (fun (r,c) ->
            let piv = M.[r,c]
            if abs piv > tol then
                M <- GaussOps.ElementaryRowScaling     M r (1.0 / piv)

                for i in 0 .. r-1 do
                    let factor = M.[i,c]
                    if abs factor > tol then
                        M <- GaussOps.ElementaryRowReplacement M i (-factor) r
        )
        M

    static member GaussElimination (A : Matrix) (b : Vector) : Vector =
        let n = A.M_Rows
        let mutable M = Matrix A
        let mutable v = Vector b
        let tol = 1e-7

        for i in 0 .. n-2 do
            let pivotRow =
                [i .. n-1]
                |> List.maxBy (fun k -> abs M.[k,i])
            if pivotRow <> i then
                M <- GaussOps.ElementaryRowInterchange M i pivotRow
                let tmp = v.[i]
                v.[i] <- v.[pivotRow]
                v.[pivotRow] <- tmp

            let piv = M.[i,i]
            if abs piv < tol then
                invalidOp $"Singulær matrix (pivot ≈ 0 på række {i})"

            for j in i+1 .. n-1 do
                let factor = M.[j,i] / piv
                // træk factor * række i fra række j
                for k in i .. n-1 do
                    M.[j,k] <- M.[j,k] - factor * M.[i,k]
                v.[j] <- v.[j] - factor * v.[i]

        let x = Vector n
        for i in n-1 .. -1 .. 0 do
            let mutable sum = v.[i]
            for k in i+1 .. n-1 do
                sum <- sum - M.[i,k] * x.[k]
            let piv = M.[i,i]
            if abs piv < tol then
                invalidOp $"Singulær matrix (pivot ≈ 0 på række {i})"
            x.[i] <- sum / piv

        x


    end
