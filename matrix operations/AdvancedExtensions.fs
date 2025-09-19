module ProjectC

open System
open LinAlgDat.Core

type AdvancedOps = class

    static member SquareSubMatrix (A : Matrix) (i : int) (j : int) : Matrix =
        if i < 0 || i >= A.M_Rows || j < 0 || j >= A.N_Cols then
            raise (ArgumentException("Row or column index out of bounds"))
        
        let newRows = A.M_Rows - 1
        let newCols = A.N_Cols - 1
        let subMatrix = Matrix(newRows, newCols)

        for row in 0..newRows-1 do
            for col in 0..newCols-1  do
                let sourceRow = if row < i then row else row + 1
                let sourceCol = if col < j then col else col + 1
                subMatrix.[row, col] <- A.[sourceRow, sourceCol]

        subMatrix

    static member Determinant (A : Matrix) : float =
        if A.M_Rows <> A.N_Cols then
            raise (ArgumentException("Matrix must be square"))
        if A.M_Rows = 1 then
            A.[0, 0]
        elif A.M_Rows = 2 then
            A.[0, 0] * A.[1, 1] - A.[0, 1] * A.[1, 0]
        else
            let mutable det = 0.0
            for j in 0..A.N_Cols-1 do
                let subMatrix = AdvancedOps.SquareSubMatrix A 0 j
                let sign = if j % 2 = 0 then 1.0 else -1.0
                det <- det + sign * A.[0, j] * AdvancedOps.Determinant subMatrix
            det

    static member VectorNorm (v : Vector) =
        let mutable n2 = 0.0
        for i in 0..v.Size-1 do
            n2 <- n2 + v.[i] * v.[i]
        sqrt n2

    static member SetColumn (A : Matrix) (v : Vector) (j : int) : Matrix =
        // Tjek at kolonne‐indekset er validt
        if j < 0 || j >= A.N_Cols then
            raise (ArgumentException($"SetColumn: Kolonne‐indeks {j} er uden for rækkevidde (0..{A.N_Cols-1})"))

        // Tjek at vektorlængden svarer til antal rækker i A
        if v.Size <> A.M_Rows then
            raise (ArgumentException($"SetColumn: Vektor‐længde ({v.Size}) matcher ikke antal rækker i A ({A.M_Rows})"))

        // Kopiér alle elementer fra v ind i j'te kolonne i A
        for i in 0 .. A.M_Rows - 1 do
            A.[i, j] <- v.[i]
        A
    
    static member GramSchmidt (A : Matrix) : Matrix * Matrix =
        if A.M_Rows < A.N_Cols then
            raise (ArgumentException("Matrix must have at least as many rows as columns"))
        
        let m = A.M_Rows
        let n = A.N_Cols
        let Q = Matrix(m, n)
        let R = Matrix(n, n)

        for j in 0 .. n-1 do
            let mutable vj = Vector(A.Column j)
            for i in 0 .. j-1 do
                let qi = Vector(Q.Column i)
                R.[i, j] <- qi * vj
                vj <- vj - R.[i, j] * qi
            
            R.[j, j] <- AdvancedOps.VectorNorm vj
            if R.[j, j] < 1e-10 then
                raise (ArgumentException "Matrix columns are not linearly independent")

            // Normalisér og indsæt i Q
            AdvancedOps.SetColumn Q (vj * (1.0 / R.[j, j])) j |> ignore

        Q, R

end