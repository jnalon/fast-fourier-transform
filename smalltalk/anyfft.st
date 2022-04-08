"---------------------------------------------------------------------------------------------------
 Fast Fourier Transform -- Gnu Smalltalk version.
 This version implements Cooley-Tukey algorithm for composite numbers (not powers of 2 only).

 José Alexandre Nalon
 --------------------------------------------------------------------------------------------------"

"---------------------------------------------------------------------------------------------------
 Small library to deal with complex numbers. Based on Gnu Smalltalk tutorial.
 --------------------------------------------------------------------------------------------------"
Number subclass: Complex [
    | realPart imagPart |

    Complex class >> new [
        ^(super new) setReal: 0 setImag: 0
    ]
    Complex class >> real: r imag: i [
        ^(super new) setReal: r setImag: i
    ]

    setReal: r setImag: i [
        realPart := r asFloatE.
        imagPart := i asFloatE.
        ^self
    ]

    real [
        ^realPart
    ]
    imag [
        ^imagPart
    ]

    + val [
        ^Complex real: (realPart + val real)
                 imag: (imagPart + val imag)
    ]

    - val [
        ^Complex real: (realPart - val real)
                 imag: (imagPart - val imag)
    ]

    * val [
        ^Complex real: (realPart * val real) - (imagPart * val imag)
                 imag: (imagPart * val real) + (realPart * val imag)
    ]

    exp [
        | realExp |
        realExp := realPart exp.
        ^Complex real: realExp * (imagPart cos)
                 imag: realExp * (imagPart sin)
    ]

    printOn: aStream [
        realPart printOn: aStream.
        aStream nextPut: $+.
        imagPart printOn: aStream.
        aStream nextPut: $i
    ]
]


"---------------------------------------------------------------------------------------------------
 Extends SmallInteger to perform formatting and compute prime factors.
 --------------------------------------------------------------------------------------------------"
SmallInteger extend [

    format [
        | s |
        s := '      ', (self asString).
        ^(s copyFrom: ((s size) - 6) to: (s size))
    ]

    "Smallest prime factor of the number. If prime itself, return its own value."
    factor [
        | rn |
        rn := self sqrt floor.
        2 to: rn do: [ :i |
            (self rem: i) = 0
                ifTrue: [ ^i ]
        ].
        ^self.
    ]

]


"---------------------------------------------------------------------------------------------------
 Extend FloatE to format output.
 --------------------------------------------------------------------------------------------------"
FloatE extend [

    format [
        | intString fracString s |
        intString := self integerPart asInteger asString.
        fracString := (self fractionPart * 1000) integerPart asInteger asString.
        fracString := ((fracString asString), '0000') copyFrom: 1 to: 4.
        s := intString, '.', fracString.
        ((s size) <= 7)
            ifTrue: [
                s := '       ', s.
                ^(s copyFrom: ((s size) - 6) to: (s size))
            ]
            ifFalse: [
                ^(s copyFrom: 1 to: 7).
            ]
    ]

]


"---------------------------------------------------------------------------------------------------
 Extends Array to perform Fourier Transforms.
 --------------------------------------------------------------------------------------------------"
Array extend [

    "Initializes an array with a given value"
    Array class >> fill: length with: value [
        | instance |
        instance := super new: length.
        1 to: length do: [ :n | instance at: n put: value ].
        ^instance
    ]

    "Copy an Array on a given range, optionally skipping some values."
    copyFrom: start to: end by: step [
        | x |
        x := OrderedCollection new.
        start to: end by: step do: [ :i |
            x add: (self at: i)
        ].
        ^(x asArray)
    ]

    "Discrete Fourier Transform directly from the definition, an algorithm that has O(N^2)
     complexity."
    directFT [
        | N X W Wk Wkn |
        N := self size.                                            "Length of the array."
        X := Array fill: N with: (Complex new).                    "Accumulate the results."
        W := (Complex real: 0 imag: (-2 * FloatE pi / N)) exp.     "Twiddle factors."
        Wk := Complex real: 1 imag: 0.
        1 to: N do: [ :k |                                         "Compute the kth coefficient."
            Wkn := Complex real: 1 imag: 0.
            1 to: N do: [ :n |                                     "  Operate the summation."
                X at: k put: ((X at: k) + ((self at: n) * Wkn)).   "    Compute every term."
                Wkn := Wkn * Wk.                                   "Update twiddle factors."
            ].
            Wk := Wk * W.
        ].
        ^X
    ]

    "Fast Fourier Transform using a recursive decimation in time algorithm for vectors of length
    different from a power of 2. This has smaller complexity than the direct FT, though the exact
    value is difficult to compute."
    recursiveFFT [
        | N N1 N2 X Xj W Wj Wkj k2 |
        N := self size.
        N1 := N factor.
        N1 = N
            ifTrue: [ ^(self directFT) ]
            ifFalse: [
                N2 := N // N1.
                X := Array fill: N with: (Complex new).
                W := (Complex real: 0 imag: (-2 * FloatE pi / N)) exp.
                Wj := Complex real: 1 imag: 0.
                1 to: N1 do: [ :j |
                    Xj := (self copyFrom: j to: N by: N1) recursiveFFT.
                    Wkj := Complex real: 1 imag: 0.
                    1 to: N do: [ :k |
                        k2 := ((k - 1) rem: N2) + 1.
                        X at: k put: ((X at: k) + (Wkj * (Xj at: k2))).
                        Wkj := Wkj * Wj.
                    ].
                    Wj := Wj * W.
                ].
                ^X
            ]
    ]

]


"--------------------------------------------------------------------------------------------------
 Main Program
 --------------------------------------------------------------------------------------------------"
| x r n lengths dtime rtime itime |
'+---------+---------+---------+---------+---------+' displayNl.
'|    N    |   N^2   | N logN  | Direct  | Recurs. |' displayNl.
'+---------+---------+---------+---------+---------+' displayNl.

sizes := Array from: #(6 12 18 30 36 100 210 ).
5 to: 10 do: [ :r |
    n := 2 raisedToInteger: r.
    x := Array new: n.
    1 to: (x size) do: [ :j | x at: j put: (Complex real: j imag: 0) ].
    dtime := Time millisecondsToRun: [ 1 to: 50 do: [ :j | x directFT ] ].
    rtime := Time millisecondsToRun: [ 1 to: 50 do: [ :j | x recursiveFFT ] ].
    '| ' display. n format display.
    ' | ' display. (n * n) format display.
    ' | ' display. (n * r) format display.
    ' | ' display. (dtime / 50 / 1000 asFloatE) format display.
    ' | ' display. (rtime / 50 / 1000 asFloatE) format display.
    ' | ' displayNl.
].

'+---------+---------+---------+---------+---------+---------+' displayNl.
