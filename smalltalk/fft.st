"---------------------------------------------------------------------------------------------------
 Fast Fourier Transform -- Gnu Smalltalk version.
 This version implements Cooley-Tukey algorithm for powers of 2 only.

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
 Extends SmallInteger to perform formatting and bit reverse.
 --------------------------------------------------------------------------------------------------"
SmallInteger extend [

    format [
        | s |
        s := '      ', (self asString).
        ^(s copyFrom: ((s size) - 6) to: (s size))
    ]

    bitReverse: bits [
        | k l |
        k := self.
        l := 0.
        1 to: bits do: [ :i |
            l := (2 * l) + (k rem: 2).
            k := k // 2.
        ].
        ^l
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

    "Fast Fourier Transform using a recursive decimation in time algorithm. This has O(N log_2(N))
     complexity."
    recursiveFFT [
        | N N2 X Xe Xo Wk Xek WXok |
        N := self size.
        (N = 1)
            ifTrue: [ ^self ]                                      "A length-1 array is its own FT."
            ifFalse: [
                N2 := N // 2.                                      "Length of the array."
                X := Array fill: N with: (Complex new).
                Xe := (self copyFrom: 1 to: N by: 2) recursiveFFT. "Transform of the even samples."
                Xo := (self copyFrom: 2 to: N by: 2) recursiveFFT. "Transform of the odd samples."
                1 to: N2 do: [ :k |
                    Wk := (Complex real: 0 imag: (-2 * (k - 1) * FloatE pi / N)) exp.  "Twiddle factors."
                    Xek := Xe at: k.
                    WXok := Wk * (Xo at: k).
                    X at: k put: (Xek + WXok).
                    X at: (k + N2) put: (Xek - WXok).
                ].
                ^X
            ]
    ]

    "Performs bit reverse ordering of the Array."
    bitReverse: bits [
        | N X |
        N := self size.
        X := Array fill: N with: (Complex new).
        0 to: N - 1 do: [ :k |
            | l |
            l := k bitReverse: bits.
            X at: (l + 1) put: (self at: (k + 1)).
        ].
        ^X
    ]

    "Fast Fourier Transform using an interative in-place decimation in time algorithm. This has
     O(N log_2(N)) complexity, and since there are less message passing, it will probably be
     marginally faster than the recursive version."
    iterativeFFT [
        | N X W Wkn r step |
        N := self size.                                        "Length of the vector."
        r := (N log) // (2 log).                               "Number of bits."
        X := self bitReverse: r.                               "Reorder according to bit-reverse."
        step := 1.                                             "Auxiliary for computation of twiddle."
        1 to: r do: [ :k |
            0 to: N - 1 by: 2 * step do: [ :l |
                W := (Complex real: 0 imag: (-1 * FloatE pi / step)) exp.  "Twiddle factors."
                Wkn := Complex real: 1 imag: 0.
                0 to: (step - 1) do: [ :n |
                    | p q |
                    p := l + n + 1.
                    q := p + step.
                    X at: q put: ((X at: p) - (Wkn * (X at: q))).          "Recombine results."
                    X at: p put: ((X at: p) + (X at: p) - (X at: q)).
                    Wkn := Wkn * W.                                        "Update twiddle factors."
                ].
            ].
            step := 2 * step.
        ].
        ^X
    ]

]


"--------------------------------------------------------------------------------------------------
 Main Program
 --------------------------------------------------------------------------------------------------"
'+---------+---------+---------+---------+---------+---------+' displayNl.
'|    N    |   N^2   | N logN  | Direct  | Recurs. | Itera.  |' displayNl.
'+---------+---------+---------+---------+---------+---------+' displayNl.

repeats := 50.
5 to: 10 do: [ :r |
    n := 2 raisedToInteger: r.
    x := Array new: n.
    1 to: (x size) do: [ :j | x at: j put: (Complex real: j imag: 0) ].
    dtime := Time millisecondsToRun: [ 1 to: repeats do: [ :j | x directFT ] ].
    rtime := Time millisecondsToRun: [ 1 to: repeats do: [ :j | x recursiveFFT ] ].
    itime := Time millisecondsToRun: [ 1 to: repeats do: [ :j | x iterativeFFT ] ].
    '| ' display. n format display.
    ' | ' display. (n * n) format display.
    ' | ' display. (n * r) format display.
    ' | ' display. (dtime / repeats / 1000 asFloatE) format display.
    ' | ' display. (rtime / repeats / 1000 asFloatE) format display.
    ' | ' display. (itime / repeats / 1000 asFloatE) format display.
    ' | ' displayNl.
].

'+---------+---------+---------+---------+---------+---------+' displayNl.
