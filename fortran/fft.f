C ----------------------------------------------------------------------
C Fast Fourier Transform -- Fortran Version
C This version implements Cooley-Tukey algorithm for powers of 2 only.
C
C José Alexandre Nalon
C ----------------------------------------------------------------------
C Since Fortran is a compiled language, it must first be translated. It
C can be done with the following command:
C
C $ gfortran -o fft fft.f
C
C This will generate a `fft` executable file in your folder. To run,
C type the command:
C
C $ ./fft
C ----------------------------------------------------------------------
      Program FFT

          External DirectFT, IterativeFFT
          Complex X(0:1023), TX(0:1023)
          Real DTIME, ITIME, TimeIt
          Integer NX, I, R, REPEATS

          Write(*,*) "+---------+---------+---------" //
     c               "+---------+---------+"
          Write(*,*) "|    N    |   N^2   | N logN  " //
     c               "| Direta  | Itera.  |"
          Write(*,*) "+---------+---------+---------" //
     c               "+---------+---------+"

          REPEATS = 500
          Do 99, R = 5, 10

              NX = 2**R
              DTIME = TimeIt(DirectFT, NX, REPEATS)
              ITIME = TimeIt(IterativeFFT, NX, REPEATS)

              Write(*,900) NX, NX**2, NX*R, DTIME, ITIME
  900         Format(" | ", I7, " | ", I7, " | ", I7," | ",
     c                F7.4, " | ", F7.4, " |")

   99     Continue

          Write(*,*) "+---------+---------+---------" //
     c               "+---------+---------+"
          Stop

      End


C ----------------------------------------------------------------------
C Measure execution time through repeated calls to a (Fast) Fourier
C Transform function.
C
C Parameters:
C   F
C     Function to be called. It should receive three arguments, as
C     described in the Fourier Transform subroutines;
C   NX
C     Number of elements in the vector on which the transform will be
C     applied;
C   REPEATS
C     Number of times the function will be called.
C
C Returns:
C   The average execution time for that function with a vector of the
C   given size.
C ----------------------------------------------------------------------
      Real Function TimeIt(F, NX, REPEATS)
          External F
          Integer NX, REPEATS

          Complex X(0:1023), TX(0:1023)
          Real T0, T1
          Integer J

          Do 10, J = 0, NX-1
              X(J) = Cmplx(J)
   10     Continue

          Call CPU_Time(T0)
          Do 20, J = 1, REPEATS
              Call F(X, TX, NX)
   20     Continue
          Call CPU_Time(T1)
          TimeIt = (T1 - T0) / REPEATS
      end


C ----------------------------------------------------------------------
C Pretty printing of an array of complex numbers, used to inspect
C results.
C
C Parameters:
C   X
C     A vector of complex numbers.
C   NX
C     Number of elements on the vector.
C ----------------------------------------------------------------------
      SubRoutine ComplexShow(X, NX)
          Complex X(0:NX-1)
          Integer NX

          Do 510, I = 0, NX
              Write(*,520) X(I)
  520         Format(" ", 2F15.8)
  510     Continue
      End


C ----------------------------------------------------------------------
C Discrete Fourier Transform directly from the definition, an algorithm
C that has O(N^2) complexity.
C
C Parameters:
C   X
C     The vector of which the DFT will be computed. Given the nature of
C     the implementation, there is no restriction on the size of the
C     vector, although it will almost always be called with a power of
C     two size to give a fair comparison.
C   TX
C     The vector that will receive the result of the transformation.
C   NX
C     The number of elements in the vector.
C ----------------------------------------------------------------------
      SubRoutine DirectFT(X, TX, NX)
          Complex X(0:NX-1), TX(0:NX-1)
          Integer NX

          Complex W, WK, WKN
          Real PI
          Integer I, K
          Parameter (PI = 3.1415926)

          W = Cmplx(Cos(2*PI/NX), -Sin(2*PI/NX))
          WK = Cmplx(1, 0)
          Do 100, K = 0, NX-1
              TX(K) = 0
              WKN = Cmplx(1, 0)
              Do 110, I = 0, NX-1
                  TX(K) = TX(K) + WKN * X(I)
                  WKN = WKN * WK
  110         Continue
              WK = WK * W
  100     Continue
          Return
      End


C ----------------------------------------------------------------------
C Bit-reversed version of an integer number.
C
C Parameters:
C   K
C     The number to be bit-reversed;
C   R
C     The number of bits to take into consideration when reversing.
C
C Returns:
C   The number k, bit-reversed according to integers with r bits.
C ----------------------------------------------------------------------
      Integer Function BitReverse(K, R)
          Integer R, K

          Integer L, M, I

          L = 0
          M = K
          Do 200, I = 0, R-1
              L = L * 2 + Mod(M, 2)
              M = Int(M / 2)
  200     Continue
          BitReverse = L
          Return
      End


C ----------------------------------------------------------------------
C Fast Fourier Transform using a recursive decimation in time algorithm.
C This has O(N log_2(N)) complexity.
C
C Parameters:
C   X
C     The vector of which the FFT will be computed. This should always
C     be called with a vector of a power of two length, or it will fail.
C     No checks on this are made.
C   TX
C     The vector that will receive the results of the computation.
C   NX
C     The number of elements in the vector.
C ----------------------------------------------------------------------
      SubRoutine IterativeFFT(X, TX, NX)
          Complex X(0:NX-1), TX(0:NX-1)
          Integer NX

          Complex W, WKN
          Real PI
          Integer BitReverse
          Integer R, K, L, N, STEP, P, Q
          Parameter (PI = 3.1415926)

          R = Int(Log(Real(NX))/Log(2.0))

          Do 310, K = 0, NX-1
              L = BitReverse(K, R)
              TX(L) = X(K)
  310     Continue

          STEP = 1
          Do 320, K = 0, R-1
              Do 330, L = 0, NX-1, 2*STEP
                  W = Cmplx(Cos(PI/STEP), -Sin(PI/STEP))
                  WKN = Cmplx(1, 0)
                  Do 340, N = 0, STEP-1
                      P = L + N
                      Q = P + STEP
                      TX(Q) = TX(P) - WKN*TX(Q)
                      TX(P) = 2*TX(P) - TX(Q)
                      WKN = WKN * W
  340             Continue
  330         Continue
              STEP = STEP*2
  320     Continue

          Return
      End
