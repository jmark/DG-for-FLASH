# include "Flash.h"
# include "constants.h"

module dgfv_utils_mod

use Hydro_data, only: N_NODES !! number of DG nodes

integer, parameter :: dp = selected_real_kind(15)
integer, parameter :: RP = selected_real_kind(15)

!! Number of DG elements per block.
integer, parameter :: N_ELEMS = NXB/N_NODES

integer, parameter :: NXE = NXB/N_NODES
integer, parameter :: NXN = N_NODES

# if NDIM > 1
integer, parameter :: NYE = NYB/N_NODES
integer, parameter :: NYN = N_NODES
# else
integer, parameter :: NYE = 1
integer, parameter :: NYN = 1 
# endif

# if NDIM > 2
integer, parameter :: NZE = NZB/N_NODES
integer, parameter :: NZN = N_NODES
# else
integer, parameter :: NZE = 1
integer, parameter :: NZN = 1 
# endif

integer, parameter :: NXE_LO = 1-K1D
integer, parameter :: NXE_HI = 1+NXE*K1D

integer, parameter :: NYE_LO = 1-K2D
integer, parameter :: NYE_HI = 1+NYE*K2D

integer, parameter :: NZE_LO = 1-K3D
integer, parameter :: NZE_HI = 1+NZE*K3D

!! face codes
integer, parameter :: ZEN = 0
integer, parameter :: NOR = 1
integer, parameter :: SOU = 2
integer, parameter :: WES = 3
integer, parameter :: EAS = 4
integer, parameter :: FRO = 5
integer, parameter :: BAC = 6

interface pp
    procedure vp
    procedure mp
    procedure bp
    procedure tp
    procedure tp5
end interface

contains

!! --------------------------------------------------------------------- !!
!! trafo operators: DG <-> FV

subroutine fill_projection_matrix(N,nodes,weights,mat)

    IMPLICIT NONE

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)
    real(dp), intent(out)   :: mat(N,N)

    call subcellfv_vanderMondeMatrix(N,N,nodes,weights,mat)
    
end subroutine

subroutine fill_Weak_Diff_Matrix(N,nodes,weights,mat)

    IMPLICIT NONE

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)
    real(dp), intent(out)   :: mat(N,N)

    real(dp) :: temp(N,N)

    integer :: i,j

    temp = 0.0_dp
    call fill_Diff_Matrix(N-1,nodes,temp)

    do j = 1,N; do i = 1,N;
        mat(i,j) = 2.0_dp * weights(j)/weights(i) * temp(j,i)
    end do; end do;

    !! correct element length
    !mat = real(N_NODES,dp)/real(N,dp)*mat
    !mat = -mat/real(N,dp)

end subroutine

subroutine fill_surface_vectors(N,nodes,weights,surfVecM,surfVecP)

    IMPLICIT NONE

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)

    real(dp), intent(out)   :: surfVecM(N)
    real(dp), intent(out)   :: surfVecP(N)

    integer :: i

    do i = 1,N
        surfVecM(i) = -2.0_dp/weights(i)*lagrange_polynomial(N,nodes,i,-1.0_dp)
        surfVecP(i) =  2.0_dp/weights(i)*lagrange_polynomial(N,nodes,i, 1.0_dp)
    end do

end subroutine

pure function lagrange_polynomial(N, nodes, j, x) result(lp)
    
    integer, intent(in) :: N
    real(dp), intent(in)    :: nodes(N)
    integer, intent(in) :: j
    real(dp), intent(in)    :: x

    real(dp)    :: lp
    integer :: i

    real(dp), parameter :: TOL = 1.0e-14

    lp = 1.0_dp

    if (abs(x-nodes(j)) < TOL) return !! Kronecker property

    do i = 1,N
        if (i == j) cycle
        lp = lp * (x - nodes(i))/(nodes(j) - nodes(i))
    end do
     
end function

pure function lagrange_basis(N, nodes, x) result(lp)
    
    integer, intent(in) :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: x

    real(dp)    :: lp(N)
    integer :: i

    real(dp), parameter :: TOL = 1.0e-14

    !! Kronecker property
    if (any(abs(x-nodes) < TOL)) then
        lp = merge(1.0_dp,0.0_dp,abs(x-nodes) < 10.0_dp*TOL)
        return
    end if

    do i = 1,N
        lp(i) = Lagrange_Polynomial(N, nodes, i, x)
    end do
     
end function

pure subroutine lagrange_VanderMondeMatrix(N, M, xs, ys, Vdm)

    integer, intent(in)     :: N, m
    real(dp), intent(in)        :: xs(N), ys(M)
    real(dp), intent(out)       :: Vdm(M,N)

    integer :: i

    do i = 1,M
        Vdm(i,:) = Lagrange_Basis(N, xs, ys(i))
    end do

end subroutine

subroutine subcellfv_midpoints(N, sfvw, sfvx, sfvb)

    integer,intent(in)      :: N             !! number of sub-cells
    real(dp),intent(out)    :: sfvw          !! width of sub-volumes in reference space
    real(dp),intent(out)    :: sfvx(N)       !! cell-centers of sub-cells in reference space
    real(dp),intent(out)    :: sfvb(N+1)     !! positions of boundaries of sub-cells in reference space 

    integer :: i

    sfvw  = 2.0_dp/N

    sfvb(1) = -1.0_dp

    do i = 1,N
      sfvb(i+1) =         sfvb(i) + sfvw
      sfvx(i)   = 0.5_dp*(sfvb(i) + sfvb(i+1))
    end do

end subroutine

subroutine subcellfv_vanderMondeMatrix(N, M, nodes, weights, Vdm)

    integer, intent(in)     :: N, M
    real(dp), intent(in)        :: nodes(N)
    real(dp), intent(in)        :: weights(N)
    real(dp), intent(out)       :: Vdm(M,N)

    real(dp) :: sfvw
    real(dp) :: sfvx(M)
    real(dp) :: sfvb(M+1)

    real(dp) :: sfvVdm(N,N)
    real(dp) :: sfvnodes(N)

    integer :: i,j,k

    Vdm = 0.0_dp

    CALL subcellfv_midpoints(M, sfvw, sfvx, sfvb)

    do i = 1,M

        sfvnodes = sfvb(i) + 0.5_dp*(nodes+1.0_dp) * (sfvb(i+1) - sfvb(i))

        call lagrange_VanderMondeMatrix(N, N, nodes, sfvnodes, sfvVdm)

        do j = 1,N
            do k = 1,N
                Vdm(i,j) = Vdm(i,j) + 0.5_dp * weights(k) * sfvVdm(k,j)
            end do
        end do

    end do

end subroutine

subroutine fill_generalized_SBP_surface_Matrix(N, nodes, weights, surfMat)

    integer, intent(in)     :: N
    real(dp), intent(in)    :: nodes(N)
    real(dp), intent(in)    :: weights(N)

    real(dp), intent(out)   :: surfMat(N_NODES,N_NODES)

    real(dp) :: matM(N_NODES,N_NODES)
    real(dp) :: matP(N_NODES,N_NODES)

    real(dp) :: Bzero(2,2)
    real(dp) :: Vface(2,N_NODES)

    real(dp) :: lagrM(N_NODES)
    real(dp) :: lagrP(N_NODES)

    integer :: i,j

    Bzero = 0.0_dp
    Bzero(1,1) = -1.0_dp
    Bzero(2,2) =  1.0_dp

    do i = 1,N_NODES
        ! lagrM(i) = 2.0_dp/REAL(N,dp)/weights(i)*lagrange_polynomial(N, nodes, i, -1.0_dp)
        ! lagrP(i) = 2.0_dp/REAL(N,dp)/weights(i)*lagrange_polynomial(N, nodes, i,  1.0_dp)
        lagrM(i) = lagrange_polynomial(N, nodes, i, -1.0_dp)
        lagrP(i) = lagrange_polynomial(N, nodes, i,  1.0_dp)
    end do

    Vface(1,:) = lagrM
    Vface(2,:) = lagrP

    surfMat = matmul(transpose(Vface),matmul(Bzero,Vface))

    do j = 1,N_NODES
    do i = 1,N_NODES
        surfMat(i,j) = surfMat(i,j)/weights(i)
    end do
    end do

    ! matM = 0.0_dp
    ! matP = 0.0_dp

    ! matM(:,1      ) = lagrM
    ! matP(:,N_NODES) = lagrP

    ! surfMat = matM - matP
    
end subroutine

pure function invert_matrix(n,o) result(c)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    !===========================================================

    integer, intent(in) :: n
    real(dp), intent(in) :: o(n,n)
    real(dp)             :: a(n,n), c(n,n)
    real(dp)             :: L(n,n), U(n,n), b(n), d(n), x(n)
    real(dp)             :: coeff
    integer             :: i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    a = o
    L = 0.0_dp
    U = 0.0_dp
    b = 0.0_dp

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0_dp
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0_dp
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0_dp
    end do

end function

pure ELEMENTAL SUBROUTINE LegendrePolynomialAndDerivative(N_in,x,L,Lder)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in   !< (IN)  polynomial degree, (N+1) CLpoints
real(dp),INTENT(IN)    :: x      !< (IN)  coordinate value in the interval [-1,1]
real(dp),INTENT(OUT)   :: L      !< (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
real(dp),INTENT(OUT)   :: Lder   !< (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iLegendre
real(dp)    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
real(dp)    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
IF(N_in .EQ. 0)THEN
  L=1.0_dp
  Lder=0.0_dp
ELSEIF(N_in .EQ. 1) THEN
  L=x
  Lder=1.0_dp
ELSE ! N_in > 1
  L_Nm2=1.0_dp
  L_Nm1=x
  Lder_Nm2=0.0_dp
  Lder_Nm1=1.0_dp
  DO iLegendre=2,N_in
    L=(real(2*iLegendre-1,dp)*x*L_Nm1 - real(iLegendre-1,dp)*L_Nm2)/real(iLegendre,dp)
    Lder=Lder_Nm2 + real(2*iLegendre-1,dp)*L_Nm1
    L_Nm2=L_Nm1
    L_Nm1=L
    Lder_Nm2=Lder_Nm1
    Lder_Nm1=Lder
  END DO !iLegendre=2,N_in
END IF ! N_in
!normalize
L=L*SQRT(real(N_in,dp)+0.5_dp)
Lder=Lder*SQRT(real(N_in,dp)+0.5_dp)
END SUBROUTINE LegendrePolynomialAndDerivative

pure SUBROUTINE LegendreGaussNodesAndWeights(N_in,xGP,wGP)

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: N_in              !< polynomial degree, (N_in+1) Gausspoints
real(dp),INTENT(OUT)          :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
real(dp),INTENT(OUT),OPTIONAL :: wGP(0:N_in)       !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, parameter          :: nIter = 10        ! max. number of newton iterations
real(dp), parameter             :: Tol   = 1.E-15_dp    ! tolerance for Newton iteration: TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
real(dp)                      :: L_Np1,Lder_Np1    ! L_{N_in+1},Lder_{N_in+1}
real(dp)                      :: dx                ! Newton step
real(dp)                      :: cheb_tmp          ! temporary variable for evaluation of chebychev node positions

!==================================================================================================================================
IF(N_in .EQ. 0) THEN
  xGP=0.0_dp
  IF(PRESENT(wGP))wGP=2.0_dp
  RETURN
ELSEIF(N_in.EQ.1)THEN
  xGP(0)=-sqrt(1.0_dp/3.0_dp)
  xGP(N_in)=-xGP(0)
  IF(PRESENT(wGP))wGP=1.0_dp
  RETURN
ELSE ! N_in>1
  cheb_tmp=2.0_dp*atan(1.0_dp)/real(N_in+1,dp) ! pi/(2N+2)
  DO iGP=0,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cheb_tmp*real(2*iGP+1,dp)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
      dx=-L_Np1/Lder_Np1
      xGP(iGP)=xGP(iGP)+dx
      IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      !WRITE(*,*) 'maximum iteration steps >10 in Newton iteration for Legendre Gausspoint'
      xGP(iGP)=-cos(cheb_tmp*real(2*iGP+1,dp)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        !WRITE(*,*) iter,xGP(iGP)    !DEBUG
        CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
        dx=-L_Np1/Lder_Np1
        xGP(iGP)=xGP(iGP)+dx
        IF(abs(dx).LT.Tol*abs(xGP(iGP))) EXIT
      END DO !iter
      !CALL abort(__STAMP__,&
      !           'ERROR: Legendre Gauss nodes could not be computed up to desired precision. Code stopped!')
    END IF ! (iter.GT.nIter)
    CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      !wGP(iGP)=2./((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1) !if Legendre not normalized
      wGP(iGP)=(2.0_dp*N_in+3)/((1.0_dp-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO !iGP
END IF ! N_in
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.
  CALL LegendrePolynomialAndDerivative(N_in+1,xGP(N_in/2),L_Np1,Lder_Np1)
  !IF(PRESENT(wGP))wGP(N_in/2)=2./(Lder_Np1*Lder_Np1) !if Legendre not normalized
  IF(PRESENT(wGP))wGP(N_in/2)=(2.0_dp*N_in+3.0_dp)/(Lder_Np1*Lder_Np1)
END IF ! (mod(N_in,2) .EQ. 0)

END SUBROUTINE

pure SUBROUTINE qAndLEvaluation(N,x,q,dq,L_N)

    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: N
    real(dp), INTENT(IN)   :: x
    real(dp), INTENT(OUT)  :: L_N,q,dq
!
    ! lokale Variablen
    INTEGER                     :: k
    real(dp)               :: L_N_1,L_N_2,L_k,dL_N,dL_N_1,dL_N_2,dL_k
!
    ! NUR FUER N>=2 !
!
    ! Initialisieren
    L_N_2 = 1.0_RP
    L_N_1 = x
    dL_N_2 = 0.0_RP
    dL_N_1 = 1.0_RP
!
    DO k=2,N
        ! Rekursionsformel verwenden
        L_N = (2.0_RP*real(k,dP)-1.0_RP)/real(k,dP) * x * L_N_1 - (real(k,dP)-1.0_RP)/real(k,dP) * L_N_2
        ! Rekursionsformel fuer die Ableitungen verwenden
        dL_N = dL_N_2 + (2.0_RP*real(k,dP)-1.0_RP) * L_N_1
        ! Werte fuer den naechsten Schritt updaten
        L_N_2 = L_N_1
        L_N_1 = L_N
        dL_N_2 = dL_N_1
        dL_N_1 = dL_N
    END DO
!
    ! Einen weiteren Schritt gehen
    k = N+1
    ! Rekursionsformel verwenden
    L_k = (2.0_RP*real(k,dP)-1.0_RP)/real(k,dP) * x * L_N - (real(k,dP)-1.0_RP)/real(k,dP) * L_N_2
    ! Rekursionsformel fuer die Ableitungen verwenden
    dL_k = dL_N_2 + (2.0_RP*real(k,dP)-1.0_RP) * L_N_1
    ! Benoetigtes Polynom und Ableitung bestimmen
    q = L_k - L_N_2
    dq = dL_k - dL_N_2
!
    RETURN
END SUBROUTINE qAndLEvaluation

! ------------------------------------------------------------------------------------------------------------- !

pure SUBROUTINE LegendreGaussLobattoNodesAndWeights(N,x,w)

    IMPLICIT NONE
    INTEGER, INTENT(IN)                         :: N
    real(dp), DIMENSION(0:N), INTENT(OUT)  :: x,w
!
    ! lokale Variablen
    INTEGER, PARAMETER            :: n_it = 4
    INTEGER                       :: j,k
    real(dp), PARAMETER      :: tol = 4.0E-16
    real(dp)                 :: q,dq,L_N,Delta

    L_N = 0.0
!
    ! Fuer den einfachsten Fall Knoten und Gewicht explizit angeben
    IF (N.EQ.1) THEN
        x(0) = -1.0_RP
        w(0) = 1.0_RP
        x(1) = 1.0_RP
        w(1) = w(0)
!
    ELSE
        ! Randwerte explizit vorgeben
        x(0) = -1.0_RP
        w(0) = 2.0_RP / ( real(N,dP)*(real(N,dP)+1) )
        x(N) = 1.0_RP
        w(N) = w(0)
        ! Ansonsten: Symmetrie ausnutzen (nur fuer die erste Haelfte berechnen)
        DO j=1,((N+1)/2 - 1)
            ! Anfangsschaetzung mittels asymptotischer Relation nach Parter:
            ! Setzte negatives Vorzeichen, um die Stuetzstellen aus dem Intervall (-1,0) zu erhalten.
            x(j) = -cos( ((real(j,dP)+0.25_RP)*PI)/real(N,dP) - (3.0_RP/(8_RP*real(N,dP)*PI)) * (1/(real(j,dP)+0.25_RP)) )
!
            ! Diese Schaetzung wird nun mithilfe des Newton-Verfahrens praezisiert
            DO k=0,n_it
                ! Dazu wird das (N+1)-ste Polynom sowie dessen Ableitung benoetigt
                ! (ausgewertet an der entspr. Stuetzstelle)
                CALL qAndLEvaluation(N,x(j),q,dq,L_N)
                ! Newton-Korrektur berechnen und anwenden
                Delta = -q/dq
                x(j) = x(j) + Delta
                ! Falls vor Abarbeiten aller festgelegten Iterationen des Newton-Verfahrens der berechnete Wert
                ! schon nahe genug an der Nullstelle liegt, wird abgebrochen.
                IF ( ABS(Delta).LE.(tol * ABS(x(j))) ) THEN
                    EXIT
                END IF
!
            END DO
!
            CALL qAndLEvaluation(N,x(j),q,dq,L_N)
            ! Nutze Symmetrie aus (Multiplikation mit -1)
            x(N-j) = -x(j)
            ! Gewichte berechnen
            w(j) = 2.0_RP / ( real(N,dP)*(real(N,dP)+1.0_RP) * (L_N*L_N) )
            ! Auch die Gewichte sind symmetrisch
            w(N-j) = w(j)
!
        END DO
!
    END IF
!
    ! Falls eine ungerade Anzahl gefordert wird, noch die "mittleren" Werte angeben
    IF (MOD(N,2).EQ.0) THEN
        CALL qAndLEvaluation(N,0.0_RP,q,dq,L_N)
        ! Berechnung der Stuetzstellen und Gewichte analog
        x(N/2) = 0.0_RP
        w(N/2) = 2.0_RP / ( real(N,dP) * (real(N,dP)+1.0_RP) * (L_N*L_N) )
    END IF
!
    RETURN
END SUBROUTINE LegendreGaussLobattoNodesAndWeights

pure SUBROUTINE BarycentricWeights(N_in,xGP,wBary)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in               !< polynomial degree
real(dp),INTENT(IN)    :: xGP(0:N_in)        !< Gauss point positions for the reference interval [-1,1]
real(dp),INTENT(OUT)   :: wBary(0:N_in)      !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,jGP
!==================================================================================================================================
wBary(:)=1.
DO iGP=1,N_in
  DO jGP=0,iGP-1
    wBary(jGP)=wBary(jGP)*(xGP(jGP)-xGP(iGP))
    wBary(iGP)=wBary(iGP)*(xGP(iGP)-xGP(jGP))
  END DO ! jGP
END DO ! iGP
wBary(:)=1./wBary(:)
END SUBROUTINE BarycentricWeights

pure SUBROUTINE fill_Diff_Matrix(N_in,xGP,D)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in              !< polynomial degree
real(dp),INTENT(IN)    :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
real(dp),INTENT(OUT)   :: D(0:N_in,0:N_in)  !< differentiation Matrix
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,iLagrange
real(dp)               :: wBary(0:N_in)
!==================================================================================================================================
CALL BarycentricWeights(N_in,xGP,wBary)
D(:,:)=0.
DO iLagrange=0,N_in
  DO iGP=0,N_in
    IF(iLagrange.NE.iGP)THEN
      D(iGP,iLagrange)=wBary(iLagrange)/(wBary(iGP)*(xGP(iGP)-xGP(iLagrange)))
      D(iGP,iGP)=D(iGP,iGP)-D(iGP,iLagrange)
    END IF ! (iLagrange.NE.iGP)
  END DO ! iGP
END DO ! iLagrange
END SUBROUTINE

subroutine tp5(title,blk)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: blk(:,:,:,:,:)

    integer :: i,k,l,q

    do q = 1,size(blk,5)
        do l = 1,size(blk,4)
            do k = 1,size(blk,3)
                write (*,'((a12),(a6),(i3),(a6),(i3),(a6),(i3))') trim(title),' | q = ', q, ' | w = ', l, ' | z = ', k
                do i = 1,size(blk,1)
                    write (*,'(99(ES12.4))') blk(i,:,k,l,q)
                end do
                write (*,*)
            end do
            write (*,*)
        end do
        write (*,*)
    end do
    write (*,*)

end subroutine

subroutine tp(title,blk)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: blk(:,:,:,:)

    integer :: i,k,l

    do l = 1,size(blk,4)
        do k = 1,size(blk,3)
            write (*,'((a12),(a6),(i3),(a6),(i3))') trim(title),' | w = ', l, ' | z = ', k
            do i = 1,size(blk,1)
                write (*,'(99(ES12.4))') blk(i,:,k,l)
            end do
            write (*,*)
        end do
        write (*,*)
    end do
    write (*,*)

end subroutine

subroutine bp(title,blk)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: blk(:,:,:)

    integer :: i,k

    do k = 1,size(blk,3)
        write (*,*) trim(title) // ' | z = ' // tostr(k)
        do i = 1,size(blk,1)
            write (*,'(99(ES12.4))') blk(i,:,k)
        end do
        write (*,*)
    end do
    write (*,*)

end subroutine

subroutine mp(title,mat)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: mat(:,:)

    integer :: i

    write (*,*) trim(title)
    do i = 1,size(mat,1)
        write (*,'(99(ES12.4))') mat(i,:)
    end do
    write (*,*)

end subroutine

subroutine vp(title,vec)

    character(len=*), intent(in)    :: title
    real(dp), intent(in)            :: vec(:)

    integer :: i

    write (*,*) trim(title)
    write (*,'(99(ES12.4))') vec
    write (*,*)

end subroutine

function tostr(num) result(str)

    integer, intent(in) :: num
    character(len=16)   :: str

    write(str,'(i0)') num

end function

end module
