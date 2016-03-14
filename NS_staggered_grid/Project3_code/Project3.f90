! user defined variables to define finite volume
! x and y direction # of cells
#define max_x 40
#define max_y 40
! x and y number of cells plus 1
#define max_xp 41
#define max_yp 41
! x and y number of cells plus 2 (to account for boundary nodes)
#define max_x2p 42
#define max_y2p 42
MODULE types
    !purpose: define data type struct
    IMPLICIT NONE
    TYPE::dat
        REAL::x,y,phi,phi_new
        INTEGER::n
    END TYPE dat
CONTAINS
    SUBROUTINE set_xy (strct,dx,dy,nx,ny)
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx-1,0:nx-1),INTENT(INOUT)::strct ! data contained from 0:nx-1 where cells 0 and nx-1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        INTEGER ::  i,j,n   ! for do loops and n is counter for cell number
        REAL    ::  xi,yi   ! x and y values for each cell
        n=1                     ! cell number 1
        DO i=1,ny-2             ! 1 to ny-2 for boundary nodes (we only are iterating through the middle values)
            yi = i*dy - dy/2.   ! y coordinate
            DO j=1,nx-2
                xi = j*dx - dx/2.       ! x coordinate
                strct(i,j)%n = n        ! input n node
                strct(i,j)%x = xi       ! x coordinate to strct
                strct(i,j)%y = yi       ! y coordinate to strct
                strct(i,j)%phi = 0.     ! phi value initialized guess
                strct(i,j)%phi_new = 0. ! phi value initialized guess for next iteration
                n=n+1                   ! count cell numbers up one
            END DO
        END DO
    ! left boundary
    strct(:,0)%x = 0.
    strct(0,0)%y = 0.
    strct(1:max_x,0)%y = RESHAPE((/ (i*dy - dy/2. ,i=1,max_x) /),(/ max_x /))
    strct(max_xp,0)%y = 1.
    ! top boundary
    strct(0,0)%x = 0.
    strct(0,1:max_y)%x = RESHAPE((/ (i*dx - dx/2. ,i=1,max_y) /),(/ max_y /))
    strct(0,max_yp)%x = 1.
    strct(0,:)%y = 0.
    ! right boundary
    strct(:,max_yp)%x = 1.
    strct(:,max_yp)%y = strct(:,0)%y
    ! bottom boundary
    strct(max_xp,:)%x = strct(0,:)%x
    strct(max_xp,:)%y = 1.
    END SUBROUTINE set_xy
END MODULE types



PROGRAM project3
    USE types !use module defined by types
    IMPLICIT NONE
    ! declare variables
    INTEGER ::  i,j,iter,k!,max_x=20,max_y=20
    REAL    ::  dx,dy,gamma,ae,aw,an,as,ap,error
    REAL    ::  Lphi,Rphi,Tphi,Bphi ! boundary condition phi values
    REAL    ::  Fe,Fw,Fn,Fs ! flux terms = rho * u where rho = 1
    REAL    ::  aet,awt,ast,ant,apt ! higher order interpolation scheme coefficients
    REAL    ::  beta        ! used for combining lower order and higher order interpolation schemes
    REAL    ::  u=2.,v=2.   ! velocity values
    TYPE(dat),DIMENSION(0:max_xp,0:max_yp)::phi ! 22 if you count edges (thin cell)
    REAL    :: TIME1,TIME2  ! for time of computation
    REAL(KIND=4),DIMENSION(2):: TIMEA  ! for time of computation
    ! set dx and dy and gamma and coefficients (without dividing by delta x between node centers)
    dx=1./REAL(max_x)
    dy=1./REAL(max_y)
    gamma = 1.
    ! initialize phi data and x,y for middle values
    CALL set_xy(phi,dx,dy,max_x2p,max_x2p)
    !! initialize BC's
    ! BC's
    Lphi = 100.
    Rphi = 0.
    Tphi = 100.
    Bphi = 0.
    ! left Boundary
    phi(:,0)%phi = Lphi
    phi(:,0)%phi_new = Lphi
    ! bottom boundary
    phi(0,:)%phi = Bphi
    phi(0,:)%phi_new = Bphi
    ! right boundary
    phi(:,max_yp)%phi = Rphi
    phi(:,max_yp)%phi_new = Rphi
    ! top boundary
    phi(max_xp,:)%phi = Tphi
    phi(max_xp,:)%phi_new = Tphi
    ! point SOR method to solve for the exact values of phi using the BC (only loop through inner values)
    ! solving using the deferred correction method
    Fe = 1.*u
    Fw = 1.*u
    Fn = 1.*v
    Fs = 1.*v
    beta = 1.
    ae = 0.
    an = 0.
    aw = Fw
    as = Fs
    ap = as + aw
    aet= -Fe/2.
    awt= Fw/2.
    ant= -Fn/2.
    ast= Fs/2.
    apt= aet + awt + ant + ast
    open(unit=5,file="output/convergence.txt")
    CALL CPU_TIME(TIME1)
    ! step 1 solve discretised momentum equations

    ! step 2 Solve pressure correction equation

    ! step 3 Correct pressure and velocities

    ! step 4 Solve all other discretised transport equations

    ! if no convergence, then iterate

    ! if converged then stop


!    DO iter=0,100000000
!        error = 0.
!        DO j=1,max_x
!            DO i=1,max_y
!                    phi(i,j)%phi_new = (aw*phi(i,j-1)%phi_new + as*phi(i-1,j)%phi_new)/ap  &
!                        -beta/ap * (&
!                         apt*phi(i  ,j  )%phi_new &
!                        -aet*phi(i  ,j+1)%phi_new &
!                        -awt*phi(i  ,j-1)%phi_new &
!                        -ant*phi(i+1,j  )%phi_new &
!                        -ast*phi(i-1,j  )%phi_new &
!                        -ap *phi(i  ,j  )%phi_new &
!                        +aw *phi(i  ,j-1)%phi_new &
!                        +as *phi(i-1,j  )%phi_new )
!                IF (ISNAN(phi(i,j)%phi_new)) THEN
!                    WRITE(*,*) "You really stink at this dude!  error on ",i,j," value of phi_new.  On iter = ",iter
!                    WRITE(*,*) phi(i,j)%phi_new
!                    STOP
!                END IF
!                error = error + (phi(i,j)%phi-phi(i,j)%phi_new)**2  ! RSS the error for each iteration
!            END DO
!        END DO
!        ! BC ghost node values
!        phi(max_yp,:)%phi_new = phi(max_y,:)%phi_new
!        phi(:,max_xp)%phi_new = phi(:,max_x)%phi_new
!        error = SQRT(error)                             ! sqrt to have RSS of error
!        write(5,*) iter,error
!        phi%phi = phi%phi_new                           ! set old phi to the new iteration guess
!        IF (error .lt. 1.E-11) THEN
!            WRITE(*,*) 'Did in ',iter,' iterations'
!            WRITE(*,*) 'With RSS error = ',error
!            EXIT                                        ! exit loop if error is small enough
!        END IF
!    END DO
    
    CALL CPU_TIME(TIME2)
    WRITE(*,*) "CPU Time = ",TIME2-TIME1
    !output
    ! user will need to specify size of 
    open(unit=9,file="output/x.txt");open(unit=10,file="output/y.txt")
    open(unit=11,file="output/phi.txt");
    100 FORMAT (max_x2p F14.6)
    WRITE( 9,100) ( phi(i,:)%x ,i=0,max_xp )
    WRITE(10,100) ( phi(i,:)%y ,i=0,max_xp )
    WRITE(11,100) ( phi(i,:)%phi ,i=0,max_xp )
    k=max_xp
    close(9);close(10);close(11);close(5)
    WRITE(*,*) "Wall Time = ",etime(TIMEA)
END PROGRAM project3
