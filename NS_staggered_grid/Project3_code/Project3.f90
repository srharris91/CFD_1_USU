! user defined variables to define finite volume
! x and y direction # of cells
#define max_x 10
#define max_y 10
! x and y number of cells plus 1
#define max_xp 11
#define max_yp 11
! x and y number of cells plus 2 (to account for boundary nodes)
#define max_x2p 12
#define max_y2p 12



PROGRAM project3
    USE types !use module defined by types
    IMPLICIT NONE
    ! declare variables
    INTEGER ::  i,j,iter!,max_x=20,max_y=20
    REAL    ::  dx,dy
    REAL    ::  Lu,Ru,Tu,Bu ! boundary condition u velocity values
    !REAL    ::  u=0.,v=0.   ! velocity values
    TYPE(dat),DIMENSION(0:max_xp,0:max_yp)::data ! 22 if you count edges (thin cell)
    REAL    ::  TIME1,TIME2  ! for time of computation
    !REAL(KIND=4),DIMENSION(2):: TIMEA  ! for time of computation
    REAL    ::  error_RSS
    ! set dx and dy and gamma and coefficients (without dividing by delta x between node centers)
    dx=1./REAL(max_x)
    dy=1./REAL(max_y)
    ! initialize data and x,y for middle values
    CALL set_xy(data,dx,dy,max_xp,max_xp)
    !! initialize BC's
    ! BC's
    Lu = 0.
    Ru = 0.
    Tu = 1.
    Bu = 0.
    data%u          = 0.  ! initialize all data
    ! left Boundary
    data(0,:)%u     = Lu
    ! bottom boundary
    data(:,0)%u     = Bu
    ! right boundary
    data(max_yp,:)%u= Ru
    ! top boundary
    data(:,max_xp)%u= Tu
    ! initialize u
    data%u_old  = data%u


    ! initialize v
    data%v_old = 0.
    !data(4,2:max_x-2)%v_old = -1.24! initialize strange value to get v to converge
    !data(38,2:max_x-2)%v_old = 1.24! initialize strange value to get v to converge
    data%v     = data%v_old
    
    ! initialize P values
    data%Pp=0.
    data%P_old=0.
    !data(4,2:max_x-2)%P_old = -5.24! initialize strange value to get p to converge
    !data(38,2:max_x-2)%P_old = 5.24! initialize strange value to get p to converge

    ! initialize velocity correction terms

    ! point SOR method to solve for the exact values of phi using the BC (only loop through inner values)
    ! solving using the deferred correction method
    CALL CPU_TIME(TIME1)
    write(*,"(12ES16.7)") (data(0:max_xp,j)%u,j=max_yp,0,-1)
    !write(*,"(12ES16.7)") (data(0:max_xp,j)%v,j=max_yp,0,-1)
    DO iter=0,20000
        ! step 1 solve discretised momentum equations
        CALL mom_uv(data,dx,dy,max_x,max_y)

        ! step 2 Solve pressure correction equation
        ! step 3 Correct pressure and velocities
        !CALL vel_correction(data,dx,dy,max_x,max_y)


        ! step 4 Solve all other discretised transport equations
        ! do we need this in this problem?

        ! if no convergence, then iterate
        error_RSS = 0.
        DO i=1,max_x
            DO j=1,max_y
                error_RSS = error_RSS + (data(i,j)%u-data(i,j)%u_old)**2
                error_RSS = error_RSS + (data(i,j)%v-data(i,j)%v_old)**2
                error_RSS = error_RSS + (data(i,j)%Pp-data(i,j)%P_old)**2
            END DO
        END DO
        error_RSS = sqrt(error_RSS)
        ! reset values
        !data%u_old = data%u
        data%v_old = data%v
        data%P_old = data%Pp
        !WRITE(*,*) "error = ",error_RSS


        ! if converged then stop
        IF (error_RSS < 0.000000000001) THEN
            WRITE(*,*) "converged on iteration and error big loop = ",iter,error_RSS
            EXIT
        ELSE
            WRITE(*,*) "iteration and error big loop = ",iter,error_RSS
        END IF
    END DO

    write(*,"(12ES16.7)") (data(0:max_xp,j)%u,j=max_yp,0,-1)


CALL CPU_TIME(TIME2)
WRITE(*,*) "CPU Time = ",TIME2-TIME1
!output
! user will need to specify size of 
open(unit= 9,file="output/x.txt")
open(unit=10,file="output/y.txt")
open(unit=11,file="output/xu.txt")
open(unit=12,file="output/yv.txt")
open(unit=13,file="output/u.txt")
open(unit=14,file="output/v.txt")
open(unit=15,file="output/P.txt")
!100 FORMAT (max_x2p F14.6)
100 FORMAT (max_x2p ES16.7)
WRITE( 9,100) ( data(:,i)%xp,i=0,max_yp )
WRITE(10,100) ( data(:,i)%yp,i=0,max_yp )
WRITE(11,100) ( data(:,i)%xu,i=0,max_yp )
WRITE(12,100) ( data(:,i)%yv,i=0,max_yp )
WRITE(13,100) ( data(:,i)%u ,i=0,max_yp )
WRITE(14,100) ( data(:,i)%v ,i=0,max_yp )
WRITE(15,100) ( data(:,i)%Pp,i=0,max_yp )
close(9);close(10);close(11);close(12);close(13);close(14);close(15)
END PROGRAM project3
