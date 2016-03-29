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
    ! left Boundary
    data(0,:)%u     = Lu
    ! bottom boundary
    data(:,0)%u     = Bu
    ! right boundary
    data(max_yp,:)%u= Ru
    ! top boundary
    data(:,max_xp)%u= Tu
    ! initialize u
    data%u_old = data%u


    ! initialize v
    data%v_old = 0.
    data%v     = 0.
    ! initialize P values
    data%Pe=1.
    data%Pn=1.
    data%Pw=1.
    data%Ps=1.
    ! point SOR method to solve for the exact values of phi using the BC (only loop through inner values)
    ! solving using the deferred correction method
    CALL CPU_TIME(TIME1)
    DO iter=0,10000
    ! step 1 solve discretised momentum equations
    CALL mom_uv(data,dx,dy,max_x,max_x)

    ! step 2 Solve pressure correction equation

    ! step 3 Correct pressure and velocities

    ! step 4 Solve all other discretised transport equations

    ! if no convergence, then iterate
    error_RSS = 0.
    DO i=1,max_x
        DO j=1,max_y
            error_RSS = error_RSS + (data(i,j)%u-data(i,j)%u_old)**2
            error_RSS = error_RSS + (data(i,j)%v-data(i,j)%v_old)**2
        END DO
    END DO
    error_RSS = sqrt(error_RSS)
    ! reset values
    data%u_old = data%u
    data%v_old = data%v
    !WRITE(*,*) "error = ",error_RSS


    ! if converged then stop
    IF (error_RSS < 0.000001) THEN
        EXIT
    ELSE
        WRITE(*,*) "iteration and error = ",iter,error_RSS
    END IF
END DO


    
CALL CPU_TIME(TIME2)
WRITE(*,*) "CPU Time = ",TIME2-TIME1
!output
! user will need to specify size of 
open(unit=9,file="output/x.txt")
open(unit=10,file="output/y.txt")
open(unit=11,file="output/u.txt")
100 FORMAT (max_x2p F14.6)
WRITE( 9,100) ( data(:,i)%x ,i=0,max_yp )
WRITE(10,100) ( data(:,i)%y ,i=0,max_yp )
WRITE(11,100) ( data(:,i)%u ,i=0,max_yp )
close(9);close(10);close(11)
END PROGRAM project3
