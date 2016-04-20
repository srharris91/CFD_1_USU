! user defined variables to define finite volume
! x and y direction # of cells
#define max_x  30
#define max_y  30
! x and y number of cells plus 1
#define max_xp  31
#define max_yp  31
! x and y number of cells plus 2 (to account for boundary nodes)
#define max_x2p  32
#define max_y2p  32



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
    REAL    ::  error_RSS,error2
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
    data(max_xp,:)%u= Ru
    ! top boundary
    data(:,max_yp)%u= Tu
    ! initialize u
    data%u_old  = data%u


    ! initialize v
    data%v_old = 0.
    data%v     = data%v_old

    ! initialize P values
    data%Pp=0.001
    data%P_old=0.00000001

    ! point SOR method to solve for the exact values of phi using the BC (only loop through inner values)
    ! solving using the deferred correction method
    CALL CPU_TIME(TIME1)
    CALL Solve_NS(data,dx,dy,max_x,max_y)
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
    open(unit=16,file="output/u_spot.txt")
    100 FORMAT (max_x2p ES16.7)
    101 FORMAT (2ES16.7)
    WRITE( 9,100) ( data(:,i)%xp,i=0,max_yp )
    WRITE(10,100) ( data(:,i)%yp,i=0,max_yp )
    WRITE(11,100) ( data(:,i)%xu,i=0,max_yp )
    WRITE(12,100) ( data(:,i)%yv,i=0,max_yp )
    WRITE(13,100) ( data(:,i)%u ,i=0,max_yp )
    WRITE(14,100) ( data(:,i)%v ,i=0,max_yp )
    WRITE(15,100) ( data(:,i)%P ,i=0,max_yp )
    DO i=0,max_xp
        IF (data(i,1)%xu <= 0.51 .AND. data(i,1)%xu>=0.49) THEN
            DO j=0,max_yp
                WRITE(16,101) data(i,j)%u,data(i,j)%yp
            END DO
        END IF
    END DO
    close(9);close(10);close(11);close(12);close(13);close(14);close(15);close(16)
END PROGRAM project3
