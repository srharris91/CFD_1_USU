! user defined variables to define finite volume
! x and y direction # of cells
#define max_x 20
#define max_y 20
! x and y number of cells plus 1
#define max_xp 21
#define max_yp 21
! x and y number of cells plus 2 (to account for boundary nodes)
#define max_x2p 22
#define max_y2p 22
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
    END SUBROUTINE set_xy
END MODULE types

PROGRAM project1
    USE types !use module defined by types
    IMPLICIT NONE
    ! declare variables
    INTEGER ::  i,j,iter!,max_x=20,max_y=20
    REAL    ::  dx,dy,gamma,deltaPhi,ae,aw,an,as,ap,error
    TYPE(dat),DIMENSION(0:max_xp,0:max_yp)::phi ! 22 if you count edges (thin cell)
    ! set dx and dy and gamma and coefficients (without dividing by delta x between node centers)
    dx=1./REAL(max_x)
    dy=1./REAL(max_y)
    gamma = 1.
    ! initialize phi data and x,y for middle values
    CALL set_xy(phi,dx,dy,max_x2p,max_x2p)
    !! initialize BC's
    ! left Boundary
    phi(:,0)%x = 0.
    phi(0,0)%y = 0.
    phi(1:max_x,0)%y = RESHAPE((/ (i*dy - dy/2. ,i=1,max_x) /),(/ max_x /))
    phi(max_xp,0)%y = 1.
    phi(:,0)%phi = 0.
    phi(:,0)%phi_new = 0.
    ! top boundary
    phi(0,0)%x = 0.
    phi(0,1:max_y)%x = RESHAPE((/ (i*dx - dx/2. ,i=1,max_y) /),(/ max_y /))
    phi(0,max_yp)%x = 1.
    phi(0,:)%y = 0.
    phi(0,:)%phi = 0.
    phi(0,:)%phi_new = 0.
    ! right boundary
    phi(:,max_yp)%x = 1.
    phi(:,max_yp)%y = phi(:,0)%y
    phi(:,max_yp)%phi = phi(:,max_yp)%y
    phi(:,max_yp)%phi_new = phi(:,max_yp)%y
    ! bottom boundary
    phi(max_xp,:)%x = phi(0,:)%x
    phi(max_xp,:)%y = 1.
    phi(max_xp,:)%phi = phi(max_xp,:)%x
    phi(max_xp,:)%phi_new = phi(max_xp,:)%x
    ! point SOR method to solve for the exact values of phi using the BC (only loop through inner 20X20 values)
    open(unit=5,file="convergence.txt")
    DO iter=0,100000
        error = 0.
        DO i=max_y,1,-1
            DO j=max_x,1,-1
                ae=dy*gamma/(abs(phi(i,j)%x-phi(i,j+1)%x)) 
                aw=dy*gamma/(abs(phi(i,j)%x-phi(i,j-1)%x)) 
                an=dx*gamma/(abs(phi(i,j)%y-phi(i+1,j)%y)) 
                as=dx*gamma/(abs(phi(i,j)%y-phi(i-1,j)%y)) 
                ap=(ae+aw+an+as) 
                deltaPhi = ae*phi(i,j+1)%phi_new + &    ! use new iteration values
                    aw*phi(i,j-1)%phi_new + &           ! use new iteration values
                    an*phi(i+1,j)%phi_new + &           ! use new iteration values
                    as*phi(i-1,j)%phi_new - &           ! use new iteration values
                    ap*phi(i,j)%phi                     ! use old iteration values
                phi(i,j)%phi_new = phi(i,j)%phi + 1.8/ap * deltaPhi     ! each iterations approximate new phi value
                error = error + deltaPhi**2             ! root sum square the error for each iteration
            END DO
        END DO
        error = SQRT(error)                             ! sqrt to have RSS of error
        write(5,*) iter,error
        phi%phi = phi%phi_new                           ! set old phi to the new iteration guess
        IF (error .lt. 1.E-14) THEN
            WRITE(*,*) 'Did in ',iter,' iterations'
            WRITE(*,*) 'With RSS error = ',error
            EXIT                                        ! exit loop if error is small enough
        END IF
    END DO
    !output
    ! user will need to specify size of 
    open(unit=9,file="x.txt");open(unit=10,file="y.txt");open(unit=11,file="phi.txt");open(unit=12,file="error.txt")
    100 FORMAT (max_x2p F10.6)
    200 FORMAT (max_x2p ES16.6)
    WRITE(9,100) ( phi(i,:)%x ,i=0,21 )
    WRITE(10,100) ( phi(i,:)%y ,i=0,21 )
    WRITE(11,100) ( phi(i,:)%phi ,i=0,21 )
    WRITE(12,200) ( abs(phi(i,:)%phi - (phi(i,:)%x*phi(i,:)%y)),i=0,21 )
    close(9);close(10);close(11);close(12);close(5);
END PROGRAM project1
