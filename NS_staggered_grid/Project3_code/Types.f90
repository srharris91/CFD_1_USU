MODULE types
    !purpose: define data type struct
    IMPLICIT NONE
    TYPE::dat
        REAL::x,y
        REAL::u,v,u_old,v_old
        REAL::AP,AE,AN,AS,AW
        !REAL::Pe,Pn,Ps,Pw,Pp,P_old
        REAL::Pp,P_old
        INTEGER::n
        REAL::uw,ue,vn,vs   ! velocity correction terms
    END TYPE dat
CONTAINS
    SUBROUTINE set_xy (strct,dx,dy,nx,ny)
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx,0:ny),INTENT(INOUT)::strct ! data contained from 0:nx-1 where cells 0 and nx-1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        INTEGER ::  i,j,n   ! for do loops and n is counter for cell number
        REAL    ::  xi,yi   ! x and y values for each cell
        n=1                     ! cell number 1
        DO i=1,ny-1             ! 1 to ny-2 for boundary nodes (we only are iterating through the middle values)
            yi = i*dy - dy/2.   ! y coordinate
            DO j=1,nx-1
                xi = j*dx - dx/2.       ! x coordinate
                strct(j,i)%n = n        ! input n node
                strct(j,i)%x = xi       ! x coordinate to strct
                strct(j,i)%y = yi       ! y coordinate to strct
                !strct(i,j)%phi = 0.     ! phi value initialized guess
                !strct(i,j)%phi_new = 0. ! phi value initialized guess for next iteration
                n=n+1                   ! count cell numbers up one
            END DO
        END DO
    ! left boundary
    strct(0,:)%x = 0.
    strct(0,0)%y = 0.
    strct(0,1:nx-1)%y = RESHAPE((/ (i*dy - dy/2. ,i=1,nx-1) /),(/ nx-1/))
    strct(0,nx+0)%y = 1.
    ! bottom boundary
    strct(0,0)%x = 0.
    strct(1:ny-1,0)%x = RESHAPE((/ (i*dx - dx/2. ,i=1,ny-1) /),(/ ny-1/))
    strct(ny+0,0)%x = 1.
    strct(:,0)%y = 0.
    ! right boundary
    strct(ny+0,:)%x = 1.
    strct(ny+0,:)%y = strct(0,:)%y
    ! top boundary
    strct(:,nx+0)%x = strct(:,0)%x
    strct(:,nx+0)%y = 1.
    END SUBROUTINE set_xy

    SUBROUTINE mom_uv(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        REAL    ::  mu = 0.01
        REAL    ::  rho= 1.
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j !loop iterators
        REAL    ::  Omega = 0.5
        DO i=1,nx
            DO j=1,ny
                mdot            =   rho*strct(i  ,j  )%ue*dy ! east face
                strct(i,j)%AE   =   max(-mdot,0.) + mu*dy/dx
                IF (j==1) THEN !BC
                    mdot            =   rho*strct(i  ,j  )%vn*dx ! north face
                    strct(i,j)%AN   =   max(-mdot,0.) + mu*dx*2./dy
                ELSE
                    mdot            =   rho*strct(i  ,j  )%vn*dx ! north face
                    strct(i,j)%AN   =   max(-mdot,0.) + mu*dx/dy
                END IF
                mdot            =   rho*strct(i  ,j  )%uw*dy ! West face
                strct(i,j)%AW   =   max( mdot,0.) + mu*dy/dx
                IF (j==ny) THEN !BC
                    mdot            =   rho*strct(i  ,j  )%vs*dx ! South face
                    strct(i,j)%AS   =   max( mdot,0.) + mu*dx*2./dy
                ELSE
                    mdot            =   rho*strct(i  ,j-1)%v*dy ! South face
                    strct(i,j)%AS   =   max( mdot,0.) + mu*dx/dy
                END IF
                strct(i,j)%AP   =   strct(i,j)%AE + &
                                    strct(i,j)%AN + &
                                    strct(i,j)%AW + &
                                    strct(i,j)%AS
            END DO
        END DO

        
        DO i=2,nx
            DO j=1,ny
                strct(i,j)%u = (1.-Omega)*strct(i,j)%u_old &
                                + &
                                (Omega/strct(i,j)%AP) &
                                * (&
                                    strct(i,j)%AE*strct(i+1,j  )%u +&
                                    strct(i,j)%AN*strct(i  ,j+1)%u +&
                                    strct(i,j)%AW*strct(i-1,j  )%u +&
                                    strct(i,j)%AS*strct(i  ,j-1)%u +&
                                    (strct(i-1,j)%P_old-strct(i+1,j)%P_old)  *&
                                    dx                              &
                                )
            END DO
        END DO

        DO i=1,nx
            DO j=2,ny
                strct(i,j)%v = (1.-Omega)*strct(i,j)%v_old &
                                + &
                                (Omega/strct(i,j)%AP) &
                                * (&
                                    strct(i,j)%AE*strct(i+1,j  )%v +&
                                    strct(i,j)%AN*strct(i  ,j+1)%v +&
                                    strct(i,j)%AW*strct(i-1,j  )%v +&
                                    strct(i,j)%AS*strct(i  ,j-1)%v +&
                                    (strct(i,j-1)%P_old-strct(i,j+1)%P_old)  *&
                                    dy                              &
                                )
            END DO
        END DO
    END SUBROUTINE mom_uv


    SUBROUTINE vel_correction(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        REAL    ::  mu = 0.01
        REAL    ::  rho= 1.
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j !loop iterators
        REAL    ::  Omega = 0.5
        DO i=1,nx
            DO j=1,ny
                strct(i,j)%uw=(strct(i-1,j  )%P_old - strct(i  ,j  )%P_old) *dy/strct(i,j)%AW
                strct(i,j)%ue=(strct(i  ,j  )%P_old - strct(i+1,j  )%P_old) *dy/strct(i,j)%AE
                strct(i,j)%vn=(strct(i  ,j  )%P_old - strct(i  ,j+1)%P_old) *dx/strct(i,j)%AN
                strct(i,j)%vs=(strct(i  ,j-1)%P_old - strct(i  ,j  )%P_old) *dx/strct(i,j)%AS
            END DO
        END DO

        DO i=1,nx
            DO j=1,ny
                strct(i,j)%Pp=strct(i,j)%P_old+(Omega/strct(i,j)%AP)&
                       *(&
                       + strct(i,j)%AE*strct(i+1,j  )%P_old&
                       + strct(i,j)%AW*strct(i-1,j  )%P_old&
                       + strct(i,j)%AN*strct(i  ,j+1)%P_old&
                       + strct(i,j)%AS*strct(i  ,j-1)%P_old&
                       !- S& ! no source terms
                       - strct(i,j)%AP*strct(i,j)%P_old&
                       )
            END DO
        END DO
    END SUBROUTINE vel_correction



END MODULE types
