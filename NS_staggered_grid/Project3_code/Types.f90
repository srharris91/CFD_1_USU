MODULE types
    !purpose: define data type struct
    IMPLICIT NONE
    ! Properties of fluid flow
    REAL    ::  Omega = 0.5 ! Relaxation factor
    REAL    ::  mu = 0.01   ! dynamic viscosity
    REAL    ::  rho= 1.     ! density
    TYPE::dat
        REAL::xu,yv,xp,yp
        REAL::u,v,u_old,v_old !u,v is in bottom left corner, or south and west sides of cell
        REAL::AP,AE,AN,AS,AW
        !REAL::Pe,Pn,Ps,Pw,Pp,P_old
        REAL::P,Pp,P_old
        INTEGER::n
        !REAL::uw,ue,vn,vs   ! velocity correction terms
    END TYPE dat
CONTAINS
    SUBROUTINE set_xy (strct,dx,dy,nx,ny)
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx,0:ny),INTENT(INOUT)::strct ! data contained from 0:nx-1 where cells 0 and nx-1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        INTEGER ::  i,j,n   ! for do loops and n is counter for cell number
        REAL    ::  xi,yi   ! x and y values for each cell
        ! left boundary
        strct(0,:)%xp= 0.
        strct(0,0)%yp= 0.
        strct(0,1:nx-1)%yp= RESHAPE((/ (i*dy - dy/2. ,i=1,nx-1) /),(/ nx-1/))
        strct(0,nx+0)%yp= 1.
        ! bottom boundary
        strct(0,0)%xp= 0.
        strct(1:ny-1,0)%xp= RESHAPE((/ (i*dx - dx/2. ,i=1,ny-1) /),(/ ny-1/))
        strct(ny+0,0)%xp= 1.
        strct(:,0)%yp= 0.
        ! right boundary
        strct(ny+0,:)%xp= 1.
        strct(ny+0,:)%yp= strct(0,:)%yp
        ! top boundary
        strct(:,nx+0)%xp= strct(:,0)%xp
        strct(:,nx+0)%yp= 1.
        ! set xu and yv to similar values
        !top
        strct(:,ny)%xu=strct(:,ny)%xp-dx/2.
        strct(:,ny)%yv=strct(:,ny)%yp-dy/2.
        !bottom
        strct(:,0)%xu=strct(:,0)%xp-dx/2.
        strct(:,0)%yv=strct(:,0)%yp-dy/2.
        !left
        strct(0:1,:)%xu=0.
        strct(0:1,:)%yv=0.
        !right
        strct(nx,:)%xu=1.
        strct(nx,:)%yv=1.
        n=1                     ! cell number 1
        DO i=1,ny-1             ! 1 to ny-2 for boundary nodes (we only are iterating through the middle values)
            yi = i*dy - dy/2.   ! y coordinate
            DO j=1,nx-1
                xi = j*dx - dx/2.       ! x coordinate
                strct(j,i)%n = n        ! input n node
                strct(j,i)%xp= xi       ! x coordinate to strct
                strct(j,i)%xu= xi-dx/2. ! x coordinate to strct
                strct(j,i)%yp= yi       ! y coordinate to strct
                strct(j,i)%yv= yi-dy/2. ! y coordinate to strct
                !strct(i,j)%phi = 0.     ! phi value initialized guess
                !strct(i,j)%phi_new = 0. ! phi value initialized guess for next iteration
                n=n+1                   ! count cell numbers up one
            END DO
        END DO
    END SUBROUTINE set_xy

    SUBROUTINE mom_uv(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        !REAL    ::  mu = 0.01
        !REAL    ::  rho= 1.
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j,iter !loop iterators
        REAL    ::  error
        !REAL    ::  Omega = 0.5
        ! DO i=1,nx
        !     DO j=1,ny
        DO i=2,nx
            DO j=1,ny
                    mdot                =   rho*(strct(i+1,j  )%u+strct(i  ,j  )%u)/2.*dy ! east face
                    strct(i,j)%AE       =   max(-mdot,0.) + mu*dy/dx
                    mdot                =   rho*(strct(i-1,j+1)%v+strct(i  ,j+1)%v)/2.*dx ! north face
                    strct(i,j)%AN       =   max(-mdot,0.) + mu*dx/dy
                    mdot                =   rho*(strct(i-1,j  )%u+strct(i  ,j  )%u)/2.*dy ! West face
                    strct(i,j)%AW       =   max( mdot,0.) + mu*dy/dx
                    mdot                =   rho*(strct(i-1,j  )%v+strct(i  ,j  )%v)/2.*dx ! south face
                    strct(i,j)%AS       =   max( mdot,0.) + mu*dx/dy
                strct(i,j)%AP           =   strct(i,j)%AE + &
                                            strct(i,j)%AN + &
                                            strct(i,j)%AW + &
                                            strct(i,j)%AS
            END DO
        END DO


        !DO iter=1,100
        DO i=2,nx
            DO j=1,ny
                strct(i,j)%u = (1.-Omega)*strct(i,j)%u_old &
                    + &
                    (Omega/strct(i,j)%AP) &
                    * (&
                    strct(i  ,j  )%AE*strct(i+1,j  )%u +&
                    strct(i  ,j  )%AN*strct(i  ,j+1)%u +&
                    strct(i  ,j  )%AW*strct(i-1,j  )%u +&
                    strct(i  ,j  )%AS*strct(i  ,j-1)%u +&
                    (strct(i-1,j)%P_old-strct(i  ,j)%P_old)  *&
                    dy                              &
                    )
            END DO
        END DO
        !error = 0.
        !DO i=1,nx
        !DO j=1,ny
        !error = error + (strct(i,j)%u-strct(i,j)%u_old)**2
        !END DO
        !END DO
        !error= sqrt(error)
        !! if converged then stop
        !IF (error < 0.0000001) THEN
        !EXIT
        !ELSE
        !WRITE(*,*) "iteration and error = ",iter,error
        !END IF
        !END DO

        !DO iter=1,10000
        !DO i=1,nx
        !DO j=2,ny
        !strct(i,j)%v = (1.-Omega)*strct(i,j)%v_old &
        !+ &
        !(Omega/strct(i,j)%AP) &
        !* (&
        !strct(i,j)%AE*strct(i+1,j  )%v +&
        !strct(i,j)%AN*strct(i  ,j+1)%v +&
        !strct(i,j)%AW*strct(i-1,j  )%v +&
        !strct(i,j)%AS*strct(i  ,j-1)%v +&
        !(strct(i,j-1)%P_old-strct(i,j  )%P_old)  *&
        !dx                              &
        !)
        !END DO
        !END DO
        !error = 0.
        !DO i=1,nx
        !DO j=1,ny
        !error = error + (strct(i,j)%u-strct(i,j)%u_old)**2
        !END DO
        !END DO
        !error= sqrt(error)
        !! if converged then stop
        !IF (error < 0.0000001) THEN
        !EXIT
        !ELSE
        !WRITE(*,*) "iteration and error = ",iter,error
        !END IF
        !END DO
    END SUBROUTINE mom_uv


    SUBROUTINE vel_correction(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        !REAL    ::  mu = 0.01
        !REAL    ::  rho= 1.
        !REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j !loop iterators
        !REAL    ::  Omega = 0.5
        DO i=1,nx
            DO j=1,ny
                strct(i,j)%Pp = (Omega/strct(i,j)%AP)&
                    *(&
                    + strct(i,j)%AE*strct(i+1,j  )%P&
                    + strct(i,j)%AW*strct(i-1,j  )%P&
                    + strct(i,j)%AN*strct(i  ,j+1)%P&
                    + strct(i,j)%AS*strct(i  ,j-1)%P&
                    - strct(i,j)%AP*strct(i,j)%P_old&
                    )

                strct(i,j)%P=strct(i,j)%P_old+strct(i,j)%Pp !(Omega/strct(i,j)%AP)&
                !*(&
                !+ strct(i,j)%AE*strct(i+1,j  )%P&
                !+ strct(i,j)%AW*strct(i-1,j  )%P&
                !+ strct(i,j)%AN*strct(i  ,j+1)%P&
                !+ strct(i,j)%AS*strct(i  ,j-1)%P&
                !- S& ! no source terms
                !- strct(i,j)%AP*strct(i,j)%P_old&
                !)
                if (ISNAN(strct(i,j)%u)) THEN
                    WRITE(*,*) "NaN on ",i,j
                    STOP
                END IF
            END DO
        END DO
        DO i=2,nx
            DO j=1,ny
                strct(i  ,j  )%u=(strct(i-1,j  )%Pp - strct(i  ,j  )%Pp) *dy/strct(i,j)%AW
            END DO
        END DO
        !strct(i  ,j  )%ue=(strct(i  ,j  )%P_old - strct(i+1,j  )%P_old) *dy/strct(i,j)%AE
        !strct(i  ,j  )%vn=(strct(i  ,j  )%P_old - strct(i  ,j+1)%P_old) *dx/strct(i,j)%AN
        DO i=1,nx
            DO j=2,ny
                strct(i  ,j  )%v=(strct(i  ,j-1)%Pp - strct(i  ,j  )%Pp) *dx/strct(i,j)%AS
            END DO
        END DO

    END SUBROUTINE vel_correction



END MODULE types
