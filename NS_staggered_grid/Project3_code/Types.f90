MODULE types
    !purpose: define data type struct
    IMPLICIT NONE
    ! Properties of fluid flow
    REAL    ::  Omega = 0.5 ! Relaxation factor
    REAL    ::  mu = 0.01   ! dynamic viscosity
    REAL    ::  rho= 1.     ! density
    TYPE::dat
        REAL::xu,yv,xp,yp
        REAL::u,v,u_old,v_old,u_orig,v_orig !u,v is in bottom left corner, or south and west sides of cell
        REAL::APu,AEu,ANu,ASu,AWu,Apv,AEv,ANv,ASv,AWv,APp,AEp,ANp,AWp,ASp
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

    SUBROUTINE mom_uv(strct,dx,dy,nx,ny,iter)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny,iter   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        !REAL    ::  mu = 0.01
        !REAL    ::  rho= 1.
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j!loop iterators
        REAL    ::  error
        !REAL    ::  Omega = 0.5
        ! DO i=1,nx
        !     DO j=1,ny
        IF (iter == 0) THEN
            DO i=1,nx
                DO j=1,ny
                    mdot                =   rho*(strct(i+1,j  )%u_orig+strct(i  ,j  )%u_orig)/2.*dy ! east face
                    strct(i,j)%AEu      =   max(-mdot,0.) + mu*dy/dx
                    mdot                =   rho*(strct(i-1,j+1)%v_orig+strct(i  ,j+1)%v_orig)/2.*dx ! north face
                    strct(i,j)%ANu      =   max(-mdot,0.) + mu*dx/dy
                    mdot                =   rho*(strct(i-1,j  )%u_orig+strct(i  ,j  )%u_orig)/2.*dy ! West face
                    strct(i,j)%AWu      =   max( mdot,0.) + mu*dy/dx
                    mdot                =   rho*(strct(i-1,j  )%v_orig+strct(i  ,j  )%v_orig)/2.*dx ! south face
                    strct(i,j)%ASu      =   max( mdot,0.) + mu*dx/dy
                    strct(i,j)%APu      =   strct(i,j)%AEu + &
                        strct(i,j)%ANu + &
                        strct(i,j)%AWu + &
                        strct(i,j)%ASu
                END DO
            END DO
        END IF


        !DO iter=1,100
        DO i=2,nx
            DO j=1,ny
                strct(i,j)%u = (1.-Omega)*strct(i,j)%u_old &
                    + &
                    (Omega/strct(i,j)%APu) &
                    * (&
                    strct(i  ,j  )%AEu*strct(i+1,j  )%u +&
                    strct(i  ,j  )%ANu*strct(i  ,j+1)%u +&
                    strct(i  ,j  )%AWu*strct(i-1,j  )%u +&
                    strct(i  ,j  )%ASu*strct(i  ,j-1)%u +&
                    (strct(i-1,j)%P_old-strct(i  ,j)%P_old)  *&
                    dy                              &
                    )
            END DO
        END DO
        IF (iter == 0) THEN
            DO i=1,nx
                DO j=1,ny
                    mdot                =   rho*(strct(i+1,j-1)%u_orig+strct(i+1,j  )%u_orig)/2.*dy ! east face
                    strct(i,j)%AEv      =   max(-mdot,0.) + mu*dy/dx
                    mdot                =   rho*(strct(i  ,j+1)%v_orig+strct(i  ,j  )%v_orig)/2.*dx ! north face
                    strct(i,j)%ANv      =   max(-mdot,0.) + mu*dx/dy
                    mdot                =   rho*(strct(i  ,j-1)%u_orig+strct(i  ,j  )%u_orig)/2.*dy ! West face
                    strct(i,j)%AWv      =   max( mdot,0.) + mu*dy/dx
                    mdot                =   rho*(strct(i  ,j-1)%v_orig+strct(i  ,j  )%v_orig)/2.*dx ! south face
                    strct(i,j)%ASv      =   max( mdot,0.) + mu*dx/dy
                    strct(i,j)%APv      =   strct(i,j)%AEv + &
                        strct(i,j)%ANv + &
                        strct(i,j)%AWv + &
                        strct(i,j)%ASv
                END DO
            END DO
        END IF
        DO i=1,nx
            DO j=2,ny
                strct(i,j)%v = (1.-Omega)*strct(i,j)%v_old &
                    + &
                    (Omega/strct(i,j)%APv) &
                    * (&
                    strct(i  ,j  )%AEv*strct(i+1,j  )%v +&
                    strct(i  ,j  )%ANv*strct(i  ,j+1)%v +&
                    strct(i  ,j  )%AWv*strct(i-1,j  )%v +&
                    strct(i  ,j  )%ASv*strct(i  ,j-1)%v +&
                    (strct(i,j-1)%P_old-strct(i  ,j)%P_old)  *&
                    dy                              &
                    )
            END DO
        END DO
    END SUBROUTINE mom_uv


    SUBROUTINE vel_correction(strct,dx,dy,nx,ny,iter)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny,iter   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        !REAL    ::  mu = 0.01
        !REAL    ::  rho= 1.
        !REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j !loop iterators
        !REAL    ::  Omega = 0.5

        IF (iter == 0) THEN
            DO i=1,nx
                DO j=1,ny
                    !IF (i==nx) THEN
                    !strct(i,j)%AEp      =   0.
                    !ELSE 
                        strct(i,j)%AEp      =   rho*dy*dy/strct(i,j)%AEu
                        !END IF
                        !IF (j==ny) THEN
                        !strct(i,j)%ANp      =   0.
                        !ELSE
                        strct(i,j)%ANp      =   rho*dx*dx/strct(i,j)%ANv
                        !END IF
                        !IF (i==1) THEN
                        !strct(i,j)%AWp      =   0.
                        !ELSE
                        strct(i,j)%AWp      =   rho*dy*dy/strct(i,j)%AWu
                        !END IF
                        !IF (j==1) THEN
                        !strct(i,j)%ASp      =   0.
                        !ELSE
                        strct(i,j)%ASp      =   rho*dx*dx/strct(i,j)%ASv
                        !END IF
                    strct(i,j)%APp          =   strct(i,j)%AEp + &
                        strct(i,j)%ANp + &
                        strct(i,j)%AWp + &
                        strct(i,j)%ASp
                END DO
            END DO
        END IF

        DO i=1,nx
            DO j=1,ny
                strct(i,j)%Pp = (Omega/strct(i,j)%APp)&
                    *(&
                    + strct(i,j)%AEp*strct(i+1,j  )%P_old&
                    + strct(i,j)%AWp*strct(i-1,j  )%P_old&
                    + strct(i,j)%ANp*strct(i  ,j+1)%P_old&
                    + strct(i,j)%ASp*strct(i  ,j-1)%P_old&
                    - strct(i,j)%APp*strct(i  ,j  )%P_old&
                    )

                strct(i,j)%P=strct(i,j)%P_old+strct(i,j)%Pp
            END DO
        END DO
        DO i=2,nx
            DO j=1,ny
                strct(i,j)%u=strct(i,j)%u + (strct(i-1,j  )%Pp - strct(i  ,j  )%Pp) *dy/strct(i,j)%AWu
            END DO
        END DO
        DO i=1,nx
            DO j=2,ny
                strct(i,j)%v=strct(i,j)%v + (strct(i  ,j-1)%Pp - strct(i  ,j  )%Pp) *dx/strct(i,j)%ASv
            END DO
        END DO

    END SUBROUTINE vel_correction



END MODULE types
