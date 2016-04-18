MODULE types
    !purpose: define data type struct
    IMPLICIT NONE
    ! Properties of fluid flow
    REAL    ::  Omega = 0.5 ! Relaxation factor
    REAL    ::  alpha = 0.04! relaxation factor for pressure correction
    REAL    ::  mu = 0.01   ! dynamic viscosity
    REAL    ::  rho= 1.     ! density
    REAL    ::  Convergence = 1.e-7
    TYPE::dat
        REAL::xu,yv,xp,yp
        REAL::u,v,u_old,v_old,u_orig,v_orig !u,v is in bottom left corner, or south and west sides of cell
        REAL::APu,AEu,ANu,ASu,AWu,Apv,AEv,ANv,ASv,AWv,APp,AEp,ANp,AWp,ASp
        !REAL::Pe,Pn,Ps,Pw,Pp,P_old
        REAL::P,Pp,P_old
        REAL::S ! source terms
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
        INTEGER,INTENT(IN)  ::  nx,ny!,iter   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        !REAL    ::  mu = 0.01
        !REAL    ::  rho= 1.
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j,iter2=0!loop iterators
        REAL    ::  error=1.,error2=1.
        !REAL    ::  Omega = 0.5
        ! DO i=1,nx
        !     DO j=1,ny
        !IF (iter2 == 0) THEN
        DO i=1,nx
            DO j=1,ny
                mdot                =   rho*(strct(i+1,j  )%u_orig+strct(i  ,j  )%u_orig)/2.*dy ! east face
                strct(i,j)%AEu      =   max(-mdot,0.) + mu*dy/dx
                mdot                =   rho*(strct(i-1,j+1)%v_orig+strct(i  ,j+1)%v_orig)/2.*dx ! north face
                IF (j==ny) THEN
                    strct(i,j)%ANu      =   max(-mdot,0.) + mu*2.*dx/dy
                ELSE
                    strct(i,j)%ANu      =   max(-mdot,0.) + mu*dx/dy
                END IF
                mdot                =   rho*(strct(i-1,j  )%u_orig+strct(i  ,j  )%u_orig)/2.*dy ! West face
                strct(i,j)%AWu      =   max( mdot,0.) + mu*dy/dx
                mdot                =   rho*(strct(i-1,j  )%v_orig+strct(i  ,j  )%v_orig)/2.*dx ! south face
                IF (j==1) THEN
                    strct(i,j)%ASu      =   max( mdot,0.) + mu*2.*dx/dy
                ELSE
                    strct(i,j)%ASu      =   max( mdot,0.) + mu*dx/dy
                END IF
                strct(i,j)%APu      =   strct(i,j)%AEu + &
                    strct(i,j)%ANu + &
                    strct(i,j)%AWu + &
                    strct(i,j)%ASu
                strct(i,j)%APu      = strct(i,j)%APu/Omega
            END DO
        END DO
        !WRITE( *,100) ( strct(:,i)%APu,i=0,11 )
        !100 FORMAT (12 ES16.7)

        

        error2 = 1.
        DO iter2=1,100000
            error2=error
            error = 0.
            DO i=2,nx
                DO j=1,ny
                    strct(i,j)%u = (1.-Omega)*strct(i,j)%u_old &
                        + &
                        (1./strct(i,j)%APu) &
                        * (&
                        strct(i  ,j  )%AEu*strct(i+1,j  )%u +   &
                        strct(i  ,j  )%ANu*strct(i  ,j+1)%u +   &
                        strct(i  ,j  )%AWu*strct(i-1,j  )%u +   &
                        strct(i  ,j  )%ASu*strct(i  ,j-1)%u +   &
                        (strct(i-1,j)%P_old-strct(i  ,j)%P_old) &
                        *dy                                      &
                        )

                    error = error + (strct(i,j)%u - strct(i,j)%u_old)**2
                END DO
            END DO
            error=sqrt(error)
            IF (abs(error - error2)<Convergence) EXIT
        END DO
        WRITE(*,*) iter2

        DO i=1,nx
            DO j=1,ny
                mdot                =   rho*(strct(i+1,j-1)%u_orig+strct(i+1,j  )%u_orig)/2.*dy ! east face
                IF (i==nx) THEN
                    strct(i,j)%AEv      =   max(-mdot,0.) + mu*2.*dy/dx
                ELSE
                    strct(i,j)%AEv      =   max(-mdot,0.) + mu*dy/dx
                END IF
                mdot                =   rho*(strct(i  ,j+1)%v_orig+strct(i  ,j  )%v_orig)/2.*dx ! north face
                strct(i,j)%ANv      =   max(-mdot,0.) + mu*dx/dy
                mdot                =   rho*(strct(i  ,j-1)%u_orig+strct(i  ,j  )%u_orig)/2.*dy ! West face
                IF (i==1) THEN
                    strct(i,j)%AWv      =   max( mdot,0.) + mu*2.*dy/dx
                ELSE
                    strct(i,j)%AWv      =   max( mdot,0.) + mu*dy/dx
                END IF
                mdot                =   rho*(strct(i  ,j-1)%v_orig+strct(i  ,j  )%v_orig)/2.*dx ! south face
                strct(i,j)%ASv      =   max( mdot,0.) + mu*dx/dy
                strct(i,j)%APv      =   strct(i,j)%AEv + &
                    strct(i,j)%ANv + &
                    strct(i,j)%AWv + &
                    strct(i,j)%ASv
                strct(i,j)%APv      = strct(i,j)%APv/Omega
            END DO
        END DO

        error2 = 1.
        DO iter2=1,100000
            error2=error
            error = 0.
            DO i=1,nx
                DO j=2,ny
                    strct(i,j)%v = (1.-Omega)*strct(i,j)%v_old &
                        + &
                        (1./strct(i,j)%APv) &
                        * (&
                        strct(i  ,j  )%AEv*strct(i+1,j  )%v +&
                        strct(i  ,j  )%ANv*strct(i  ,j+1)%v +&
                        strct(i  ,j  )%AWv*strct(i-1,j  )%v +&
                        strct(i  ,j  )%ASv*strct(i  ,j-1)%v +&
                        (strct(i,j-1)%P_old-strct(i  ,j)%P_old)  *&
                        dy                              &
                        )
                    error = error + (strct(i,j)%v - strct(i,j)%v_old)**2
                END DO
            END DO
            error=sqrt(error)
            IF (abs(error - error2)<Convergence) EXIT
        END DO
        WRITE(*,*) iter2
    END SUBROUTINE mom_uv


    SUBROUTINE vel_correction(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny!,iter   ! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        !REAL    ::  mu = 0.01
        !REAL    ::  rho= 1.
        !REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j,iter2=0 !loop iterators
        REAL    :: error,error2
        REAL    ::  S_sum
        !REAL    ::  Omega = 0.5

        !IF (iter == 0) THEN
            DO i=1,nx
                DO j=1,ny
                    IF (i==nx) THEN
                        strct(i,j)%AEp      =   0.
                    ELSE 
                        strct(i,j)%AEp      =   rho*dy*dy/strct(i,j)%APu
                    END IF
                    IF (j==ny) THEN
                        strct(i,j)%ANp      =   0.
                    ELSE
                        strct(i,j)%ANp      =   rho*dx*dx/strct(i,j)%APv
                    END IF
                    IF (i==1) THEN
                        strct(i,j)%AWp      =   0.
                    ELSE
                        strct(i,j)%AWp      =   rho*dy*dy/strct(i,j)%APu
                    END IF
                    IF (j==1) THEN
                        strct(i,j)%ASp      =   0.
                    ELSE
                        strct(i,j)%ASp      =   rho*dx*dx/strct(i,j)%APv
                    END IF
                    strct(i,j)%APp          =   strct(i,j)%AEp + &
                        strct(i,j)%ANp + &
                        strct(i,j)%AWp + &
                        strct(i,j)%ASp
                    strct(i,j)%APp      = strct(i,j)%APp/(Omega) ! overrelax because it is linear
                END DO
            END DO
        !END IF

        DO iter2=1,100000
            error2=error
            error=0.
            S_sum = 0.
            DO i=1,nx
                DO j=1,ny
                    strct(i,j)%S = &
                          (rho*strct(i+1,j  )%u_old-rho*strct(i,j)%u_old)*dy&     ! source terms
                        + (rho*strct(i  ,j+1)%v_old-rho*strct(i,j)%v_old)*dx
                    strct(i,j)%Pp = (1./strct(i,j)%APp)&
                        *(&
                        + strct(i,j)%AEp*strct(i+1,j  )%P_old&
                        + strct(i,j)%AWp*strct(i-1,j  )%P_old&
                        + strct(i,j)%ANp*strct(i  ,j+1)%P_old&
                        + strct(i,j)%ASp*strct(i  ,j-1)%P_old&
                        - strct(i,j)%S                       &
                        - strct(i,j)%APp*strct(i  ,j  )%P_old&
                        )

                    strct(i,j)%P=strct(i,j)%P_old+alpha*strct(i,j)%Pp
                    error = error + (strct(i,j)%P - strct(i,j)%P_old)**2
                    S_sum = S_sum + strct(i,j)%S**2
                    IF (ISNAN(strct(i,j)%Pp)) THEN
                        WRITE(*,*) "error on ",i,j
                        STOP
                    END IF
                END DO
            END DO
            IF (abs(error - error2)<Convergence) THEN
            !IF (abs(S_sum)<Convergence) THEN
                !WRITE(*,*) strct%S
                EXIT
            END IF
        END DO
        WRITE(*,*) iter2,S_sum
        DO i=2,nx
            DO j=1,ny
                strct(i,j)%u=strct(i,j)%u + (strct(i-1,j  )%Pp - strct(i  ,j  )%Pp) *dy/strct(i,j)%APu
            END DO
        END DO
        DO i=1,nx
            DO j=2,ny
                strct(i,j)%v=strct(i,j)%v + (strct(i  ,j-1)%Pp - strct(i  ,j  )%Pp) *dx/strct(i,j)%APv
            END DO
        END DO

    END SUBROUTINE vel_correction



END MODULE types
