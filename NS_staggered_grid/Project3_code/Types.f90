MODULE types
    !purpose: define data type struct
    IMPLICIT NONE
    ! Properties of fluid flow
    REAL    ::  Omega = 0.6 ! Relaxation factor
    REAL    ::  OmegaP= 1.4 ! Relaxation factor for pressure correction
    REAL    ::  alpha = 0.4 ! relaxation factor for pressure correction
    REAL    ::  mu = 0.01   ! dynamic viscosity
    REAL    ::  rho= 1.     ! density
    REAL    ::  Convergence = 1.e-34
    REAL    ::  Convergence2= 1.e-14
    TYPE::dat
        REAL::xu,yv,xp,yp
        REAL::u,v,u_old,v_old !u,v is in bottom left corner, or south and west sides of cell
        REAL::APu,AEu,ANu,ASu,AWu,Apv,AEv,ANv,ASv,AWv,APp,AEp,ANp,AWp,ASp
        REAL::P,Pp,P_old
        REAL::S ! source terms
        INTEGER::n
    END TYPE dat
CONTAINS
    SUBROUTINe set_xy (strct,dx,dy,nx,ny)
        real,intent(in)     ::  dx,dy
        integer,intent(in)  ::  nx,ny   ! size of strct in x and y directions 
        type(dat),dimension(0:,0:),intent(inout)::strct ! data contained from 0:nx-1 where cells 0 and nx-1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        integer ::  i,j,n   ! for do loops and n is counter for cell number
        real    ::  xi,yi   ! x and y values for each cell
        ! left boundary
        strct(0,:)%xp= 0.
        strct(0,0)%yp= 0.
        strct(0,1:ny-1)%yp= reshape((/ (i*dy - dy/2. ,i=1,ny-1) /),(/ ny-1/))
        strct(0,ny+0)%yp= 1.
        ! bottom boundary
        strct(0,0)%xp= 0.
        strct(1:nx-1,0)%xp= reshape((/ (i*dx - dx/2. ,i=1,nx-1) /),(/ nx-1/))
        strct(nx+0,0)%xp= 1.
        strct(:,0)%yp= 0.
        ! right boundary
        strct(nx+0,:)%xp= 1.
        strct(nx+0,:)%yp= strct(0,:)%yp
        ! top boundary
        strct(:,ny+0)%xp= strct(:,0)%xp
        strct(:,ny+0)%yp= 1.
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
    SUBROUTINe set_xy2 (strct,dx,dy,nx,ny)
        real,intent(in)     ::  dx,dy
        integer,intent(in)  ::  nx,ny   ! size of strct in x and y directions 
        type(dat),dimension(0:,0:),intent(inout)::strct ! data contained from 0:nx-1 where cells 0 and nx-1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        integer ::  i,j,n   ! for do loops and n is counter for cell number
        real    ::  xi,yi   ! x and y values for each cell
        ! left boundary
        strct(0,:)%xp= 0.
        strct(0,0)%yp= 0.
        strct(0,1:ny-1)%yp= reshape((/ (i*dy - dy/2. ,i=1,ny-1) /),(/ ny-1/))
        strct(0,ny+0)%yp= 1.
        ! bottom boundary
        strct(0,0)%xp= 0.
        strct(1:nx-1,0)%xp= reshape((/ (i*dx - dx/2. ,i=1,nx-1) /),(/ nx-1/))
        strct(nx+0,0)%xp= 10.
        strct(:,0)%yp= 0.
        ! right boundary
        strct(nx+0,:)%xp= 10.
        strct(nx+0,:)%yp= strct(0,:)%yp
        ! top boundary
        strct(:,ny+0)%xp= strct(:,0)%xp
        strct(:,ny+0)%yp= 1.
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
        strct(nx,:)%xu=10.
        strct(nx,:)%yv=10.
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
        INTEGER,INTENT(IN)  ::  nx,ny! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j,iter=0!loop iterators
        REAL    ::  error=1.,error2=1.

        ! mdot and Au values
        DO i=1,nx
            DO j=1,ny
                mdot                =   rho*(strct(i+1,j  )%u_old+strct(i  ,j  )%u_old)/2.*dy ! east face
                strct(i,j)%AEu      =   max(-mdot,0.) + mu*dy/dx
                mdot                =   rho*(strct(i-1,j+1)%v_old+strct(i  ,j+1)%v_old)/2.*dx ! north face
                IF (j==ny) THEN
                    strct(i,j)%ANu      =   max(-mdot,0.) + mu*2.*dx/dy
                ELSE
                    strct(i,j)%ANu      =   max(-mdot,0.) + mu*dx/dy
                END IF
                mdot                =   rho*(strct(i-1,j  )%u_old+strct(i  ,j  )%u_old)/2.*dy ! West face
                strct(i,j)%AWu      =   max( mdot,0.) + mu*dy/dx
                mdot                =   rho*(strct(i-1,j  )%v_old+strct(i  ,j  )%v_old)/2.*dx ! south face
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

        ! mdot and Av values
        DO i=1,nx
            DO j=1,ny
                mdot                =   rho*(strct(i+1,j-1)%u_old+strct(i+1,j  )%u_old)/2.*dy ! east face
                IF (i==nx) THEN
                    strct(i,j)%AEv      =   max(-mdot,0.) + mu*2.*dy/dx
                ELSE
                    strct(i,j)%AEv      =   max(-mdot,0.) + mu*dy/dx
                END IF
                mdot                =   rho*(strct(i  ,j+1)%v_old+strct(i  ,j  )%v_old)/2.*dx ! north face
                strct(i,j)%ANv      =   max(-mdot,0.) + mu*dx/dy
                mdot                =   rho*(strct(i  ,j-1)%u_old+strct(i  ,j  )%u_old)/2.*dy ! West face
                IF (i==1) THEN
                    strct(i,j)%AWv      =   max( mdot,0.) + mu*2.*dy/dx
                ELSE
                    strct(i,j)%AWv      =   max( mdot,0.) + mu*dy/dx
                END IF
                mdot                =   rho*(strct(i  ,j-1)%v_old+strct(i  ,j  )%v_old)/2.*dx ! south face
                strct(i,j)%ASv      =   max( mdot,0.) + mu*dx/dy
                strct(i,j)%APv      =   strct(i,j)%AEv + &
                    strct(i,j)%ANv + &
                    strct(i,j)%AWv + &
                    strct(i,j)%ASv
                strct(i,j)%APv      = strct(i,j)%APv/Omega
            END DO
        END DO

        ! solve u-momentum
        error2 = 1.
        DO iter=1,100000
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
            IF (abs(error - error2)<Convergence) EXIT   ! error stops changing convergence
        END DO
        WRITE(*,*) iter

        ! solve v-momentum
        error2 = 1.
        DO iter=1,100000
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
            IF (abs(error - error2)<Convergence) EXIT   ! error stops changing convergence
        END DO
        WRITE(*,*) iter
    END SUBROUTINE mom_uv


    SUBROUTINE vel_correction(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        INTEGER ::  i,j,iter=0 !loop iterators
        REAL    :: error,error2
        REAL    ::  S_sum

            DO i=1,nx
                DO j=1,ny
                    IF (i==nx) THEN
                        strct(i,j)%AEp      =   0.
                    ELSE 
                        strct(i,j)%AEp      =   rho*dy*dy/strct(i+1,j)%APu
                    END IF
                    IF (j==ny) THEN
                        strct(i,j)%ANp      =   0.
                    ELSE
                        strct(i,j)%ANp      =   rho*dx*dx/strct(i,j+1)%APv
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
                END DO
            END DO

        DO iter=1,100000
            error2=error
            error=0.
            S_sum = 0.
            DO i=1,nx
                DO j=1,ny
                    strct(i,j)%S = &     ! source terms
                          (rho*strct(i+1,j  )%u_old-rho*strct(i,j)%u_old)*dy&
                        + (rho*strct(i  ,j+1)%v_old-rho*strct(i,j)%v_old)*dx
                    strct(i,j)%Pp = strct(i,j)%Pp + (OmegaP/strct(i,j)%APp)&
                        *(&
                        + strct(i,j)%AEp*strct(i+1,j  )%Pp&
                        + strct(i,j)%AWp*strct(i-1,j  )%Pp&
                        + strct(i,j)%ANp*strct(i  ,j+1)%Pp&
                        + strct(i,j)%ASp*strct(i  ,j-1)%Pp&
                        - strct(i,j)%S                       &
                        - strct(i,j)%APp*strct(i  ,j  )%Pp&
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
            IF (abs(error - error2)<Convergence) THEN   ! error stops changing convergence
            !IF (abs(S_sum)<Convergence) THEN
                strct%P_old = strct%P
                EXIT
            END IF
        END DO
        WRITE(*,*) iter,S_sum ! output iterations along with RSS of source term
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
    SUBROUTINE Solve_NS(strct,dx,dy,nx,ny)
        ! requires uniform grid of dx and dy spacing
        REAL,INTENT(IN)     ::  dx,dy
        INTEGER,INTENT(IN)  ::  nx,ny! size of strct in x and y directions 
        TYPE(dat),DIMENSION(0:nx+1,0:ny+1),INTENT(INOUT)::strct ! data contained from 0:nx+1 where cells 0 and nx+1 are boundary nodes (cell volume approaches 0 on boundary nodes)
        REAL    ::  mdot ! temporary value for mass flow values
        INTEGER ::  i,j,iter=0!loop iterators
        REAL    ::  error2=1.,error_RSS=0.

    DO iter=0,20000
        ! step 1 solve discretised momentum equations
        CALL mom_uv(strct,dx,dy,nx,ny)

!WRITE(*,100) ( data(:,i)%u ,i=0,max_yp )
        ! step 2 Solve pressure correction equation
        ! step 3 Correct pressure and velocities
        CALL vel_correction(strct,dx,dy,nx,ny)


        ! step 4 Solve all other discretised transport equations

        ! if no convergence, then iterate
        error2 = error_RSS
        error_RSS = 0.
        DO i=1,nx
            DO j=1,ny
                error_RSS = error_RSS + (strct(i,j)%u-strct(i,j)%u_old)**2
                error_RSS = error_RSS + (strct(i,j)%v-strct(i,j)%v_old)**2
                error_RSS = error_RSS + (strct(i,j)%P-strct(i,j)%P_old)**2
            END DO
        END DO
        error_RSS = sqrt(error_RSS)
        ! reset values
        strct%u_old  = strct%u
        strct%v_old  = strct%v
        strct%P_old  = strct%P
        !WRITE(*,*) "error = ",error_RSS


        ! if converged then stop
        IF (abs(error_RSS-error2) <= Convergence2) THEN
            WRITE(*,*) "converged on iteration and error big loop = ",iter,abs(error_RSS-error2)
            !WRITE(*,100) ( data(:,i)%S,i=0,max_yp )
            EXIT
        ELSE
            WRITE(*,*) "iteration and error big loop = ",iter,abs(error_RSS-error2)
        END IF
    END DO



    END SUBROUTINE Solve_NS



END MODULE types
