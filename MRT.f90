! FOR Multiple relaxation time.
! Specifc model is D3Q15
! Paper and derivation from 
! Dhumieres, Dominique & Ginzburg, Irina & Krafczyk, Manfred & Lallemand, Pierre & Luo, Li-Shi & Bushnell, Dennis. (2002). Multiple-Relaxation-Time Lattice Boltzmann Models in 3D. 

Module MRT
    implicit none

    ! Transformation matrix
    real(8) :: M(0:14, 0:14)
    ! Inverse M
    real(8) :: M_inv(0:14, 0:14)
    ! Diagonal Relaxation matrix
    real(8) :: S(0:14)

    ! Need lattice info from main

    real(8) :: dt
    real(8) :: c_lattice
    real(8) :: cs2

    ! Integer lattice directions for streaming
    integer :: cx(0:14), cy(0:14), cz(0:14)

    ! Discrete velocity vectors used in momentum/forcing
    real(8) :: ex(0:14), ey(0:14), ez(0:14)

    ! Lattice weights
    real(8) :: w(0:14)

    ! Mean/reference density used in MRT equilibrium
    real(8) :: rho_0


contains 

    subroutine initialize_MRT_matrices(tau, S1, S2, S4, S14, rho0_in)
        implicit none

        real(8), intent(in) :: tau, S1, S2, S4, S14, rho0_in

        integer :: c2

        integer :: i, alpha
        real(8) :: norm_sq

        ! Store reference density
        rho_0 = rho0_in
        
        ! M( 0,:) = [ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]
        
        ! M( 1,:) = [-2, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1]

        ! M( 2,:) = [16, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1]

        ! M( 3,:) = [ 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1]

        ! M( 4,:) = [ 0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1]

        ! M( 5,:) = [ 0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1]

        ! M( 6,:) = [ 0,  0,  0, -4,  4,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1]

        ! M( 7,:) = [ 0,  0,  0,  0,  0,  1, -1,  1,  1,  1,  1, -1, -1, -1, -1]
        
        ! M( 8,:) = [ 0,  0,  0,  0,  0, -4,  4,  1,  1,  1,  1, -1, -1, -1, -1]

        ! M( 9,:) = [ 0,  2,  2, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0]

        ! M(10,:) = [ 0,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0]

        ! M(11,:) = [ 0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1]

        ! M(12,:) = [ 0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1]

        ! M(13,:) = [ 0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1, -1,  1]

        ! M(14,:) = [ 0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1]

        ! Ordering from paper
        ! |m> = (rho, e, epsilon, j_x, q_x, j_y, q_y, j_z, q_z, 3p_xx, p_ww, p_xy, p_yz, p_zx, m_xyz)^T

        ! HOWEVER, main code uses diff ordering than paper column wise so we will reconstruct M from lattice
        ! Eq are in paper
        do i = 0, 14

            c2 = cx(i)*cx(i) + cy(i)*cy(i) + cz(i)*cz(i)

            M(0,i)  = 1.0d0
            M(1,i)  = dble(c2 - 2)
            M(2,i)  = 0.5d0*dble(15*c2*c2 - 55*c2 + 32)

            M(3,i)  = dble(cx(i))
            M(4,i)  = 0.5d0*dble(5*c2 - 13) * dble(cx(i))

            M(5,i)  = dble(cy(i))
            M(6,i)  = 0.5d0*dble(5*c2 - 13) * dble(cy(i))

            M(7,i)  = dble(cz(i))
            M(8,i)  = 0.5d0*dble(5*c2 - 13) * dble(cz(i))

            M(9,i)  = dble(3*cx(i)*cx(i) - c2)
            M(10,i) = dble(cy(i)*cy(i) - cz(i)*cz(i))

            M(11,i) = dble(cx(i)*cy(i))
            M(12,i) = dble(cy(i)*cz(i))
            M(13,i) = dble(cx(i)*cz(i))

            M(14,i) = dble(cx(i)*cy(i)*cz(i))

        end do
        

        ! Calcualte M_inverse
        do i = 0, 14
            ! First get the norm squared for row
            norm_sq = 0.0d0
            do alpha = 0, 14
                norm_sq = norm_sq + M(i, alpha) ** 2
            end do
            
            ! Get inverse
            do alpha = 0, 14
                M_inv(alpha, i) = M(i, alpha) / norm_sq               
            end do
            
        end do

        ! Set relaxation times S
        ! For conserved quatntites, these should be fixed
        S(0) = 0.0d0
        S(3) = 0.0d0
        S(5) = 0.0d0
        S(7) = 0.0d0

        ! Energy modes 
        ! S(1) is usually optimized so should be in input
        ! usually between 1 - 1.5, where 1 is more stable and 1.5 is more agressive (more disspative to less)
        ! This physically controls the bulk viscosity
        S(1) = S1 
        S(2) = S2

        !Energy flux, keep these equal for isotropic behavior
        S(4) = S4
        S(6) = S4
        S(8) = S4

        ! Stress tensor
        S(9)  = 1.0d0/tau
        S(10) = 1.0d0/tau
        S(11) = 1.0d0/tau
        S(12) = 1.0d0/tau
        S(13) = 1.0d0/tau
    
        ! 3rd order moment
        S(14) = S14

    end subroutine

    subroutine MRT_collision(f_in, f_out, ux, uy, uz, rho)
        implicit none
        real(8), intent(in)  :: f_in(0:14)
        real(8), intent(out) :: f_out(0:14)
        real(8), intent(in)  :: ux, uy, uz, rho
        real(8) :: m_temp(0:14), m_eq(0:14)
        integer :: i, j
        ! Before everyting need to set relaxation rates only once
            ! Usually this is predetermiend but for basic MRT, denisty (mass) and momenta are conserved (predetermined to 0)
            ! initialize is done in Main 

        ! 1. Transform f into moment space with M
        ! 2. Get equilibrium moments
        ! 3. Compute Collision in moment space
        ! 4. Transform m back to f (moment to velocity space)

        ! 1.
        do i = 0, 14
            m_temp(i) = 0.0d0
            do j = 0, 14
                m_temp(i) = m_temp(i) + M(i, j) * f_in(j)
            end do
        end do

        ! 2. 
        call MRT_equilibrium_moments(ux, uy, uz, rho, m_eq)

        ! 3. 
        ! m_new = m - S *( m - m_eq)
        do i = 0, 14
            m_temp(i) = m_temp(i) - S(i) * (m_temp(i) - m_eq(i))
        end do

        !4. 
        do j = 0, 14
            f_out(j) = 0.0d0
            do i = 0, 14
                f_out(j) = f_out(j) + M_inv(j, i) * m_temp(i)
            end do
        end do


    end subroutine MRT_collision

    subroutine MRT_IMB_collision(f_in, f_out, os_out, ux, uy, uz, rho, upx, upy, upz, B)
        implicit none

        real(8), intent(in)  :: f_in(0:14)
        real(8), intent(out) :: f_out(0:14)
        real(8), intent(out) :: os_out(0:14)

        real(8), intent(in) :: ux, uy, uz, rho
        real(8), intent(in) :: upx, upy, upz
        real(8), intent(in) :: B

        real(8) :: m_old(0:14)
        real(8) :: m_post(0:14)
        real(8) :: m_new(0:14)

        real(8) :: m_eq_f(0:14)
        real(8) :: m_eq_p(0:14)

        real(8) :: os_mom(0:14)

        integer :: i, j

        ! 1. Transform incoming distributions to moment space
        do i = 0, 14
            m_old(i) = 0.0d0
            do j = 0, 14
                m_old(i) = m_old(i) + M(i,j) * f_in(j)
            end do
        end do

        ! 2. Fluid equilibrium moments:
        !    m_eq_f = m_eq(rho, u_fluid)
        call MRT_equilibrium_moments(ux, uy, uz, rho, m_eq_f)

        ! 3. Particle equilibrium moments:
        !    m_eq_p = m_eq(rho, U_particle)
        call MRT_equilibrium_moments(upx, upy, upz, rho, m_eq_p)

        ! 4. MRT collision in moment space
        !
        !    m_post = m_old - S * (m_old - m_eq_f)
        do i = 0, 14
            m_post(i) = m_old(i) - S(i) * (m_old(i) - m_eq_f(i))
        end do

        ! 5. IMB correction in moment space
        !
        !    m_new = m_post + B * (m_eq_p - m_eq_f)
        do i = 0, 14
            m_new(i) = m_post(i) + B * (m_eq_p(i) - m_eq_f(i))
        end do

        ! 6. os term for hydrodynamic force
        !
        ! Current code uses:
        !
        !   os = collision_delta + peq - feq
        !
        ! MRT equivalent:
        !
        !   os_mom = (m_post - m_old) + (m_eq_p - m_eq_f)
        !
        ! Notice B is NOT included here because the main code already
        ! multiplies the force by B.
        do i = 0, 14
            os_mom(i) = (m_post(i) - m_old(i)) + (m_eq_p(i) - m_eq_f(i))
        end do

        ! 7. Transform corrected moments back to velocity space
        do j = 0, 14
            f_out(j)  = 0.0d0
            os_out(j) = 0.0d0

            do i = 0, 14
                f_out(j)  = f_out(j)  + M_inv(j,i) * m_new(i)
                os_out(j) = os_out(j) + M_inv(j,i) * os_mom(i)
            end do
        end do

    end subroutine MRT_IMB_collision

    ! Following Guo forcing in moment space
    subroutine MRT_collision_body_forcing(f_in, f_out, ux, uy, uz, rho, Fx, Fy, Fz)
        implicit none
        real(8), intent(in)  :: f_in(0:14)
        real(8), intent(out) :: f_out(0:14)
        real(8), intent(out)  :: ux, uy, uz, rho ! this is out here bc equilibrium needs updated values from forcing
        real(8) :: m_temp(0:14), m_eq(0:14)

        ! Force density comp
        real(8), intent(in) :: Fx, Fy, Fz

        real(8) :: Rg(0:14)   ! Guo forcing in vel space
        real(8) :: Cmom(0:14) !moment spaceing forcing term

        real(8) :: momx, momy, momz
        real(8) :: e_dot_u, e_dot_F, u_dot_F

        integer :: i, j

        ! 1. Compute density and raw momentum from f_in
        ! 2. Compute velocity with Guo half-force correction
        ! 3. Transform f into moment space with M
        ! 4. Get equilibrium moments with force corrected velocity
        ! 5. Compute Guo forcing in velocity space
        ! 6. Transform force to moment space
        ! 7. Compute Collision in moment space w/ forcing
        ! 8. Transform m back to f (moment to velocity space)

        ! 1. 
        rho  = 0.0d0
        momx = 0.0d0
        momy = 0.0d0
        momz = 0.0d0

        do i = 0, 14
            rho  = rho  + f_in(i)
            momx = momx + f_in(i)*ex(i)
            momy = momy + f_in(i)*ey(i)
            momz = momz + f_in(i)*ez(i)
        end do


        ! 2.  
        ! With guo half force correction
        !. rho * u = sum_i f_i e_i + 0.5 *dt * F
        ux = (momx + 0.5d0*dt*Fx)/rho
        uy = (momy + 0.5d0*dt*Fy)/rho
        uz = (momz + 0.5d0*dt*Fz)/rho
        
        ! 3. 
        ! Same as normal
        do i = 0, 14
            m_temp(i) = 0.0d0
            do j = 0, 14
                m_temp(i) = m_temp(i) + M(i,j)*f_in(j)
            end do
        end do 

        ! 4. 
        ! Same but with new force corrected vel
        call MRT_equilibrium_moments(ux, uy, uz, rho, m_eq)

        ! 5. 
        u_dot_F = ux*Fx + uy*Fy + uz*Fz

        do i = 0, 14

            e_dot_u = dble(cx(i))*ux + dble(cy(i))*uy + dble(cz(i))*uz
            e_dot_F = dble(cx(i))*Fx + dble(cy(i))*Fy + dble(cz(i))*Fz

            Rg(i) = w(i) * ((e_dot_F - u_dot_F)/cs2 + (e_dot_u*e_dot_F)/(cs2*cs2) )

        end do

        ! 6.
        do i = 0, 14
            Cmom(i) = 0.0d0
            do j = 0, 14
                Cmom(i) = Cmom(i) + M(i,j)*Rg(j)
            end do
        end do         

        ! 7. 
        ! Collision in moment space with forcing
        do i = 0, 14
            m_temp(i) = m_temp(i) - S(i)*(m_temp(i) - m_eq(i)) + dt*(1.0d0 - 0.5d0*S(i))*Cmom(i)
        end do

        ! 8. Back to vel space
        do j = 0, 14
            f_out(j) = 0.0d0
            do i = 0, 14
                f_out(j) = f_out(j) + M_inv(j,i)*m_temp(i)
            end do
        end do


    end subroutine MRT_collision_body_forcing
 
    subroutine MRT_equilibrium_moments(ux, uy, uz, rho, m_eq_out)
        implicit none
        real(8), intent(in) :: ux, uy, uz, rho
        real(8), intent(out) :: m_eq_out(0:14)
        real(8) :: jx, jy, jz
        real(8) :: j_sqr

        ! Momentum comp
        jx = rho * ux
        jy = rho * uy
        jz = rho * uz

        j_sqr = jx*jx + jy*jy + jz*jz

        ! Equilibrium moments from paper

        ! Conserved moments so these are exact
        m_eq_out(0) = rho 
        m_eq_out(3) = jx 
        m_eq_out(5) = jy
        m_eq_out(7) = jz

        ! Energy modes
        m_eq_out(1) = -rho + (1.0d0/rho_0) * j_sqr
        m_eq_out(2) = -rho
        ! m_eq_out(2) = rho - 5.0d0*j_sqr/rho_0
        
        ! Energy flux modes
        m_eq_out(4) = -(7.0d0/3.0d0) * jx
        m_eq_out(6) = -(7.0d0/3.0d0) * jy
        m_eq_out(8) = -(7.0d0/3.0d0) * jz

        ! Stress tensor eq
        ! This factor of 3 is redundant but it is written this way explicitly to not confuse eq.
        m_eq_out(9)  = 3.0d0 * (1.0d0/(3.0d0*rho_0)) * (2.0d0 * jx * jx - (jy * jy + jz * jz))
        m_eq_out(10) = (1.0d0/rho_0) * (jy*jy - jz*jz)

        m_eq_out(11) = (1.0d0/rho_0) * jx * jy 
        m_eq_out(12) = (1.0d0/rho_0) * jy * jz 
        m_eq_out(13) = (1.0d0/rho_0) * jx * jz
        
        m_eq_out(14) = 0

    end subroutine MRT_equilibrium_moments


    subroutine initialize_MRT_lattice(cx_in, cy_in, cz_in, w_in)
        ! This is for getting information from main that is required by MRT forcing
        ! Very future goal, it would be better if seprate module called D3Q15 lattice so both modules can read
        implicit none

        integer, intent(in) :: cx_in(0:14)
        integer, intent(in) :: cy_in(0:14)
        integer, intent(in) :: cz_in(0:14)

        real(8), intent(in) :: w_in(0:14)

        integer :: i

        dt        = 1.0d0
        c_lattice = 1.0d0
        cs2       = 1.0d0/3.0d0

        do i = 0, 14
            cx(i) = cx_in(i)
            cy(i) = cy_in(i)
            cz(i) = cz_in(i)

            ex(i) = dble(cx(i))
            ey(i) = dble(cy(i))
            ez(i) = dble(cz(i))

            w(i) = w_in(i)
        end do
    end subroutine initialize_MRT_lattice
 
end module MRT
