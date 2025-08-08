! T-PFC (Phase Field Crystal) model, for single-component 2D system with transverse interactions
! 2D hexagonal/triangular lattice structure
! double precision, data output

! pseudo-spetral method (using FFTW3)
! use the exponential propagation scheme of Cross et al., Chaos 4, 607 (1994)
! with minor corrections of expanding around small alpha_1 to get cf_1 and cf2_1 when alpha_1 -> 0
! with the transient time

!# define ini_homogeneous ! homogeneous initial condition
# define ini_nucleus_crystal_cir ! initial crystalline nucleus (circular crystallite) in homogeneous medium

!# define save_conf

# define open_MP ! use OpenMP for multi-threading


module global_variables
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  
  save ! save all the variables in the module

  integer, parameter :: ly=512, lx=512
  
  real(C_DOUBLE), pointer :: psi(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: psiq(:,:)
  real(8) :: qx(lx),qy(ly/2+1)
  real(8) :: eps,g,q02,dx,dy,psi0,time,pi
  real(8) :: alpha1,alpha2
  integer :: i_alpha1,i_alpha2
  double precision, allocatable :: q2(:,:),exp_1(:,:),cf_1(:,:),cf2_1(:,:)
  real(8) :: TOL,err_s ! for iteration of predictor-corrector
  integer :: nmbr_eval

! for FFTW parameters
  type(C_PTR) :: plan2d_f,plan2d_b
  real(8) :: scale2d_b

  character(:), allocatable :: filename ! to be automatically allocated

end module global_variables


Program Tpfc_2D
  use global_variables
  implicit none

  type(C_PTR) :: p_alloc_psi,p_alloc_psiq ! for allocating aligned memory in fftw
  real(8) :: psi0_sol,psi0_liq
  real(8) :: qx0,qy0,q0
  real(8) :: dt,dti,dtime
  real(8) :: tmax,time0,t_i
  real(8) :: alpha_1,alpha_dt
  real(8) :: A0a,q_1,q_2,q_3
  real(8) :: noise,noise0
  real(8) :: ran2,x,y
  integer :: i,j,iter,idum,istep
  integer :: nend,raout,nout_conf,n_i,n_0
  integer :: n_dx,iter0,i_psi0
  
# ifdef ini_nucleus_crystal_cir
  ! initial crystalline nucleus (circular crystallite) in homogeneous medium
  real(8) :: r0,r02,psi0_s,psi0_l,noise_n
# endif
  
# ifdef open_MP
  ! use multi-threaded fftw with openMP
  integer iret,nthreads,omp_get_max_threads
# endif

# ifdef ini_homogeneous
  ! for homogeneous initial condition
  filename='th512eps01g05al3al0psi00' ! 512x512,eps=0.1,g=0.5,alpha1=3,alpha2=0,psi0=0,noise=noise0=0.1; fragmented multigrains
# endif

# ifdef ini_nucleus_crystal_cir
  ! initial crystalline nucleus (circular crystallite) in homogeneous medium
  filename='tn512eps01g05al3al0_r83' ! 512x512, eps=0.1, g=0.5, alpha1=3, alpha2=0, r0=3*lx/8, psi0_s=-0.089431,psi0_l=-0.11847, noise=noise_n=0
                                     ! surface cusp instability of self-rotating grain, and grain fragmentation at large t
# endif
  
  eps=0.1d0
  g=0.5d0
  
! parameters for transverse interaction
  alpha1=3
  if(abs(alpha1) < 1.e-8) then
     i_alpha1=0 ! if alpha1=0
  else
     i_alpha1=1 ! if alpha1 nonzero
  endif
  
  alpha2=0
  if(abs(alpha2) < 1.e-8) then
     i_alpha2=0 ! if alpha2=0
  else
     i_alpha2=1 ! if alpha2 nonzero
  endif
  
  pi=DACOS(-1.0d0)

  q0=1
  q02=1
  qx0=(sqrt(3.d0)/2)*q0
  qy0=q0

! setup of average density psi0
! psi0=-0.0894314330862110,-0.1184741052584242 for solid-liquid coexistence, if eps=0.1, g=0.5 (from original PFC phase diagram)
  psi0_sol=-0.089431d0
  psi0_liq=-0.11847d0

  psi0=0
  
  n_dx=8 ! n_dx: number of grid points per lattice period
  dy=2*pi/(n_dx*qy0)
  dx=dy
  
  noise=0.1d0
  !noise=0
  noise0=noise
  
# ifdef ini_nucleus_crystal_cir
  ! initial crystalline nucleus (circular crystallite) in homogeneous medium
  r0=3*lx/8 ! radius of initial grain/crystallite: number of grid points
  
  psi0_s=psi0_sol ! solid state (inside the grain)
  psi0_l=psi0_liq ! liquid state (outside the grain)
  
  i_psi0=0 ! reset n0 from grain lattice initialization, not using pre-fixed psi0
  !i_psi0=1 ! average density psi0 pre-fixed
  
  !noise_n=0.1d0 ! initial random noise inside nucleus
  noise_n=0
  noise=0 ! no initial noise outside nucleus
  noise0=noise
  
  r02=r0*r0*dx*dy
# endif
  
! discretization of time

  time0=0
  iter0=1
  t_i=1 ! initial transient time
!  t_i=0 ! no transient time
  dti=0.01d0 ! used for initial transient up to t=t_i

  dt=0.5d0 ! note: t_i/dt be integer

!  tmax=1000 ! maximum time
  tmax=10000
  
  n_i=NINT(t_i/dti)
  n_0=NINT(t_i/dt)-n_i ! for output
  nend=NINT((tmax-t_i)/dt)+n_i
  raout=NINT(100/dt)
  nout_conf=raout*50

! nmbr_eval: number of iterations used in the predictor-corrector algorithm
  nmbr_eval=1 ! no predictor-corrector iteration
!  nmbr_eval=2 ! one predictor-corrector iteration
!  nmbr_eval=100 ! maximum steps of iterations
  TOL=1.d-3
  err_s=TOL/10

  idum=-796921 ! random number seed; <0 for ran2()

  p_alloc_psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi, psi, [ly,lx])
  p_alloc_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiq, psiq, [ly/2+1,lx])
  
  allocate(q2(ly/2+1,lx))
  allocate(exp_1(ly/2+1,lx),cf_1(ly/2+1,lx),cf2_1(ly/2+1,lx))
  
! for FFTW parameters
  scale2d_b=1.d0/(ly*lx) ! for backward FFT
  do j=1,ly/2+1
    qy(j)= 2.d0*pi*(j-1)/(ly*dy)
  enddo 
  do i=1,lx
    if (i.le.lx/2+1) then
      qx(i)=2.d0*pi*(i-1)/(lx*dx)
    else
      qx(i)=-2.d0*pi*(lx-i+1)/(lx*dx)
    endif
  enddo

! initialization for FFTW (out of place; create the plans before initializing the inputs)
  plan2d_f=fftw_plan_dft_r2c_2d(lx,ly,psi,psiq,FFTW_MEASURE) ! note the reverse order of dimensions
  plan2d_b=fftw_plan_dft_c2r_2d(lx,ly,psiq,psi,FFTW_MEASURE)

! write parameters in a separate file
  open(unit=2,file=filename//'_para.dat')
  write(2,*) 'For rescaled 2D T-PFC equation'
# ifdef ini_homogeneous
  write(2,*) 'def ini_homogeneous: for homogeneous initial condition'
# endif
# ifdef ini_nucleus_crystal_cir
  write(2,*) 'def ini_nucleus_crystal_cir: initial crystalline nucleus (circular crystallite) in homogeneous medium'
  write(2,*) 'r0=',r0, '  noise_n(%)=',noise_n
  write(2,*) 'psi0_s=',psi0_s, '  psi0_l=',psi0_l
  write(2,*) 'i_psi0=',i_psi0
# endif
  write(2,*)
  write(2,*) 'lx=',lx, '  ly=',ly
  write(2,*) 'eps=',eps, '  g=',g
  write(2,*) 'alpha1=',alpha1, '  alpha2=',alpha2
  write(2,*) 'psi0_sol=',psi0_sol, '  psi0_liq=',psi0_liq
  write(2,*) 'q0=',q0, '  psi0=',psi0
  write(2,*) 'qx0=',qx0, '  qy0=',qy0
  write(2,*)
  write(2,*) 'dt=',dt, '  dti=',dti
  write(2,*) 'dx=',dx, '  dy=',dy
  write(2,*) 'n_dx=', n_dx
  write(2,*) 'tmax=',tmax, '  t_transient=',t_i
  write(2,*) 'nend=',nend, '  n_transient=',n_i, '  raout=',raout
  write(2,*) '  data output at t=t_transient, and every ',raout*dt
# ifdef save_conf
  write(2,*) '  save conf file every ',nout_conf*dt
# else
  write(2,*) '  ndef save_conf'
# endif
  write(2,*) 'time0=',time0, '  iter0=',iter0
  write(2,*) 'noise=',noise,'  noise0=',noise0
  write(2,*) '  random number seed: idum=',idum
  write(2,*)
  write(2,*) 'For the predictor-corrector method: '
  write(2,*) '    Maximum number of iterations:',nmbr_eval
  write(2,*) '    TOL=',TOL,'  err_s=',err_s
# ifdef open_MP 
  write(2,*) 'use openMP for multi-threading'
# else
  write(2,*) 'no openMP'
# endif
  close(2)

# ifdef open_MP
  ! use multi-threaded fftw with openMP
  nthreads=omp_get_max_threads()
  iret=fftw_init_threads()
  if(iret.eq.0) then
     write(*,*) 'iret=0: error during thread initialization'
     stop
  endif
  call fftw_plan_with_nthreads(nthreads)
  open(unit=2,file=filename//'_para.dat',position='append')
  write(2,*) 'use multi-threaded fftw; nthreads=',nthreads
  close(2)
# endif

  y=ran2(idum) ! initialize ran2
  
! set up the initial conditions
!
# ifdef ini_homogeneous
  ! homogeneous initial condition
  do j=1,ly
     do i=1,lx
        if(abs(psi0).lt.1.d-10) then
           psi(j,i)=psi0+noise0*(ran2(idum)-0.5d0) ! if psi0 ~ 0
        else
           !psi(j,i)=psi0*(1+noise*(ran2(idum)-0.5d0)) ! noise percentage
           psi(j,i)=psi0+noise*(ran2(idum)-0.5d0)
        endif
     enddo
  enddo
# endif

# ifdef ini_nucleus_crystal_cir
  ! initial condition of crystalline nucleus (circular crystallite) in homogeneous medium
  ! centered at (lx/2+0.5,ly/2+0.5)
  A0a=(g-3*psi0_s)**2-15*(-eps+psi0_s*(3*psi0_s-2*g))
  if(A0a > 0) then
     A0a=(g-3*psi0_s+sqrt(A0a))/15 ! one-mode solution
  else
     A0a=0.1d0
  endif
  do i=1,lx
     x=(i-lx/2-0.5d0)*dx
     do j=1,ly
        y=(j-ly/2-0.5d0)*dy
        if(x*x+y*y <= r02) then ! within the circular nucleus
           q_1=-qx0*x-qy0*y/2
           q_2=qy0*y
           q_3=qx0*x-qy0*y/2
           psi(j,i)=psi0_s+2*A0a*(cos(q_1)+cos(q_2)+cos(q_3)) ! reconstruct density field psi(y,x)
        else ! outside the nucleus
           psi(j,i)=psi0_l ! liquid/homogeneous state
        endif
     enddo
  enddo
  if(i_psi0==0) psi0=sum(psi)*scale2d_b ! =<psi>, reset average density of the whole system
  open(unit=2,file=filename//'_para.dat',position='append')
  write(2,*)
  write(2,*) 'psi0=',psi0
  write(2,*) 'initial crystalline nucleus amplitude A0a=',A0a
  close(2)
     
! add noise
  do i=1,lx
     x=(i-lx/2-0.5d0)*dx
     do j=1,ly
        y=(j-ly/2-0.5d0)*dy
        if(x*x+y*y <= r02) then ! within the circular nucleus (with noise_n)
           psi(j,i)=psi(j,i)*(1+noise_n*(ran2(idum)-0.5d0)) ! noise_n percentage
           !psi(j,i)=psi(j,i)+noise_n*(ran2(idum)-0.5d0)
        else ! outside the nucleus (with noise or noise0)
           if(abs(psi0_l).lt.1.d-10) then
              psi(j,i)=psi(j,i)+noise0*(ran2(idum)-0.5d0) ! if psi0_l ~ 0
           else
              psi(j,i)=psi(j,i)*(1+noise*(ran2(idum)-0.5d0)) ! noise percentage
              !psi(j,i)=psi(j,i)+noise*(ran2(idum)-0.5d0)
           endif
        endif
     enddo
  enddo
# endif
  
! set initial average value of psi
  call fftw_execute_dft_r2c(plan2d_f,psi,psiq) ! forward FFT
  psiq(1,1)=psi0*lx*ly ! <psi>=psi0
! backward FFT (note: the input array is overwritten for all c2r transforms)
  call fftw_execute_dft_c2r(plan2d_b,psiq,psi)
  psi=psi*scale2d_b

! get initial values of psiq
  call fftw_execute_dft_r2c(plan2d_f,psi,psiq) ! forward FFT

  time=time0
  iter=iter0-1
  if(iter0.le.n_i) then
     dtime=dti ! for initial transient
  else
     dtime=dt
  endif

  do i=1,lx
     do j=1,ly/2+1
        q2(j,i)=qx(i)*qx(i)+qy(j)*qy(j)
        alpha_1=-q2(j,i)*(-eps+(q2(j,i)-q02)**2) ! conserved dynamics
        alpha_dt=alpha_1*dtime
        exp_1(j,i)=exp(alpha_dt)
        if(abs(alpha_dt).lt.2.d-5) then
           cf_1(j,i)=dtime*(1.d0+0.5d0*alpha_dt*(1.d0+alpha_dt/3.d0))
           cf2_1(j,i)=0.5d0*dtime*(1.d0+alpha_dt*(1.d0+0.25d0*alpha_dt)/3.d0)
        else
           cf_1(j,i)=(exp_1(j,i)-1)/alpha_1
           cf2_1(j,i)=(exp_1(j,i)-(1+alpha_dt))/(alpha_1*alpha_dt)
        endif
     enddo
  enddo

  open(unit=9,file=filename//'_out.dat',position='append')
  write(9,*) 't=',time
  write(9,*) '<n(psi)>=',sum(psi)*scale2d_b,real(psiq(1,1))*scale2d_b
  write(9,*) 'n(psi)=',maxval(psi(1:ly,1:lx)),minval(psi(1:ly,1:lx))
  write(9,*)
  close(9)

! start running

  do iter=iter0,nend

    if(iter.le.n_i) then ! for initial transient only
       time=iter*dti
    else
       time=(iter-n_i)*dt+t_i ! NINT(time/dt)=iter+n_0
       if(iter.eq.n_i+1) then
          dtime=dt
          do i=1,lx
             do j=1,ly/2+1
                alpha_1=-q2(j,i)*(-eps+(q2(j,i)-q02)**2) ! conserved dynamics
                alpha_dt=alpha_1*dtime
                exp_1(j,i)=exp(alpha_dt)
                if(abs(alpha_dt).lt.2.d-5) then
                   cf_1(j,i)=dtime*(1.d0+0.5d0*alpha_dt*(1.d0+alpha_dt/3.d0))
                   cf2_1(j,i)=0.5d0*dtime*(1.d0+alpha_dt*(1.d0+0.25d0*alpha_dt)/3.d0)
                else
                   cf_1(j,i)=(exp_1(j,i)-1)/alpha_1
                   cf2_1(j,i)=(exp_1(j,i)-(1+alpha_dt))/(alpha_1*alpha_dt)
                endif
             enddo
          enddo
       endif
    endif

    call TPFC_eq(iter,istep)
    
    if((iter.eq.n_i).or.(iter.gt.n_i.and.mod(iter+n_0,raout).eq.0).or.(iter.eq.nend)) then

      open(unit=9,file=filename//'_out.dat',position='append')
      write(9,*) 't=',time,  ' istep=',istep
      write(9,*) '<n(psi)>=',sum(psi)*scale2d_b,real(psiq(1,1))*scale2d_b
      write(9,*) 'n(psi)=',maxval(psi(1:ly,1:lx)),minval(psi(1:ly,1:lx))
      write(9,*)
      close(9)

# ifdef save_conf
      ! save conf file
      if((iter.gt.n_i.and.mod(iter+n_0,nout_conf).eq.0).or.(iter.eq.nend)) then
         open(unit=11,file=filename//'_conf.dat',form='unformatted')
         write(11) psi
         close(11)
         open(unit=9,file=filename//'_out.dat',position='append')
         write(9,*) 't=',time, '  iter=',iter,' for ouput of conf file'
         write(9,*)
         close(9)
      endif
# endif

    endif ! for raout

  enddo ! for iter=iter0,nend

! output psi data file at t=tmax 
  open(unit=20,file=filename//'_psi.dat')
  do i=1,lx
     do j=1,ly
        write(20,*) psi(j,i) ! psi(j,i) data
     enddo
  enddo
  close(20)
  
  deallocate(q2,exp_1,cf_1,cf2_1)
  call fftw_free(p_alloc_psi)
  call fftw_free(p_alloc_psiq)

  call fftw_destroy_plan(plan2d_f)
  call fftw_destroy_plan(plan2d_b)

end Program Tpfc_2D


! main routine that solves the T-PFC equation
subroutine TPFC_eq(iter,istep)
  use global_variables
  implicit none

  integer :: i,j,iter,istep
  real(8) :: max_conv_psi,conv,abs_psi

  real(C_DOUBLE), pointer :: psi_0(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: psi0q(:,:),psiq_tmp(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Nonlin1_q(:,:),sigN_psiq(:,:)
! for allocating aligned memory in fftw (need different pointers for different arrays
! so that the memory can be released at the end of each call (otherwise the used memory will accumulate))
  type(C_PTR) :: p_alloc_psi0,p_alloc_psi0q,p_alloc_psiq_tmp,p_alloc_Nonlin1_q,p_alloc_sigN_psiq

  p_alloc_psi0=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi0, psi_0, [ly,lx])
  p_alloc_psi0q=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi0q, psi0q, [ly/2+1,lx])
  p_alloc_psiq_tmp=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiq_tmp, psiq_tmp, [ly/2+1,lx])
  p_alloc_Nonlin1_q=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin1_q, Nonlin1_q, [ly/2+1,lx])
  p_alloc_sigN_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_sigN_psiq, sigN_psiq, [ly/2+1,lx])
  
! calculate the nonlinear terms
  call nonlin_q(Nonlin1_q)

! First step: predictor values of psiq

  psi0q=exp_1*psiq+cf_1*Nonlin1_q
  if(nmbr_eval.gt.1) sigN_psiq=-cf2_1*Nonlin1_q ! for use in iteration

!$OMP PARALLEL DO
  do i=lx/2+2,lx
     psi0q(1,i)=conjg(psi0q(1,lx-i+2)) ! real data treatment for j=1 (qy=0) in 2D fftw, to avoid divergence at late time
  enddo
!$OMP END PARALLEL DO

  psiq_tmp=psi0q
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiq_tmp,psi_0)
  psi_0=psi_0*scale2d_b

  psi=psi_0 ! for nmbr_eval=1
  psiq=psi0q ! for nmbr_eval=1

! for use in iteration
  if(nmbr_eval.gt.1) psi0q=psi0q+sigN_psiq

! Second step: corrector
  do istep=2,nmbr_eval   ! get the corrector values of psiq by functional iteration
                         ! note: psi0q remain unchanged during the iteration

! for the nonlinear terms (note: here psi,psiq are the updated values)
     call nonlin_q(Nonlin1_q)

     psiq=psi0q+cf2_1*Nonlin1_q

!$OMP PARALLEL DO
     do i=lx/2+2,lx
        psiq(1,i)=conjg(psiq(1,lx-i+2)) ! real data treatment for j=1 (qy=0) in 2D fftw, to avoid divergence at late time
     enddo
!$OMP END PARALLEL DO

     psiq_tmp=psiq
! backward FFT (note: the input array is overwritten)
     call fftw_execute_dft_c2r(plan2d_b,psiq_tmp,psi)
     psi=psi*scale2d_b

! real space (psi) convergence
     if(nmbr_eval.gt.2) then

     max_conv_psi=0.d0
     do i=1,lx
        do j=1,ly
           abs_psi=abs(psi(j,i))
           if(abs_psi.gt.1.d5) then
              open(unit=2,file=filename//'_out.dat',position='append')
              write(2,*) 'psi diverges at t=',time, ' iter=',iter,' istep=',istep
              write(2,*) 'j=',j,' i=',i,' psi=',psi(j,i)
              close(2)
              stop
           endif
           if(abs_psi.le.err_s) then
              conv=abs(psi(j,i)-psi_0(j,i))
           else
              conv=abs((psi(j,i)-psi_0(j,i))/psi(j,i))
           endif
           if(conv.gt.max_conv_psi) max_conv_psi=conv
        enddo
     enddo

     if(max_conv_psi.lt.TOL) then
        call fftw_free(p_alloc_psi0)
        call fftw_free(p_alloc_psi0q)
        call fftw_free(p_alloc_psiq_tmp)
        call fftw_free(p_alloc_Nonlin1_q)
        call fftw_free(p_alloc_sigN_psiq)
        return
     endif

     psi_0=psi

     endif ! if nmbr_eval>2

  enddo ! for istep=2,nmbr_eval

  istep=istep-1
  if(nmbr_eval.gt.2) then
     open(unit=2,file=filename//'_out.dat',position='append')
     write(2,*)
     write(2,*) 'Exceeding maximum iterations=',nmbr_eval
     write(2,*) 'iter=',iter,'   t=',time,'   TOL=',TOL
     write(2,*) 'max_conv_psi=',max_conv_psi
     write(2,*) 'psi=',maxval(psi(1:ly,1:lx)),minval(psi(1:ly,1:lx))
     write(2,*)
     close(2)
  endif

  call fftw_free(p_alloc_psi0)
  call fftw_free(p_alloc_psi0q)
  call fftw_free(p_alloc_psiq_tmp)
  call fftw_free(p_alloc_Nonlin1_q)
  call fftw_free(p_alloc_sigN_psiq)
  return
end subroutine TPFC_eq


! calculate the nonlinear terms
subroutine nonlin_q(Nonlin1_q)
  use global_variables
  implicit none

  complex(C_DOUBLE_COMPLEX) :: Nonlin1_q(ly/2+1,lx)
  real(C_DOUBLE), pointer :: Nonlin1(:,:),dx_psi(:,:),dy_psi(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: dx_psiq(:,:),dy_psiq(:,:)
  type(C_PTR) :: p_alloc_Nonlin1 ! for allocating aligned memory in fftw
  type(C_PTR) :: p_alloc_dx_psi,p_alloc_dy_psi,p_alloc_dx_psiq,p_alloc_dy_psiq
  
  integer :: i,j

  p_alloc_Nonlin1=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin1, Nonlin1, [ly,lx])
  
  p_alloc_dx_psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dx_psi, dx_psi, [ly,lx])
  p_alloc_dy_psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dy_psi, dy_psi, [ly,lx])

  p_alloc_dx_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dx_psiq, dx_psiq, [ly/2+1,lx])
  p_alloc_dy_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dy_psiq, dy_psiq, [ly/2+1,lx])

  Nonlin1=psi*psi*(psi-g) ! original PFC nonlinear terms
  call fftw_execute_dft_r2c(plan2d_f,Nonlin1,Nonlin1_q) ! forward FFT
! for conserved dynamics
  Nonlin1_q=-q2*Nonlin1_q

! terms from transverse interaction

  if((i_alpha1.ne.0).or.(i_alpha2.ne.0)) then ! if alpha1 or alpha2 nonzero

  if(i_alpha2==0) then ! if alpha2=0
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        dy_psiq(j,i)=-alpha1*q2(j,i)*psiq(j,i) ! alpha1 \nabla^2 psi
        dx_psiq(j,i)=CMPLX(0.d0,qx(i))*dy_psiq(j,i) ! alpha1 \partial_x \nabla^2 psi
        dy_psiq(j,i)=CMPLX(0.d0,qy(j))*dy_psiq(j,i) ! alpha1 \partial_y \nabla^2 psi
     enddo
  enddo
!$OMP END PARALLEL DO
  else ! if i_alpha2=1, alpha2 nonzero
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        dy_psiq(j,i)=q2(j,i)*(alpha2*q2(j,i)-alpha1)*psiq(j,i) ! (alpha1 \nabla^2 + alpha2 \nabla^4) psi
        dx_psiq(j,i)=CMPLX(0.d0,qx(i))*dy_psiq(j,i) ! \partial_x (alpha1 \nabla^2 + alpha2 \nabla^4) psi
        dy_psiq(j,i)=CMPLX(0.d0,qy(j))*dy_psiq(j,i) ! \partial_y (alpha1 \nabla^2 + alpha2 \nabla^4) psi
     enddo
  enddo
!$OMP END PARALLEL DO
  endif
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,dx_psiq,dx_psi) ! = \partial_x (alpha1 \nabla^2 + alpha2 \nabla^4) psi
  dx_psi=dx_psi*scale2d_b
  call fftw_execute_dft_c2r(plan2d_b,dy_psiq,dy_psi) ! = \partial_y (alpha1 \nabla^2 + alpha2 \nabla^4) psi
  dy_psi=dy_psi*scale2d_b
  
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        dx_psiq(j,i)=CMPLX(0.d0,qx(i))*psiq(j,i) ! \partial_x psi
        dy_psiq(j,i)=CMPLX(0.d0,qy(j))*psiq(j,i) ! \partial_y psi
     enddo
  enddo
!$OMP END PARALLEL DO
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,dx_psiq,Nonlin1) ! = \partial_x psi
  Nonlin1=Nonlin1*dy_psi ! = (\partial_x psi)*[\partial_y (alpha1 \nabla^2 + alpha2 \nabla^4) psi]
  call fftw_execute_dft_c2r(plan2d_b,dy_psiq,dy_psi) ! = \partial_y psi
  Nonlin1=(Nonlin1-dy_psi*dx_psi)*scale2d_b ! = ...-(\partial_y psi)*[\partial_x (alpha1 \nabla^2 + alpha2 \nabla^4) psi]

  call fftw_execute_dft_r2c(plan2d_f,Nonlin1,dx_psiq) ! forward FFT
  Nonlin1_q=Nonlin1_q+dx_psiq

  endif ! if(i_alpha1.ne.0.or.i_alpha2.ne.0)

  call fftw_free(p_alloc_Nonlin1)
  call fftw_free(p_alloc_dx_psi)
  call fftw_free(p_alloc_dy_psi)
  call fftw_free(p_alloc_dx_psiq)
  call fftw_free(p_alloc_dy_psiq)
  
  return
end subroutine nonlin_q
