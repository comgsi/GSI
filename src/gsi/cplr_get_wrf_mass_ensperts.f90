module get_wrf_mass_ensperts_mod
use abstract_get_wrf_mass_ensperts_mod
  use kinds, only : i_kind
  type, extends(abstract_get_wrf_mass_ensperts_class) :: get_wrf_mass_ensperts_class
  contains
    procedure, pass(this) :: get_wrf_mass_ensperts => get_wrf_mass_ensperts_wrf
    procedure, pass(this) :: ens_spread_dualres_regional => ens_spread_dualres_regional_wrf
    procedure, pass(this) :: general_read_wrf_mass
    procedure, pass(this) :: parallel_read_wrf_mass_step1
    procedure, pass(this) :: parallel_read_wrf_mass_step2
    procedure, pass(this) :: general_read_wrf_mass2
! --- CAPS ---
    procedure, pass(this) :: get_fv3_mass_ensperts => get_fv3_mass_ensperts_fv3
    procedure, pass(this) :: general_read_fv3_mass
    procedure, pass(this) :: general_read_wrf_mass3
! --- CAPS ---
    procedure, nopass :: fill_regional_2d
  end type get_wrf_mass_ensperts_class
contains
  subroutine get_wrf_mass_ensperts_wrf(this,en_perts,nelen,ps_bar)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    get_wrf_mass_ensperts  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in netcdf format, for use
  !           with hybrid ensemble option.  ensemble spread is also written out as
  !           a byproduct for diagnostic purposes.
  !
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2011-08-31  todling - revisit en_perts (single-prec) in light of extended bundle
  !   2012-02-08  kleist  - add extra dimension to en_perts for 4d application 
  !   (currently use placeholder of value 1, since regional 4d application not 
  !   2017-07-30  Hu  - added code to read in multiple-time level ensemble forecast to
  !                     get 4D peerturbations
  !   2019-04-22  cTong - add CAPS radar DA option
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  
      use kinds, only: r_kind,i_kind,r_single
      use constants, only: zero,one,half,zero_single,rd_over_cp,one_tenth
      use mpimod, only: mpi_comm_world,ierror,mype,npe
      use hybrid_ensemble_parameters, only: n_ens,grd_ens,ens_fast_read
      use hybrid_ensemble_parameters, only: ntlevs_ens,ensemble_path
      use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
      use gsi_bundlemod, only: gsi_bundlecreate
      use gsi_bundlemod, only: gsi_grid
      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      use gsi_bundlemod, only: gsi_bundledestroy
      use gsi_bundlemod, only: gsi_gridcreate
      use mpeu_util, only: getindex
      use guess_grids,   only: ntguessig,ifilesig
      use gsi_4dvar,     only: nhr_assimilation
! --- CAPS --->
      use caps_radaruse_mod, only: l_use_log_qx, l_use_log_qx_pval, l_use_log_nt, cld_nt_updt ! chenll
      use caps_radaruse_mod, only: l_use_dbz_caps 
      use caps_radaruse_mod, only: l_cvpnr, cvpnr_pval                                            
! <--- CAPS ---
   

      implicit none
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      type(gsi_bundle),allocatable, intent(inout) :: en_perts(:,:)
      integer(i_kind), intent(in   ):: nelen
      real(r_single),dimension(:,:,:),allocatable:: ps_bar
  
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: u,v,tv,cwmr,oz,rh
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2):: ps
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)::w,qr,qi,qg,qs,qni,qnc,qnr
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)::dbz
  
      real(r_single),pointer,dimension(:,:,:):: w3
      real(r_single),pointer,dimension(:,:):: w2
      real(r_kind),pointer,dimension(:,:,:):: x3
      real(r_kind),pointer,dimension(:,:):: x2
      type(gsi_bundle):: en_bar 
      type(gsi_grid):: grid_ens
      real(r_kind):: bar_norm,sig_norm,kapr,kap1

      integer(i_kind):: i,j,k,n,mm1,istatus
      integer(i_kind):: ic2,ic3,i_radar_qr,i_radar_qg
      integer(i_kind):: its,ite, it

      character(255) filelists(ntlevs_ens)
      character(255) filename

      logical :: do_radar
      logical :: do_ens_fast_read
      
      ! Variables used only by the ensemble fast read
      integer(i_kind) :: iope
      logical :: bad_input
      real(r_kind),dimension(:,:,:),allocatable :: gg_u,gg_v,gg_tv,gg_rh
      real(r_kind),dimension(:,:),allocatable :: gg_ps

! --- CAPS --->
     !add by chenll
      integer,parameter   :: USEZG = 1
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig,n_ens):: qr_4d_tmp,qs_4d_tmp,qg_4d_tmp

      if ( l_use_dbz_caps ) then ! JP: Set flag to run this for CAPS only
         do n=1,n_ens,1
            do k=1,grd_ens%nsig
               do i=1,grd_ens%lon2
                  do j=1,grd_ens%lat2
                     qr_4d_tmp(j,i,k,n)=zero
                     qs_4d_tmp(j,i,k,n)=zero
                     qg_4d_tmp(j,i,k,n)=zero
                  end do
               end do
            end do
         end do
      end if
! <--- CAPS ----

      call gsi_gridcreate(grid_ens,grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)
      call gsi_bundlecreate(en_bar,grid_ens,'ensemble',istatus,names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
      if(istatus/=0) then
         write(6,*)' get_wrf_mass_ensperts_netcdf: trouble creating en_bar bundle'
         call stop2(999)
      endif
  
      if(ntlevs_ens > 1) then
         do i=1,ntlevs_ens
            write(filelists(i),'("filelist",i2.2)')ifilesig(i)
         enddo
         its=1
         ite=ntlevs_ens
      else
         write(filelists(1),'("filelist",i2.2)')nhr_assimilation
         its=ntguessig
         ite=ntguessig
      endif

      do it=its,ite
         if (mype == 0) write(*,*) 'ensemble file==',it,its,ite,ntlevs_ens,n_ens
         if(ntlevs_ens > 1) then
            open(10,file=trim(filelists(it)),form='formatted',err=30)
         else
            open(10,file=trim(filelists(1)),form='formatted',err=30)
         endif


  !
  ! INITIALIZE ENSEMBLE MEAN ACCUMULATORS
         en_bar%values=zero
  
         do n=1,n_ens
            en_perts(n,it)%valuesr4 = zero
         enddo
  
       if ( l_use_dbz_caps ) then ! CAPS
          mm1=mype+1
          kap1=rd_over_cp+one
          kapr=one/rd_over_cp
       else
  !    Determine if qr and qg are control variables for radar data assimilation,
         i_radar_qr=0
         i_radar_qg=0
         i_radar_qr=getindex(cvars3d,'qr')
         i_radar_qg=getindex(cvars3d,'qg')
         do_radar=i_radar_qr > 0 .and. i_radar_qg > 0

         mm1=mype+1
         kap1=rd_over_cp+one
         kapr=one/rd_over_cp

         ! If ens_fast_read is requested, check whether we really can use it.
         do_ens_fast_read = ens_fast_read
         can_ens_fast_read: if( do_ens_fast_read ) then ! make sure we can
         if(n_ens>npe) then
            do_ens_fast_read=.false.
130         format('Disabling ens_fast_read because number of ensemble members (',I0,') is greater than number of MPI ranks (',I0,').')
            if(mype==0) then
               write(6,130) n_ens,npe
            endif
         endif
         if(do_radar) then
            do_ens_fast_read=.false.
            if(mype==0) then
               write(6,'(A)') 'Disabling ens_fast_read because "radar mode" is in use (qg and qr are control variables).  Fast read is not yet implemented in "radar mode."'
            endif
         endif
         endif can_ens_fast_read
         if(do_ens_fast_read .and. mype==0) then
            write(6,'(I0,A)') mype,': will read ensemble data in parallel (ens_fast_read=.true.)'
         endif
      
  !
  ! If we're doing ens_fast_read, then this loop reads data.
         ens_parallel_read: if(do_ens_fast_read) then
         if(mype==0) then
            write(0,*) 'Will use ens_fast_read to read ARW ensemble.'
         endif
         ens_read_loop: do n=1,n_ens
            read(10,'(a)',err=20,end=20)filename
            filename=trim(ensemble_path) // trim(filename)
            iope=(n-1)*npe/n_ens
            if(mype==iope) then
               allocate(gg_u(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
               allocate(gg_v(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
               allocate(gg_tv(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
               allocate(gg_rh(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
               allocate(gg_ps(grd_ens%nlat,grd_ens%nlon))
               bad_input=.false.
               call this%parallel_read_wrf_mass_step1(filename,gg_ps,gg_tv,gg_u,gg_v,gg_rh)
            endif
         end do ens_read_loop
         
         call MPI_Barrier(mpi_comm_world,ierror)
         end if ens_parallel_read

         rewind(10)

       end if
  !
  ! LOOP OVER ENSEMBLE MEMBERS 
         ens_main_loop: do n=1,n_ens
  !
  ! DEFINE INPUT FILE NAME
  ! 
  ! READ OR SCATTER ENEMBLE MEMBER DATA
            read(10,'(a)',err=20,end=20)filename
            filename=trim(ensemble_path) // trim(filename)
            scatter_or_read: if(do_ens_fast_read) then
                ! Scatter data from the parallel read.
                iope=(n-1)*npe/n_ens
                if(mype==iope) then
                   write(0,'(I0,A,I0,A)') mype,': scatter member ',n,' to other ranks...'
                   call this%parallel_read_wrf_mass_step2(mype,iope,&
                        ps,u,v,tv,rh,cwmr,oz, &
                        gg_ps,gg_tv,gg_u,gg_v,gg_rh)
                else
                   call this%parallel_read_wrf_mass_step2(mype,iope,&
                        ps,u,v,tv,rh,cwmr,oz)
                endif
             else
                if (mype == 0) then
                   write(6,'(a,a)') 'CALL READ_WRF_MASS_ENSPERTS FOR ENS DATA : ',trim(filename)
                endif
                if ( l_use_dbz_caps) then ! CAPS
                  ! call this%general_read_wrf_mass3(filename,ps,u,v,tv,rh,cwmr,oz,qr,qs,qg,qnr,mype)
                   call this%general_read_wrf_mass3(filename,ps,u,v,tv,rh,cwmr,oz,qr,qs,qg,qnr,w,mype)
                else
                   if( do_radar )then
                      call this%general_read_wrf_mass2(filename,ps,u,v,tv,rh,cwmr,oz,w,dbz,qs,qg,qi,qr,qnc,qni,qnr,mype) 
                   else
                      call this%general_read_wrf_mass(filename,ps,u,v,tv,rh,cwmr,oz,mype) 
                   end if
                end if
             endif scatter_or_read

             call MPI_Barrier(mpi_comm_world,ierror)
  
  ! SAVE ENSEMBLE MEMBER DATA IN COLUMN VECTOR
            member_data_loop: do ic3=1,nc3d
  
               call gsi_bundlegetpointer(en_perts(n,it),trim(cvars3d(ic3)),w3,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' for ensemble member ',n
                  call stop2(999)
               end if
               call gsi_bundlegetpointer(en_bar,trim(cvars3d(ic3)),x3,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' for en_bar'
                  call stop2(999)
               end if
  
               select case (trim(cvars3d(ic3)))
  
                  case('sf','SF')
     
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = u(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+u(j,i,k)
                           end do
                        end do
                     end do
  
                  case('vp','VP')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = v(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+v(j,i,k)
                           end do
                        end do
                     end do
  
                  case('t','T')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = tv(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+tv(j,i,k)
                           end do
                        end do
                     end do
  
                  case('q','Q')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = rh(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+rh(j,i,k)
                           end do
                        end do
                     end do

               case('w','W')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = w(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+w(j,i,k)
                        end do
                     end do
                  end do

               case('qr','QR')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           if (l_use_log_qx) then ! CAPS
                               if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'log transform for qr : from member-->',n
                               if (qr(j,i,k) <= 1.0E-5_r_kind) then  !Originally Gang used 5.0E-5
                                  qr(j,i,k) = 1.0E-5_r_kind         !Rong Kong
                               end if
                               if (l_use_log_qx_pval .gt. 0.0_r_kind ) then ! CVpq
                                  qr(j,i,k) =((qr(j,i,k)**l_use_log_qx_pval)-1)/l_use_log_qx_pval  !chenll
                                  qr_4d_tmp(j,i,k,n)=qr(j,i,k) !chenll
                               else  ! CVlogq
                                  qr(j,i,k) = log(qr(j,i,k))
                               end if
                           end if
                           w3(j,i,k) = qr(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qr(j,i,k)
                        end do
                     end do
                  end do

               case('qs','QS')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           if (l_use_log_qx) then
                               if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'log transform for qs : from member-->',n
                               if (qs(j,i,k) <= 1.0E-5_r_kind) then
                                   qs(j,i,k) = 1.0E-5_r_kind
                               end if
                               if (l_use_log_qx_pval .gt. 0.0_r_kind ) then ! CVpq
                                  qs(j,i,k)=((qs(j,i,k)**l_use_log_qx_pval)-1)/l_use_log_qx_pval !chenll
                                  qs_4d_tmp(j,i,k,n)=qs(j,i,k) !chenll
                               else  ! CVlogq
                                  qs(j,i,k) = log(qs(j,i,k))
                               end if
                           end if
                           w3(j,i,k) = qs(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qs(j,i,k)
                        end do
                     end do
                  end do

               case('qi','QI')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qi(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qi(j,i,k)
                        end do
                     end do
                  end do

               case('qnr','QNR')
! --- CAPS ---                 
                  if ( l_use_dbz_caps) then
                    if ( cld_nt_updt .gt. 0 ) then
                       do k=1,grd_ens%nsig
                           do i=1,grd_ens%lon2
                               do j=1,grd_ens%lat2
                                   if (l_use_log_nt) then
                                       if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'updating qnr in hybrid analysis and log transform for qnr : from member-->',n
                                       if (qnr(j,i,k) < one) then
                                           qnr(j,i,k) = one
                                       end if
                                       qnr(j,i,k) = log(qnr(j,i,k))
                                   else
                                       if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'updating qnr in hybrid analysis : NO log-transform to member-->',n
                                   end if
                                   if (l_cvpnr) then
                                       if (qnr(j,i,k) < one) then
                                           qnr(j,i,k) = one
                                       end if
                                       qnr(j,i,k) = ((qnr(j,i,k)**cvpnr_pval)-1)/cvpnr_pval
                                   endif
                                   w3(j,i,k) = qnr(j,i,k)
                                   x3(j,i,k)=x3(j,i,k)+qnr(j,i,k)
                               end do
                           end do
                       end do
                    else
                       write(6,*)'qnr is not analyzed in hybrid analysis : member-->',n
                    end if
! --- CAPS ---
                  else
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = qnr(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+qnr(j,i,k)
                           end do
                        end do
                     end do
                  end if

               case('qnc','QNC')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qnc(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qnc(j,i,k)
                        end do
                     end do
                  end do

               case('qni','QNI')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = qni(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qni(j,i,k)
                        end do
                     end do
                  end do

               case('dbz','DBZ')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w3(j,i,k) = dbz(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+dbz(j,i,k)
                        end do
                     end do
                  end do

               case('qg','QG')

                  do k=1,grd_ens%nsig
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           if (l_use_log_qx) then ! CAPS
                               if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'log transform for qg : from member-->',n
                               if (qg(j,i,k) <= 1.0E-5_r_kind) then
                                   qg(j,i,k) = 1.0E-5_r_kind
                               end if
                               if (l_use_log_qx_pval .gt. 0.0_r_kind ) then ! CVpq
                                  qg(j,i,k)=((qg(j,i,k)**l_use_log_qx_pval)-1)/l_use_log_qx_pval !chenll
                                  qg_4d_tmp(j,i,k,n)=qg(j,i,k) !chenll
                               else  ! CVlogq
                                  qg(j,i,k) = log(qg(j,i,k))
                               end if
                           end if
                           w3(j,i,k) = qg(j,i,k)
                           x3(j,i,k)=x3(j,i,k)+qg(j,i,k)
                        end do
                     end do
                  end do
  
                  case('oz','OZ')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = oz(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+oz(j,i,k)
                           end do
                        end do
                     end do
  
                  case('cw','CW', 'ql', 'QL')
  
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = cwmr(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+cwmr(j,i,k)
                           end do
                        end do
                     end do
  
               end select
            end do member_data_loop
  
            member_mass_loop: do ic2=1,nc2d
     
               call gsi_bundlegetpointer(en_perts(n,it),trim(cvars2d(ic2)),w2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for ensemble member ',n
                  call stop2(999)
               end if
               call gsi_bundlegetpointer(en_bar,trim(cvars2d(ic2)),x2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for en_bar'
                  call stop2(999)
               end if
  
               select case (trim(cvars2d(ic2)))
  
                  case('ps','PS')
  
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = ps(j,i)
                           x2(j,i)=x2(j,i)+ps(j,i)
                        end do
                     end do
  
                  case('sst','SST')
  ! IGNORE SST IN HYBRID for now
  
                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = zero
                           x2(j,i)=zero
                        end do
                     end do
  
               end select
            end do member_mass_loop
         enddo ens_main_loop

! if ( 1 == 0) then  ! CAPS
!   call calc_num_rms_caps(qr_4d_tmp,qs_4d_tmp,qg_4d_tmp,l_use_log_qx_pval)
! end if

  !
  ! CALCULATE ENSEMBLE MEAN
         bar_norm = one/float(n_ens)
         en_bar%values=en_bar%values*bar_norm
  
  ! Copy pbar to module array.  ps_bar may be needed for vertical localization
  ! in terms of scale heights/normalized p/p
         pbar_loop: do ic2=1,nc2d
   
            if(trim(cvars2d(ic2)) == 'ps'.or.trim(cvars2d(ic2)) == 'PS') then
  
               call gsi_bundlegetpointer(en_bar,trim(cvars2d(ic2)),x2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for en_bar to get ps_bar'
                  call stop2(999)
               end if
   
               do i=1,grd_ens%lon2
                  do j=1,grd_ens%lat2
                     ps_bar(j,i,1)=x2(j,i)
                  end do
               end do
               exit
            end if
         end do pbar_loop
  
         call mpi_barrier(mpi_comm_world,ierror)
  !
  ! CALCULATE ENSEMBLE SPREAD
         call this%ens_spread_dualres_regional(mype,en_perts,nelen,en_bar)
         call mpi_barrier(mpi_comm_world,ierror)
  !
  ! CONVERT ENSEMBLE MEMBERS TO ENSEMBLE PERTURBATIONS
         sig_norm=sqrt(one/max(one,n_ens-one))
  
         do n=1,n_ens
            do i=1,nelen
               en_perts(n,it)%valuesr4(i)=(en_perts(n,it)%valuesr4(i)-en_bar%values(i))*sig_norm
            end do
         end do

     enddo ! it 4d loop
  !
     call gsi_bundledestroy(en_bar,istatus)
     if(istatus/=0) then
        write(6,*)' in get_wrf_mass_ensperts_netcdf: trouble destroying en_bar bundle'
               call stop2(999)
            endif

     if(allocated(gg_u)) deallocate(gg_u)
     if(allocated(gg_v)) deallocate(gg_v)
     if(allocated(gg_tv)) deallocate(gg_tv)
     if(allocated(gg_rh)) deallocate(gg_rh)
     if(allocated(gg_ps)) deallocate(gg_ps)

  return
30 write(6,*) 'get_wrf_mass_ensperts_netcdf: open filelist failed '
   call stop2(555)
20 write(6,*) 'get_wrf_mass_ensperts_netcdf: read WRF-ARW ens failed ',n
   call stop2(555)

  end subroutine get_wrf_mass_ensperts_wrf

  subroutine general_read_wrf_mass(this,filename,g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz,mype)
    use kinds, only: r_kind,i_kind,r_single
    use hybrid_ensemble_parameters, only: grd_ens
    implicit none
  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                    g_u,g_v,g_tv,g_rh,g_cwmr,g_oz
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2),intent(out):: g_ps
      character(255),intent(in):: filename
      integer,intent(in) :: mype

      real(r_kind),dimension(:,:,:),allocatable :: gg_u,gg_v,gg_tv,gg_rh
      real(r_kind),dimension(:,:),allocatable :: gg_ps

      if(mype==0) then
         allocate(gg_u(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
         allocate(gg_v(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
         allocate(gg_tv(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
         allocate(gg_rh(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
         allocate(gg_ps(grd_ens%nlat,grd_ens%nlon))
         call this%parallel_read_wrf_mass_step1(filename,gg_ps,gg_tv,gg_u,gg_v,gg_rh)
         call this%parallel_read_wrf_mass_step2(mype,0, &
              g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz, &
              gg_ps,gg_tv,gg_u,gg_v,gg_rh)
         deallocate(gg_u,gg_v,gg_tv,gg_rh,gg_ps)
      else
         call this%parallel_read_wrf_mass_step2(mype,0, &
              g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz)
      endif
  end subroutine general_read_wrf_mass
  
  subroutine parallel_read_wrf_mass_step1(this,filename,gg_ps,gg_tv,gg_u,gg_v,gg_rh)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    general_read_wrf_mass  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in "wrfout" netcdf format
  !           for use with hybrid ensemble option. 
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2010-09-10  parrish, modify so ensemble variables are read in the same way as in
  !               subroutines convert_netcdf_mass and read_wrf_mass_binary_guess.
  !               There were substantial differences due to different opinion about what
  !               to use for surface pressure.  This issue should be resolved by coordinating
  !               with Ming Hu (ming.hu@noaa.gov).  At the moment, these changes result in
  !               agreement to single precision between this input method and the guess input
  !               procedure when the same file is read by both methods.
  !   2012-03-12  whitaker:  read data on root, distribute with scatterv.
  !                          remove call to general_reload.
  !                          simplify, fix memory leaks, reduce memory footprint.
  !                          use genqsat, remove genqsat2_regional.
  !                          replace bare 'stop' statements with call stop2(999).
  !   2017-03-23  Hu      - add code to use hybrid vertical coodinate in WRF MASS core
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  
      use netcdf, only: nf90_nowrite
      use netcdf, only: nf90_open,nf90_close
      use netcdf, only: nf90_inq_dimid,nf90_inquire_dimension
      use netcdf, only: nf90_inq_varid,nf90_inquire_variable,nf90_get_var
      use kinds, only: r_kind,r_single,i_kind
      use gridmod, only: nsig,eta1_ll,pt_ll,aeta1_ll,eta2_ll,aeta2_ll
      use constants, only: zero,one,fv,zero_single,rd_over_cp_mass,one_tenth,h300
      use hybrid_ensemble_parameters, only: grd_ens,q_hyb_ens
      use netcdf_mod, only: nc_check
      use mpimod, only: mype

      implicit none
  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this

      character(255),intent(in):: filename

      real(r_kind),dimension(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig) :: gg_u,gg_v,gg_tv,gg_rh
      real(r_kind),dimension(grd_ens%nlat,grd_ens%nlon):: gg_ps

  !
  ! Declare local parameters
      real(r_kind),parameter:: r0_01 = 0.01_r_kind
      real(r_kind),parameter:: r10   = 10.0_r_kind
      real(r_kind),parameter:: r100  = 100.0_r_kind
  !
  !   Declare local variables
      real(r_single),allocatable,dimension(:):: temp_1d
      real(r_single),allocatable,dimension(:,:):: temp_2d,temp_2d2
      real(r_single),allocatable,dimension(:,:,:):: temp_3d
      real(r_kind),allocatable,dimension(:):: p_top
      real(r_kind),allocatable,dimension(:,:):: q_integral,q_integralc4h
      real(r_kind),allocatable,dimension(:,:,:):: tsn,qst,prsl
      integer(i_kind),allocatable,dimension(:):: dim,dim_id
  
      integer(i_kind):: nx,ny,nz,i,j,k,d_max,file_id,var_id,ndim
      integer(i_kind):: Time_id,s_n_id,w_e_id,b_t_id,s_n_stag_id,w_e_stag_id,b_t_stag_id
      integer(i_kind):: Time_len,s_n_len,w_e_len,b_t_len,s_n_stag_len,w_e_stag_len,b_t_stag_len
      integer(i_kind) iderivative
  
      real(r_kind):: deltasigma
      real(r_kind) psfc_this_dry,psfc_this
      real(r_kind) work_prslk,work_prsl
  
      logical ice

      character(len=24),parameter :: myname_ = 'general_read_wrf_mass'


      associate( this => this ) ! eliminates warning for unused dummy argument
                                ! needed for binding
      end associate

  !
  ! OPEN ENSEMBLE MEMBER DATA FILE
30 format(I0,': read ',A)
      write(6,*) mype,trim(filename)
      call nc_check( nf90_open(trim(filename),nf90_nowrite,file_id),&
          myname_,'open '//trim(filename) )
  !
  ! WRF FILE DIMENSIONS
      call nc_check( nf90_inq_dimid(file_id,'Time',Time_id),&
          myname_,'inq_dimid Time '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'south_north',s_n_id),&
          myname_,'inq_dimid south_north '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'west_east',w_e_id),&
          myname_,'inq_dimid west_east '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'bottom_top',b_t_id),&
          myname_,'inq_dimid bottom_top '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'south_north_stag',s_n_stag_id),&
          myname_,'inq_dimid south_north_stag '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'west_east_stag',w_e_stag_id),&
          myname_,'inq_dimid west_east_stag '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'bottom_top_stag',b_t_stag_id),&
          myname_,'inq_dimid bottom_top_stag '//trim(filename) )
  
      d_max=max(Time_id, s_n_id, w_e_id, b_t_id, s_n_stag_id, w_e_stag_id, b_t_stag_id)
      allocate(dim(d_max))
      dim(:)=-999
  
      call nc_check( nf90_inquire_dimension(file_id,Time_id,len=Time_len),&
          myname_,'inquire_dimension Time '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_id,len=s_n_len),&
          myname_,'inquire_dimension south_north '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_id,len=w_e_len),&
          myname_,'inquire_dimension west_east '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_id,len=b_t_len),&
          myname_,'inquire_dimension bottom_top '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_stag_id,len=s_n_stag_len),&
          myname_,'inquire_dimension south_north_stag '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_stag_id,len=w_e_stag_len),&
          myname_,'inquire_dimension west_east_stag '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_stag_id,len=b_t_stag_len),&
          myname_,'inquire_dimension bottom_top_stag '//trim(filename) )
  
      nx=w_e_len
      ny=s_n_len
      nz=b_t_len
      if (nx /= grd_ens%nlon .or. ny /= grd_ens%nlat .or. nz /= grd_ens%nsig) then
       print *,trim(filename)//': ','incorrect grid size in netcdf file'
       print *,trim(filename)//': ','nx,ny,nz,nlon,nlat,nsig',nx,ny,nz,grd_ens%nlon,grd_ens%nlat,grd_ens%nsig
       call stop2(999)
      endif
  
      dim(Time_id)=Time_len
      dim(s_n_id)=s_n_len
      dim(w_e_id)=w_e_len
      dim(b_t_id)=b_t_len
      dim(s_n_stag_id)=s_n_stag_len
      dim(w_e_stag_id)=w_e_stag_len
      dim(b_t_stag_id)=b_t_stag_len
  !
  ! READ PERTURBATION POTENTIAL TEMPERATURE (K)
  !    print *,trim(filename)//': ', 'read T ',filename
      call nc_check( nf90_inq_varid(file_id,'T',var_id),&
          myname_,'inq_varid T '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable T '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable T '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var T '//trim(filename) )
      allocate(tsn(dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))))
      tsn = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
  
  !  READ MU, MUB, P_TOP  (construct psfc as done in gsi--gives different result compared to PSFC)
  
      call nc_check( nf90_inq_varid(file_id,'P_TOP',var_id),&
          myname_,'inq_varid P_TOP '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable P_TOP '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable P_TOP '//trim(filename) )
      allocate(temp_1d(dim(dim_id(1))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_1d),&
          myname_,'get_var P_TOP '//trim(filename) )
      allocate(p_top(dim(dim_id(1))))
      do i=1,dim(dim_id(1))
         p_top(i)=temp_1d(i)
      enddo
      deallocate(dim_id)
  
      call nc_check( nf90_inq_varid(file_id,'MUB',var_id),&
          myname_,'inq_varid MUB '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable MUB '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable MUB '//trim(filename) )
      allocate(temp_2d(dim(dim_id(1)),dim(dim_id(2))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_2d),&
          myname_,'get_var MUB '//trim(filename) )
      deallocate(dim_id)
  
      call nc_check( nf90_inq_varid(file_id,'MU',var_id),&
          myname_,'inq_varid MU '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable MU '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable MU '//trim(filename) )
      allocate(temp_2d2(dim(dim_id(1)),dim(dim_id(2))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_2d2),&
          myname_,'get_var MU '//trim(filename) )
  
      do j=1,dim(dim_id(2))
         do i=1,dim(dim_id(1))
            temp_2d2(i,j)=temp_2d(i,j)+temp_2d2(i,j)+temp_1d(1)
            gg_ps(j,i)=temp_2d2(i,j)
         enddo
      enddo
      print *,trim(filename)//': ','min/max ps',minval(gg_ps),maxval(gg_ps)
      deallocate(temp_2d,temp_2d2,temp_1d,dim_id)
  
  !
  ! READ U (m/s)
      !print *,trim(filename)//': ', 'read U ',filename
      call nc_check( nf90_inq_varid(file_id,'U',var_id),&
          myname_,'inq_varid U '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable U '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable U '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var U '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))
         do j=1,dim(dim_id(2))
            do i=1,dim(dim_id(1))-1
               gg_u(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i+1,j,k))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
  !
  ! READ V (m/s)
      !print *,trim(filename)//': ', 'read V ',filename
      call nc_check( nf90_inq_varid(file_id,'V',var_id),&
          myname_,'inq_varid V '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable V '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable V '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var V '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))
         do j=1,dim(dim_id(2))-1
            do i=1,dim(dim_id(1))
               gg_v(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i,j+1,k))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,trim(filename)//': ','min/max u',minval(gg_u),maxval(gg_u)
      print *,trim(filename)//': ','min/max v',minval(gg_v),maxval(gg_v)
  !
  ! READ QVAPOR (kg/kg)
      !print *,trim(filename)//': ', 'read QVAPOR ',filename
      call nc_check( nf90_inq_varid(file_id,'QVAPOR',var_id),&
          myname_,'inq_varid QVAPOR '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QVAPOR '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QVAPOR '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QVAPOR '//trim(filename) )
      gg_rh = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id,dim)
  
      call nc_check( nf90_close(file_id),&
          myname_,'close '//trim(filename) )
  !
  ! CALCULATE TOTAL POTENTIAL TEMPERATURE (K)
      !print *,trim(filename)//': ', 'calculate total temperature ',filename
      do i=1,nx
         do j=1,ny
            do k=1,nz
              tsn(j,i,k)=tsn(j,i,k)+h300
            enddo
         enddo
      enddo   
  !
  ! INTEGRATE {1 + WATER VAPOR} TO CONVERT DRY AIR PRESSURE    
      !print *,trim(filename)//': ', 'integrate 1 + q vertically ',filename
      allocate(q_integral(ny,nx))
      allocate(q_integralc4h(ny,nx))
      q_integral(:,:)=one
      q_integralc4h=0.0_r_single
      do i=1,nx
         do j=1,ny
            do k=1,nz
               deltasigma=eta1_ll(k)-eta1_ll(k+1)
               q_integral(j,i)=q_integral(j,i)+deltasigma*gg_rh(j,i,k)
               q_integralc4h(j,i)=q_integralc4h(j,i)+(eta2_ll(k)-eta2_ll(k+1))*gg_rh(j,i,k)
            enddo
         enddo
      enddo
  !
  ! CONVERT WATER VAPOR MIXING RATIO TO SPECIFIC HUMIDITY
      do i=1,nx
         do j=1,ny
            do k=1,nz
               gg_rh(j,i,k)=gg_rh(j,i,k)/(one+gg_rh(j,i,k))
            enddo
         enddo
      enddo
  
  !  obtaining psfc as done in subroutine read_wrf_mass_netcdf_guess
      do i=1,nx
         do j=1,ny
            psfc_this_dry=r0_01*gg_ps(j,i)
            psfc_this=(psfc_this_dry-pt_ll)*q_integral(j,i)+pt_ll+q_integralc4h(j,i)
            gg_ps(j,i)=one_tenth*psfc_this  ! convert from mb to cb
         end do
      end do
  !
  ! CONVERT POTENTIAL TEMPERATURE TO VIRTUAL TEMPERATURE
      !print *,trim(filename)//': ', 'convert potential temp to virtual temp ',filename
      allocate(prsl(ny,nx,nz))
      do k=1,nz
         do i=1,nx
            do j=1,ny
               work_prsl  = one_tenth*(aeta1_ll(k)*(r10*gg_ps(j,i)-pt_ll)+&
                                       aeta2_ll(k) + pt_ll)
               prsl(j,i,k)=work_prsl
               work_prslk = (work_prsl/r100)**rd_over_cp_mass
               ! sensible temp from pot temp
               tsn(j,i,k)     = tsn(j,i,k)*work_prslk
               ! virtual temp from sensible temp
               gg_tv(j,i,k) = tsn(j,i,k) * (one+fv*gg_rh(j,i,k))
               ! recompute sensible temp from virtual temp
               tsn(j,i,k)= gg_tv(j,i,k)/(one+fv*max(zero,gg_rh(j,i,k)))
            end do
         end do
      end do
      print *,trim(filename)//': ','min/max tv',minval(gg_tv),maxval(gg_tv)
  
  !
  ! CALCULATE PSEUDO RELATIVE HUMIDITY IF USING RH VARIABLE
      if (.not.q_hyb_ens) then
         allocate(qst(ny,nx,nz))
         ice=.true. 
         iderivative=0
         call genqsat(qst,tsn,prsl,ny,nx,nsig,ice,iderivative)
         do k=1,nz
            do i=1,nx
               do j=1,ny
                  gg_rh(j,i,k)=gg_rh(j,i,k)/qst(j,i,k)
               enddo
            enddo
         enddo
         print *,trim(filename)//': ','min/max rh',minval(gg_rh),maxval(gg_rh)
         deallocate(qst)
      else
         print *,trim(filename)//': ','min/max q',minval(gg_rh),maxval(gg_rh)
      end if
  
  ! DEALLOCATE REMAINING TEMPORARY STORAGE
      deallocate(tsn,prsl,q_integral,p_top)
  
  return       
  end subroutine parallel_read_wrf_mass_step1

  subroutine parallel_read_wrf_mass_step2(this,mype,iope, &
       g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz, &
       gg_ps,gg_tv,gg_u,gg_v,gg_rh)

    use hybrid_ensemble_parameters, only: grd_ens
    use mpimod, only: mpi_comm_world,ierror,mpi_rtype
    use kinds, only: r_kind,r_single,i_kind
    implicit none

  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                    g_u,g_v,g_tv,g_rh,g_cwmr,g_oz
      integer(i_kind), intent(in) :: mype, iope
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2),intent(out):: g_ps

      ! The gg_ arrays are only sent by the rank doing I/O (mype==iope)
      real(r_kind),optional,dimension(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig) :: &
           gg_u,gg_v,gg_tv,gg_rh
      real(r_kind),optional,dimension(grd_ens%nlat,grd_ens%nlon):: gg_ps

  ! Declare local variables
      real(r_kind),allocatable,dimension(:):: wrk_send_2d
      integer(i_kind) :: k

  ! transfer data from root to subdomains on each task
  ! scatterv used, since full grids exist only on root task.
    allocate(wrk_send_2d(grd_ens%itotsub))
  ! first PS (output from fill_regional_2d is a column vector with a halo)
    if(mype==iope) call this%fill_regional_2d(gg_ps,wrk_send_2d)
    call mpi_scatterv(wrk_send_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
    g_ps,grd_ens%ijn_s(mype+1),mpi_rtype,iope,mpi_comm_world,ierror)       
  ! then TV,U,V,RH
    do k=1,grd_ens%nsig
       if (mype==iope) then
          call this%fill_regional_2d(gg_tv(:,:,k),wrk_send_2d)
       endif
       call mpi_scatterv(wrk_send_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_tv(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,iope,mpi_comm_world,ierror)       
       if (mype==iope) call this%fill_regional_2d(gg_u(1,1,k),wrk_send_2d)
       call mpi_scatterv(wrk_send_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_u(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,iope,mpi_comm_world,ierror)       
       if (mype==iope) call this%fill_regional_2d(gg_v(1,1,k),wrk_send_2d)
       call mpi_scatterv(wrk_send_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_v(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,iope,mpi_comm_world,ierror)       
       if (mype==iope) call this%fill_regional_2d(gg_rh(1,1,k),wrk_send_2d)
       call mpi_scatterv(wrk_send_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_rh(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,iope,mpi_comm_world,ierror)       
    enddo
  ! for now, don't do anything with oz, cwmr
    g_oz = 0.; g_cwmr = 0.
    deallocate(wrk_send_2d)
  end subroutine parallel_read_wrf_mass_step2  

  subroutine general_read_wrf_mass2(this,filename,g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz,&
                                    g_w,g_dbz,g_qs,g_qg,g_qi,g_qr,g_qnc,g_qni,g_qnr,mype)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    general_read_wrf_mass  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in "wrfout" netcdf format
  !           for use with hybrid ensemble option.
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2010-09-10  parrish, modify so ensemble variables are read in the same way
  !   as in
  !               subroutines convert_netcdf_mass and
  !               read_wrf_mass_binary_guess.
  !               There were substantial differences due to different opinion
  !               about what
  !               to use for surface pressure.  This issue should be resolved by
  !               coordinating
  !               with Ming Hu (ming.hu@noaa.gov).  At the moment, these changes
  !               result in
  !               agreement to single precision between this input method and
  !               the guess input
  !               procedure when the same file is read by both methods.
  !   2012-03-12  whitaker:  read data on root, distribute with scatterv.
  !                          remove call to general_reload.
  !                          simplify, fix memory leaks, reduce memory
  !                          footprint.
  !                          use genqsat, remove genqsat2_regional.
  !                          replace bare 'stop' statements with call
  !                          stop2(999).
  !   2017-03-23  Hu      - add code to use hybrid vertical coodinate in WRF
  !   MASS core
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block

      use netcdf, only: nf90_nowrite
      use netcdf, only: nf90_open,nf90_close
      use netcdf, only: nf90_inq_dimid,nf90_inquire_dimension
      use netcdf, only: nf90_inq_varid,nf90_inquire_variable,nf90_get_var
      use kinds, only: r_kind,r_single,i_kind
      use gridmod, only: nsig,eta1_ll,pt_ll,aeta1_ll,eta2_ll,aeta2_ll
      use constants, only: zero,one,fv,zero_single,rd_over_cp_mass,one_tenth,h300,rd,r1000
      use constants, only: r0_01,r10,r100
      use hybrid_ensemble_parameters, only: grd_ens,q_hyb_ens
      use mpimod, only: mpi_comm_world,ierror,mpi_rtype
      use netcdf_mod, only: nc_check
      use wrf_vars_mod, only : w_exist, dbz_exist
      use obsmod,only: if_model_dbz
      use setupdbz_lib,only: hx_dart

      implicit none
  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                    g_u,g_v,g_tv,g_rh,g_cwmr,g_oz, &
                                                    g_w,g_dbz,g_qs,g_qg,g_qi,g_qr, &
                                                    g_qnc,g_qni,g_qnr
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2),intent(out):: g_ps
      character(24),intent(in):: filename
  !
  !   Declare local variables
      real(r_single),allocatable,dimension(:):: temp_1d
      real(r_single),allocatable,dimension(:,:):: temp_2d,temp_2d2
      real(r_single),allocatable,dimension(:,:,:):: temp_3d
      real(r_kind),allocatable,dimension(:):: p_top
      real(r_kind),allocatable,dimension(:,:):: q_integral,gg_ps,q_integralc4h
      real(r_kind),allocatable,dimension(:,:,:):: tsn,qst,prsl,&
       gg_u,gg_v,gg_tv,gg_rh
      real(r_kind),allocatable,dimension(:,:,:):: gg_w,gg_qr,gg_qi,gg_qg,gg_qs,&
                                                  gg_dbz,gg_rho,gg_cwmr,gg_qnc,gg_qni,gg_qnr
      real(r_kind),allocatable,dimension(:):: wrk_fill_2d
      integer(i_kind),allocatable,dimension(:):: dim,dim_id

      integer(i_kind):: nx,ny,nz,i,j,k,d_max,file_id,var_id,ndim,mype
      integer(i_kind):: Time_id,s_n_id,w_e_id,b_t_id,s_n_stag_id,w_e_stag_id,b_t_stag_id
      integer(i_kind):: Time_len,s_n_len,w_e_len,b_t_len,s_n_stag_len,w_e_stag_len,b_t_stag_len
      integer(i_kind) iderivative

      real(r_kind):: deltasigma
      real(r_kind) psfc_this_dry,psfc_this
      real(r_kind) work_prslk,work_prsl

      logical ice

      character(len=24),parameter :: myname_ = 'general_read_wrf_mass2'


  !
  ! OPEN ENSEMBLE MEMBER DATA FILE
    if (mype==0) then ! only read data on root proc
      allocate(gg_u(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_v(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_tv(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_rh(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_ps(grd_ens%nlat,grd_ens%nlon))
      if( w_exist ) allocate(gg_w(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      if( dbz_exist ) allocate(gg_dbz(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qs(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qi(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qg(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_rho(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_cwmr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qnc(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qni(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qnr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      call nc_check( nf90_open(trim(filename),nf90_nowrite,file_id),&
          myname_,'open '//trim(filename) )
  !
  ! WRF FILE DIMENSIONS
      call nc_check( nf90_inq_dimid(file_id,'Time',Time_id),&
          myname_,'inq_dimid Time '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'south_north',s_n_id),&
          myname_,'inq_dimid south_north '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'west_east',w_e_id),&
          myname_,'inq_dimid west_east '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'bottom_top',b_t_id),&
          myname_,'inq_dimid bottom_top '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'south_north_stag',s_n_stag_id),&
          myname_,'inq_dimid south_north_stag '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'west_east_stag',w_e_stag_id),&
          myname_,'inq_dimid west_east_stag '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'bottom_top_stag',b_t_stag_id),&
          myname_,'inq_dimid bottom_top_stag '//trim(filename) )

      d_max=max(Time_id, s_n_id, w_e_id, b_t_id, s_n_stag_id, w_e_stag_id, b_t_stag_id)
      allocate(dim(d_max))
      dim(:)=-999
  
      call nc_check( nf90_inquire_dimension(file_id,Time_id,len=Time_len),&
          myname_,'inquire_dimension Time '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_id,len=s_n_len),&
          myname_,'inquire_dimension south_north '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_id,len=w_e_len),&
          myname_,'inquire_dimension west_east '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_id,len=b_t_len),&
          myname_,'inquire_dimension bottom_top '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_stag_id,len=s_n_stag_len),&
          myname_,'inquire_dimension south_north_stag '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_stag_id,len=w_e_stag_len),&
          myname_,'inquire_dimension west_east_stag '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_stag_id,len=b_t_stag_len),&
          myname_,'inquire_dimension bottom_top_stag '//trim(filename) )

      nx=w_e_len
      ny=s_n_len
      nz=b_t_len
      if (nx /= grd_ens%nlon .or. ny /= grd_ens%nlat .or. nz /= grd_ens%nsig) then
       print *,'incorrect grid size in netcdf file'
       print *,'nx,ny,nz,nlon,nlat,nsig',nx,ny,nz,grd_ens%nlon,grd_ens%nlat,grd_ens%nsig
       call stop2(999)
      endif

      dim(Time_id)=Time_len
      dim(s_n_id)=s_n_len
      dim(w_e_id)=w_e_len
      dim(b_t_id)=b_t_len
      dim(s_n_stag_id)=s_n_stag_len
      dim(w_e_stag_id)=w_e_stag_len
      dim(b_t_stag_id)=b_t_stag_len
  !
  ! READ PERTURBATION POTENTIAL TEMPERATURE (K)
      call nc_check( nf90_inq_varid(file_id,'T',var_id),&
          myname_,'inq_varid T '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable T '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable T '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var T '//trim(filename) )
      allocate(tsn(dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))))
      tsn = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)

  !  READ MU, MUB, P_TOP  (construct psfc as done in gsi--gives different result
  !  compared to PSFC)

      call nc_check( nf90_inq_varid(file_id,'P_TOP',var_id),&
          myname_,'inq_varid P_TOP '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable P_TOP '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable P_TOP '//trim(filename) )
      allocate(temp_1d(dim(dim_id(1))))

      call nc_check( nf90_get_var(file_id,var_id,temp_1d),&
          myname_,'get_var P_TOP '//trim(filename) )
      allocate(p_top(dim(dim_id(1))))
      do i=1,dim(dim_id(1))
         p_top(i)=temp_1d(i)
      enddo
      deallocate(dim_id)

      call nc_check( nf90_inq_varid(file_id,'MUB',var_id),&
          myname_,'inq_varid MUB '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable MUB '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable MUB '//trim(filename) )
      allocate(temp_2d(dim(dim_id(1)),dim(dim_id(2))))

      call nc_check( nf90_get_var(file_id,var_id,temp_2d),&
          myname_,'get_var MUB '//trim(filename) )
      deallocate(dim_id)

      call nc_check( nf90_inq_varid(file_id,'MU',var_id),&
          myname_,'inq_varid MU '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable MU '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable MU '//trim(filename) )
      allocate(temp_2d2(dim(dim_id(1)),dim(dim_id(2))))

      call nc_check( nf90_get_var(file_id,var_id,temp_2d2),&
          myname_,'get_var MU '//trim(filename) )

      do j=1,dim(dim_id(2))
         do i=1,dim(dim_id(1))
            temp_2d2(i,j)=temp_2d(i,j)+temp_2d2(i,j)+temp_1d(1)
            gg_ps(j,i)=temp_2d2(i,j)
         enddo
      enddo
      print *,'min/max ps',minval(gg_ps),maxval(gg_ps)
      deallocate(temp_2d,temp_2d2,temp_1d,dim_id)

  !
  ! READ U (m/s)
      call nc_check( nf90_inq_varid(file_id,'U',var_id),&
          myname_,'inq_varid U '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable U '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable U '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var U '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))
         do j=1,dim(dim_id(2))
            do i=1,dim(dim_id(1))-1
               gg_u(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i+1,j,k))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
  !
  ! READ V (m/s)
      call nc_check( nf90_inq_varid(file_id,'V',var_id),&
          myname_,'inq_varid V '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable V '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable V '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var V '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))
         do j=1,dim(dim_id(2))-1
            do i=1,dim(dim_id(1))
               gg_v(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i,j+1,k))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max u',minval(gg_u),maxval(gg_u)
      print *,'min/max v',minval(gg_v),maxval(gg_v)

  if( w_exist )then
  !
  ! READ W (m/s)
      call nc_check( nf90_inq_varid(file_id,'W',var_id),&
          myname_,'inq_varid W '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable W '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable W '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var W '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))-1
         do j=1,dim(dim_id(2))
            do i=1,dim(dim_id(1))
               gg_w(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i,j,k+1))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max w',minval(gg_w),maxval(gg_w)
  end if

  !
  ! READ QR (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'QRAIN',var_id),&
          myname_,'inq_varid QR '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QR '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QR '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QR '//trim(filename) )

      gg_qr = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qr',minval(gg_qr),maxval(gg_qr)

  !
  ! READ QS (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'QSNOW',var_id),&
          myname_,'inq_varid QS '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QS '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QS '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QS '//trim(filename) )

      gg_qs = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qs',minval(gg_qs),maxval(gg_qs)

  !
  ! READ QI (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'QICE',var_id),&
          myname_,'inq_varid QI '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QI '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QI '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QI '//trim(filename) )

      gg_qi = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qi',minval(gg_qi),maxval(gg_qi)

  !
  ! READ QG (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'QGRAUP',var_id),&
          myname_,'inq_varid QG '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QG '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QG '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QG '//trim(filename) )

      gg_qg = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qg',minval(gg_qg),maxval(gg_qg)

  !
  ! READ QNC
      call nc_check( nf90_inq_varid(file_id,'QNCLOUD',var_id),&
          myname_,'inq_varid QNC '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QNC '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QNC '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QNC '//trim(filename) )

      gg_qnc = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qnc',minval(gg_qnc),maxval(gg_qnc)

  !
  ! READ QNI
      call nc_check( nf90_inq_varid(file_id,'QNICE',var_id),&
          myname_,'inq_varid QNI '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QNI '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QNI '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QNI '//trim(filename) )

      gg_qni = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qni',minval(gg_qni),maxval(gg_qni)

  !
  ! READ QNR
      call nc_check( nf90_inq_varid(file_id,'QNRAIN',var_id),&
          myname_,'inq_varid QNR '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QNR '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QNR '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QNR '//trim(filename) )

      gg_qnr = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qnr',minval(gg_qnr),maxval(gg_qnr)

  !
  ! READ QC (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'QCLOUD',var_id),&
          myname_,'inq_varid QC '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QC '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QC '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QC '//trim(filename) )

      gg_cwmr = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max qc',minval(gg_cwmr),maxval(gg_cwmr)

  if( if_model_dbz .and. dbz_exist ) then
  !
  ! READ Reflectivity (dBZ)
      call nc_check( nf90_inq_varid(file_id,'REFL_10CM',var_id),&
          myname_,'inq_varid dBZ '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable dBZ '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable dBZ '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var dBZ '//trim(filename) )

      gg_dbz = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      where( gg_dbz < 0.0_r_kind )
        gg_dbz = 0.0_r_kind
      end where
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max dBZ',minval(gg_dbz),maxval(gg_dbz)
  end if

  !
  ! READ QVAPOR (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'QVAPOR',var_id),&
          myname_,'inq_varid QVAPOR '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QVAPOR '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QVAPOR '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QVAPOR '//trim(filename) )
      gg_rh = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id,dim)

      call nc_check( nf90_close(file_id),&
          myname_,'close '//trim(filename) )
  !
  ! CALCULATE TOTAL POTENTIAL TEMPERATURE (K)
      !print *, 'calculate total temperature ',filename
      do i=1,nx
         do j=1,ny
            do k=1,nz
              tsn(j,i,k)=tsn(j,i,k)+h300
            enddo
         enddo
      enddo
  !
  ! INTEGRATE {1 + WATER VAPOR} TO CONVERT DRY AIR PRESSURE
      allocate(q_integral(ny,nx))
      allocate(q_integralc4h(ny,nx))
      q_integral(:,:)=one
      q_integralc4h=0.0_r_single
      do i=1,nx
         do j=1,ny
            do k=1,nz
               deltasigma=eta1_ll(k)-eta1_ll(k+1)
               q_integral(j,i)=q_integral(j,i)+deltasigma*gg_rh(j,i,k)
               q_integralc4h(j,i)=q_integralc4h(j,i)+(eta2_ll(k)-eta2_ll(k+1))*gg_rh(j,i,k)
            enddo
         enddo
      enddo
  !
  ! CONVERT WATER VAPOR MIXING RATIO TO SPECIFIC HUMIDITY
      do i=1,nx
         do j=1,ny
            do k=1,nz
               gg_rh(j,i,k)=gg_rh(j,i,k)/(one+gg_rh(j,i,k))
            enddo
         enddo
      enddo

  !  obtaining psfc as done in subroutine read_wrf_mass_netcdf_guess
      do i=1,nx
         do j=1,ny
            psfc_this_dry=r0_01*gg_ps(j,i)
            psfc_this=(psfc_this_dry-pt_ll)*q_integral(j,i)+pt_ll+q_integralc4h(j,i)
            gg_ps(j,i)=one_tenth*psfc_this  ! convert from mb to cb
         end do
      end do
  !
  ! CONVERT POTENTIAL TEMPERATURE TO VIRTUAL TEMPERATURE
      allocate(prsl(ny,nx,nz))
      do k=1,nz
         do i=1,nx
            do j=1,ny
               work_prsl  = one_tenth*(aeta1_ll(k)*(r10*gg_ps(j,i)-pt_ll)+&
                                       aeta2_ll(k) + pt_ll)
               prsl(j,i,k)=work_prsl
               work_prslk = (work_prsl/r100)**rd_over_cp_mass
               ! sensible temp from pot temp
               tsn(j,i,k)     = tsn(j,i,k)*work_prslk
               ! virtual temp from sensible temp
               gg_tv(j,i,k) = tsn(j,i,k) * (one+fv*gg_rh(j,i,k))
               ! recompute sensible temp from virtual temp
               tsn(j,i,k)= gg_tv(j,i,k)/(one+fv*max(zero,gg_rh(j,i,k)))
            end do
         end do
      end do
      print *,'min/max tv',minval(gg_tv),maxval(gg_tv)

  if( dbz_exist .and. (.not. if_model_dbz) )then
     gg_rho = (prsl/(gg_tv*rd))*r1000  
      do k=1,nz
        do i=1,nx
          do j=1,ny
            call hx_dart(gg_qr(j,i,k),gg_qg(j,i,k),gg_qs(j,i,k),gg_rho(j,i,k),tsn(j,i,k),gg_dbz(j,i,k),.false.)
          enddo
        enddo
      enddo
  end if

  !
  ! CALCULATE PSEUDO RELATIVE HUMIDITY IF USING RH VARIABLE
      if (.not.q_hyb_ens) then
         allocate(qst(ny,nx,nz))
         ice=.true.
         iderivative=0
         call genqsat(qst,tsn,prsl,ny,nx,nsig,ice,iderivative)
         do k=1,nz
            do i=1,nx
               do j=1,ny
                  gg_rh(j,i,k)=gg_rh(j,i,k)/qst(j,i,k)
               enddo
            enddo
         enddo
         print *,'min/max rh',minval(gg_rh),maxval(gg_rh)
         deallocate(qst)
      else
         print *,'min/max q',minval(gg_rh),maxval(gg_rh)
      end if

  ! DEALLOCATE REMAINING TEMPORARY STORAGE
      deallocate(tsn,prsl,q_integral,p_top)
    endif ! done netcdf read on root

  ! transfer data from root to subdomains on each task
  ! scatterv used, since full grids exist only on root task.
    allocate(wrk_fill_2d(grd_ens%itotsub))
  ! first PS (output from fill_regional_2d is a column vector with a halo)
    if(mype==0) call this%fill_regional_2d(gg_ps,wrk_fill_2d)
    call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
    g_ps,grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
  ! then TV,U,V,RH
    do k=1,grd_ens%nsig
       if (mype==0) call this%fill_regional_2d(gg_tv(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_tv(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_u(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_u(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_v(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_v(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_rh(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_rh(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if(w_exist)then
         if (mype==0) call this%fill_regional_2d(gg_w(1,1,k),wrk_fill_2d)
         call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
         g_w(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       end if
       if(dbz_exist)then
         if (mype==0) call this%fill_regional_2d(gg_dbz(1,1,k),wrk_fill_2d)
         call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
         g_dbz(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       end if
       if (mype==0) call this%fill_regional_2d(gg_qr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qs(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qs(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qi(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qi(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qg(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qg(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_cwmr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_cwmr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qnc(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qnc(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qni(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qni(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call this%fill_regional_2d(gg_qnr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qnr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    enddo
  ! for now, don't do anything with oz, cwmr
    g_oz = 0.0_r_kind
    deallocate(wrk_fill_2d)
    if (mype==0) deallocate(gg_u,gg_v,gg_tv,gg_rh,gg_ps,gg_dbz,gg_w,&
                            gg_qr,gg_qs,gg_qi,gg_qg,gg_cwmr,gg_qnc, &
                            gg_qni,gg_qnr)

  return
  end subroutine general_read_wrf_mass2

! --- CAPS ---
  subroutine general_read_wrf_mass3 (this,filename,g_ps,g_u,g_v,g_tv,g_rh,g_cwmr,g_oz, &
                                                 g_qr,g_qs,g_qg,g_qnr,g_w,mype)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    general_read_wrf_mass3  read arw model ensemble members
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract: read ensemble members from the arw model in "wrfout" netcdf format
  !           for use with hybrid ensemble option. 
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2010-09-10  parrish, modify so ensemble variables are read in the same way as in
  !               subroutines convert_netcdf_mass and read_wrf_mass_binary_guess.
  !               There were substantial differences due to different opinion about what
  !               to use for surface pressure.  This issue should be resolved by coordinating
  !               with Ming Hu (ming.hu@noaa.gov).  At the moment, these changes result in
  !               agreement to single precision between this input method and the guess input
  !               procedure when the same file is read by both methods.
  !   2012-03-12  whitaker:  read data on root, distribute with scatterv.
  !                          remove call to general_reload.
  !                          simplify, fix memory leaks, reduce memory footprint.
  !                          use genqsat, remove genqsat2_regional.
  !                          replace bare 'stop' statements with call stop2(999).
  !   2017-03-23  Hu      - add code to use hybrid vertical coodinate in WRF MASS core
  !   2019-04-22  cTong   - add CAPS radar DA option
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  
      use netcdf, only: nf90_nowrite
      use netcdf, only: nf90_open,nf90_close
      use netcdf, only: nf90_inq_dimid,nf90_inquire_dimension
      use netcdf, only: nf90_inq_varid,nf90_inquire_variable,nf90_get_var
      use kinds, only: r_kind,r_single,i_kind
      use gridmod, only: nsig,eta1_ll,pt_ll,aeta1_ll,eta2_ll,aeta2_ll
      use constants, only: zero,one,fv,zero_single,rd_over_cp_mass,one_tenth,h300
      use hybrid_ensemble_parameters, only: grd_ens,q_hyb_ens
      use mpimod, only: mpi_comm_world,ierror,mpi_rtype
      use netcdf_mod, only: nc_check
      use caps_radaruse_mod, only: l_use_log_qx,l_use_log_qx_pval     !chenll

      implicit none
  !
  ! Declare passed variables
      class(get_wrf_mass_ensperts_class), intent(inout) :: this
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                    g_u,g_v,g_tv,g_rh,g_cwmr,g_oz, &
                                                    g_qr,g_qs,g_qg,g_qnr
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2),intent(out):: g_ps
      real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out)::g_w
      character(255),intent(in):: filename
  !
  ! Declare local parameters
      real(r_kind),parameter:: r0_01 = 0.01_r_kind
      real(r_kind),parameter:: r10   = 10.0_r_kind
      real(r_kind),parameter:: r100  = 100.0_r_kind
  !
  !   Declare local variables
      real(r_single),allocatable,dimension(:):: temp_1d
      real(r_single),allocatable,dimension(:,:):: temp_2d,temp_2d2
      real(r_single),allocatable,dimension(:,:,:):: temp_3d
      real(r_kind),allocatable,dimension(:):: p_top
      real(r_kind),allocatable,dimension(:,:):: q_integral,gg_ps,q_integralc4h
      real(r_kind),allocatable,dimension(:,:,:):: tsn,qst,prsl,&
       gg_u,gg_v,gg_tv,gg_rh
      real(r_kind),allocatable,dimension(:,:,:):: gg_qr,gg_qs,gg_qg
      real(r_kind),allocatable,dimension(:,:,:):: gg_qnr
      real(r_kind),allocatable,dimension(:,:,:):: gg_w
      real(r_kind),allocatable,dimension(:):: wrk_fill_2d
      integer(i_kind),allocatable,dimension(:):: dim,dim_id
  
      integer(i_kind):: nx,ny,nz,i,j,k,d_max,file_id,var_id,ndim,mype
      integer(i_kind):: Time_id,s_n_id,w_e_id,b_t_id,s_n_stag_id,w_e_stag_id,b_t_stag_id
      integer(i_kind):: Time_len,s_n_len,w_e_len,b_t_len,s_n_stag_len,w_e_stag_len,b_t_stag_len
      integer(i_kind) iderivative
  
      real(r_kind):: deltasigma
      real(r_kind) psfc_this_dry,psfc_this
      real(r_kind) work_prslk,work_prsl
  
      logical ice

      character(len=24),parameter :: myname_ = 'general_read_wrf_mass3'


  !
  ! OPEN ENSEMBLE MEMBER DATA FILE
    if (mype==0) then ! only read data on root proc
      allocate(gg_u(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_v(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_tv(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_rh(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qs(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qg(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_qnr(grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_w  (grd_ens%nlat,grd_ens%nlon,grd_ens%nsig))
      allocate(gg_ps(grd_ens%nlat,grd_ens%nlon))
      call nc_check( nf90_open(trim(filename),nf90_nowrite,file_id),&
          myname_,'open '//trim(filename) )
  !
  ! WRF FILE DIMENSIONS
      call nc_check( nf90_inq_dimid(file_id,'Time',Time_id),&
          myname_,'inq_dimid Time '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'south_north',s_n_id),&
          myname_,'inq_dimid south_north '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'west_east',w_e_id),&
          myname_,'inq_dimid west_east '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'bottom_top',b_t_id),&
          myname_,'inq_dimid bottom_top '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'south_north_stag',s_n_stag_id),&
          myname_,'inq_dimid south_north_stag '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'west_east_stag',w_e_stag_id),&
          myname_,'inq_dimid west_east_stag '//trim(filename) )
      call nc_check( nf90_inq_dimid(file_id,'bottom_top_stag',b_t_stag_id),&
          myname_,'inq_dimid bottom_top_stag '//trim(filename) )
  
      d_max=max(Time_id, s_n_id, w_e_id, b_t_id, s_n_stag_id, w_e_stag_id, b_t_stag_id)
      allocate(dim(d_max))
      dim(:)=-999
  
      call nc_check( nf90_inquire_dimension(file_id,Time_id,len=Time_len),&
          myname_,'inquire_dimension Time '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_id,len=s_n_len),&
          myname_,'inquire_dimension south_north '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_id,len=w_e_len),&
          myname_,'inquire_dimension west_east '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_id,len=b_t_len),&
          myname_,'inquire_dimension bottom_top '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_stag_id,len=s_n_stag_len),&
          myname_,'inquire_dimension south_north_stag '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_stag_id,len=w_e_stag_len),&
          myname_,'inquire_dimension west_east_stag '//trim(filename) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_stag_id,len=b_t_stag_len),&
          myname_,'inquire_dimension bottom_top_stag '//trim(filename) )
  
      nx=w_e_len
      ny=s_n_len
      nz=b_t_len
      if (nx /= grd_ens%nlon .or. ny /= grd_ens%nlat .or. nz /= grd_ens%nsig) then
       print *,'incorrect grid size in netcdf file'
       print *,'nx,ny,nz,nlon,nlat,nsig',nx,ny,nz,grd_ens%nlon,grd_ens%nlat,grd_ens%nsig
       call stop2(999)
      endif
  
      dim(Time_id)=Time_len
      dim(s_n_id)=s_n_len
      dim(w_e_id)=w_e_len
      dim(b_t_id)=b_t_len
      dim(s_n_stag_id)=s_n_stag_len
      dim(w_e_stag_id)=w_e_stag_len
      dim(b_t_stag_id)=b_t_stag_len
  !
  ! READ PERTURBATION POTENTIAL TEMPERATURE (K)
  !    print *, 'read T ',filename
      call nc_check( nf90_inq_varid(file_id,'T',var_id),&
          myname_,'inq_varid T '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable T '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable T '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var T '//trim(filename) )
      allocate(tsn(dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))))
      tsn = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
  
  !  READ MU, MUB, P_TOP  (construct psfc as done in gsi--gives different result compared to PSFC)
  
      call nc_check( nf90_inq_varid(file_id,'P_TOP',var_id),&
          myname_,'inq_varid P_TOP '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable P_TOP '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable P_TOP '//trim(filename) )
      allocate(temp_1d(dim(dim_id(1))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_1d),&
          myname_,'get_var P_TOP '//trim(filename) )
      allocate(p_top(dim(dim_id(1))))
      do i=1,dim(dim_id(1))
         p_top(i)=temp_1d(i)
      enddo
      deallocate(dim_id)
  
      call nc_check( nf90_inq_varid(file_id,'MUB',var_id),&
          myname_,'inq_varid MUB '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable MUB '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable MUB '//trim(filename) )
      allocate(temp_2d(dim(dim_id(1)),dim(dim_id(2))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_2d),&
          myname_,'get_var MUB '//trim(filename) )
      deallocate(dim_id)
  
      call nc_check( nf90_inq_varid(file_id,'MU',var_id),&
          myname_,'inq_varid MU '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable MU '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable MU '//trim(filename) )
      allocate(temp_2d2(dim(dim_id(1)),dim(dim_id(2))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_2d2),&
          myname_,'get_var MU '//trim(filename) )
  
      do j=1,dim(dim_id(2))
         do i=1,dim(dim_id(1))
            temp_2d2(i,j)=temp_2d(i,j)+temp_2d2(i,j)+temp_1d(1)
            gg_ps(j,i)=temp_2d2(i,j)
         enddo
      enddo
      print *,'min/max ps',minval(gg_ps),maxval(gg_ps)
      deallocate(temp_2d,temp_2d2,temp_1d,dim_id)
  
  !
  ! READ U (m/s)
      !print *, 'read U ',filename
      call nc_check( nf90_inq_varid(file_id,'U',var_id),&
          myname_,'inq_varid U '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable U '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable U '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var U '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))
         do j=1,dim(dim_id(2))
            do i=1,dim(dim_id(1))-1
               gg_u(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i+1,j,k))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
  !
  ! READ V (m/s)
      !print *, 'read V ',filename
      call nc_check( nf90_inq_varid(file_id,'V',var_id),&
          myname_,'inq_varid V '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable V '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable V '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var V '//trim(filename) )
  !
  ! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))
         do j=1,dim(dim_id(2))-1
            do i=1,dim(dim_id(1))
               gg_v(j,i,k)=.5*(temp_3d(i,j,k)+temp_3d(i,j+1,k))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id)
      print *,'min/max u',minval(gg_u),maxval(gg_u)
      print *,'min/max v',minval(gg_v),maxval(gg_v)
  !
  ! READ QVAPOR (kg/kg)
      !print *, 'read QVAPOR ',filename
      call nc_check( nf90_inq_varid(file_id,'QVAPOR',var_id),&
          myname_,'inq_varid QVAPOR '//trim(filename) )
  
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QVAPOR '//trim(filename) )
      allocate(dim_id(ndim))
  
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QVAPOR '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
  
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QVAPOR '//trim(filename) )
      gg_rh = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      deallocate(temp_3d)
      deallocate(dim_id)
  
  !
  ! CALCULATE TOTAL POTENTIAL TEMPERATURE (K)
      !print *, 'calculate total temperature ',filename
      do i=1,nx
         do j=1,ny
            do k=1,nz
              tsn(j,i,k)=tsn(j,i,k)+h300
            enddo
         enddo
      enddo   
  !
  ! INTEGRATE {1 + WATER VAPOR} TO CONVERT DRY AIR PRESSURE    
      !print *, 'integrate 1 + q vertically ',filename
      allocate(q_integral(ny,nx))
      allocate(q_integralc4h(ny,nx))
      q_integral(:,:)=one
      q_integralc4h=0.0_r_single
      do i=1,nx
         do j=1,ny
            do k=1,nz
               deltasigma=eta1_ll(k)-eta1_ll(k+1)
               q_integral(j,i)=q_integral(j,i)+deltasigma*gg_rh(j,i,k)
               q_integralc4h(j,i)=q_integralc4h(j,i)+(eta2_ll(k)-eta2_ll(k+1))*gg_rh(j,i,k)
            enddo
         enddo
      enddo
  !
  ! CONVERT WATER VAPOR MIXING RATIO TO SPECIFIC HUMIDITY
      do i=1,nx
         do j=1,ny
            do k=1,nz
               gg_rh(j,i,k)=gg_rh(j,i,k)/(one+gg_rh(j,i,k))
            enddo
         enddo
      enddo
  
  !  obtaining psfc as done in subroutine read_wrf_mass_netcdf_guess
      do i=1,nx
         do j=1,ny
            psfc_this_dry=r0_01*gg_ps(j,i)
            psfc_this=(psfc_this_dry-pt_ll)*q_integral(j,i)+pt_ll+q_integralc4h(j,i)
            gg_ps(j,i)=one_tenth*psfc_this  ! convert from mb to cb
         end do
      end do
  !
  ! CONVERT POTENTIAL TEMPERATURE TO VIRTUAL TEMPERATURE
      !print *, 'convert potential temp to virtual temp ',filename
      allocate(prsl(ny,nx,nz))
      do k=1,nz
         do i=1,nx
            do j=1,ny
               work_prsl  = one_tenth*(aeta1_ll(k)*(r10*gg_ps(j,i)-pt_ll)+&
                                       aeta2_ll(k) + pt_ll)
               prsl(j,i,k)=work_prsl
               work_prslk = (work_prsl/r100)**rd_over_cp_mass
               ! sensible temp from pot temp
               tsn(j,i,k)     = tsn(j,i,k)*work_prslk
               ! virtual temp from sensible temp
               gg_tv(j,i,k) = tsn(j,i,k) * (one+fv*gg_rh(j,i,k))
               ! recompute sensible temp from virtual temp
               tsn(j,i,k)= gg_tv(j,i,k)/(one+fv*max(zero,gg_rh(j,i,k)))
            end do
         end do
      end do
      print *,'min/max tv',minval(gg_tv),maxval(gg_tv)
  
  !
  ! CALCULATE PSEUDO RELATIVE HUMIDITY IF USING RH VARIABLE
      if (.not.q_hyb_ens) then
         allocate(qst(ny,nx,nz))
         ice=.true. 
         iderivative=0
         call genqsat(qst,tsn,prsl,ny,nx,nsig,ice,iderivative)
         do k=1,nz
            do i=1,nx
               do j=1,ny
                  gg_rh(j,i,k)=gg_rh(j,i,k)/qst(j,i,k)
               enddo
            enddo
         enddo
         print *,'min/max rh',minval(gg_rh),maxval(gg_rh)
         deallocate(qst)
      else
         print *,'min/max q',minval(gg_rh),maxval(gg_rh)
      end if
  
!------------------------------------------------------------------------------!
!      Read ensemble cloud hydrometer mixing ratios
!
! READ QRAIN (kg/kg)
      print *, 'read QRAIN ',filename
      call nc_check( nf90_inq_varid(file_id,'QRAIN',var_id),&
          myname_,'inq_varid QRAIN '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QRAIN '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QRAIN '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
    
      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QRAIN '//trim(filename) )
      gg_qr = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      print *,'min/max QRIAN',minval(gg_qr),maxval(gg_qr)
      deallocate(temp_3d)
      deallocate(dim_id)

! READ QSNOW (kg/kg)
      !print *, 'read QSNOW ',filename
      call nc_check( nf90_inq_varid(file_id,'QSNOW',var_id),&
          myname_,'inq_varid QSNOW '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QSNOW '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QSNOW '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QSNOW '//trim(filename) )
      gg_qs = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      print *,'min/max QSNOW',minval(gg_qs),maxval(gg_qs)
      deallocate(temp_3d)
      deallocate(dim_id)

! READ QGRAUPEL (kg/kg)
      !print *, 'read QGRAUP ',filename
      call nc_check( nf90_inq_varid(file_id,'QGRAUP',var_id),&
          myname_,'inq_varid QGRAUP '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QGRAUP '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QGRAUP '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QGRAUP '//trim(filename) )
      gg_qg = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      print *,'min/max QGRAUP',minval(gg_qg),maxval(gg_qg)
      deallocate(temp_3d)
      deallocate(dim_id)

! READ QNRAIN (#/kg)
      !print *, 'read QNRAIN ',filename
      call nc_check( nf90_inq_varid(file_id,'QNRAIN',var_id),&
          myname_,'inq_varid QNRAIN '//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable QNRAIN '//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable QNRAIN '//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var QNRAIN '//trim(filename) )
      gg_qnr = reshape(temp_3d,(/dim(dim_id(2)),dim(dim_id(1)),dim(dim_id(3))/),order=(/2,1,3/))
      print *,'min/max QNRAIN',minval(gg_qnr),maxval(gg_qnr)
      deallocate(temp_3d)
      deallocate(dim_id)

! READ w (m/s)
      print *, 'read W ',filename
      call nc_check( nf90_inq_varid(file_id,'W',var_id),&
          myname_,'inq_varid W'//trim(filename) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable W'//trim(filename) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable W'//trim(filename) )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var W'//trim(filename) )
!
! INTERPOLATE TO MASS GRID
      do k=1,dim(dim_id(3))-1
         do j=1,dim(dim_id(2))
            do i=1,dim(dim_id(1))
                gg_w(j,i,k)=0.5*(temp_3d(i,j,k)+temp_3d(i,j,k+1))
            enddo
         enddo
      enddo
      deallocate(temp_3d)
      deallocate(dim_id,dim)
      print *,'min/max W',minval(gg_w),maxval(gg_w)

      call nc_check( nf90_close(file_id),&
          myname_,'close '//trim(filename) )

  ! DEALLOCATE REMAINING TEMPORARY STORAGE
      deallocate(tsn,prsl,q_integral,p_top)
    endif ! done netcdf read on root
  
  ! transfer data from root to subdomains on each task
  ! scatterv used, since full grids exist only on root task.
    allocate(wrk_fill_2d(grd_ens%itotsub))
  ! first PS (output from fill_regional_2d is a column vector with a halo)
    if(mype==0) call this%fill_regional_2d(gg_ps,wrk_fill_2d)
    call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
    g_ps,grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
  ! then TV,U,V,RH
  !               ,and qr,qs,qg, qnr, and w
    do k=1,grd_ens%nsig
       if (mype==0) call this%fill_regional_2d(gg_tv(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_tv(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
       if (mype==0) call this%fill_regional_2d(gg_u(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_u(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
       if (mype==0) call this%fill_regional_2d(gg_v(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_v(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)       
       if (mype==0) call this%fill_regional_2d(gg_rh(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_rh(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
 
       if (mype==0) call fill_regional_2d(gg_qr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call fill_regional_2d(gg_qs(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qs(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call fill_regional_2d(gg_qg(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qg(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
       if (mype==0) call fill_regional_2d(gg_qnr(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_qnr(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)

       if (mype==0) call fill_regional_2d(gg_w(1,1,k),wrk_fill_2d)
       call mpi_scatterv(wrk_fill_2d,grd_ens%ijn_s,grd_ens%displs_s,mpi_rtype, &
       g_w(1,1,k),grd_ens%ijn_s(mype+1),mpi_rtype,0,mpi_comm_world,ierror)

    enddo
  ! for now, don't do anything with oz, cwmr
    g_oz = 0.; g_cwmr = 0.
    deallocate(wrk_fill_2d)
    if (mype==0) deallocate(gg_u,gg_v,gg_tv,gg_rh,gg_ps)
    if (mype==0) deallocate(gg_qr,gg_qs,gg_qg)
    if (mype==0) deallocate(gg_qnr)
    if (mype==0) deallocate(gg_w)
  
  return       
  end subroutine general_read_wrf_mass3
! --- CAPS ---

  subroutine fill_regional_2d(fld_in,fld_out)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    fill_regional_2d
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract:  create a column vector for the subdomain (including halo)
  ! from global 2d grid.
  !
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2012-03-12  whitaker, remove nx,ny,itotsub from argument list.
  !
  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
    use kinds, only: r_kind,i_kind
    use hybrid_ensemble_parameters, only: grd_ens
    implicit none
    real(r_kind),dimension(grd_ens%nlat,grd_ens%nlon)::fld_in
    real(r_kind),dimension(grd_ens%itotsub)::fld_out
    integer(i_kind):: i,j,k
    do k=1,grd_ens%itotsub
       i=grd_ens%ltosj_s(k)
       j=grd_ens%ltosi_s(k)
       fld_out(k)=fld_in(j,i)
    enddo
  return 
  end subroutine fill_regional_2d

  subroutine ens_spread_dualres_regional_wrf(this,mype,en_perts,nelen,en_bar)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    ens_spread_dualres_regional
  !   prgmmr: mizzi            org: ncar/mmm            date: 2010-08-11
  !
  ! abstract:
  !
  !
  ! program history log:
  !   2010-08-11  parrish, initial documentation
  !   2011-04-05  parrish - add pseudo-bundle capability
  !   2011-08-31  todling - revisit en_perts (single-prec) in light of extended bundle
  !   2019-04-22  cTong   - add CAPS radar DA option
  !
  !   input argument list:
  !     en_bar - ensemble mean
  !      mype  - current processor number
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
  !
    use kinds, only: r_single,r_kind,i_kind
    use hybrid_ensemble_parameters, only: n_ens,grd_ens,grd_anl,p_e2a,uv_hyb_ens, &
                                          regional_ensemble_option,write_ens_sprd
    use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info,general_sube2suba
    use constants, only:  zero,two,half,one
    use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
    use gsi_bundlemod, only: gsi_bundlecreate
    use gsi_bundlemod, only: gsi_grid
    use gsi_bundlemod, only: gsi_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer
    use gsi_bundlemod, only: gsi_bundledestroy
    use gsi_bundlemod, only: gsi_gridcreate
    implicit none

    class(get_wrf_mass_ensperts_class), intent(inout) :: this
    type(gsi_bundle),OPTIONAL,intent(in):: en_bar
    integer(i_kind),intent(in):: mype
    type(gsi_bundle),allocatable, intent(in   ) :: en_perts(:,:)
    integer(i_kind), intent(in   ):: nelen
  
    type(gsi_bundle):: sube,suba
    type(gsi_grid):: grid_ens,grid_anl
    real(r_kind) sp_norm,sig_norm_sq_inv
    type(sub2grid_info)::se,sa
    integer(i_kind) k
  
    integer(i_kind) i,n,ic3
    logical regional
    integer(i_kind) num_fields,inner_vars,istat,istatus
    logical,allocatable::vector(:)
    real(r_kind),pointer,dimension(:,:,:):: st,vp,tv,rh,oz,cw
    real(r_kind),pointer,dimension(:,:):: ps
! --- CAPS ---
    real(r_kind),pointer,dimension(:,:,:):: qr,qs,qg
    real(r_kind),pointer,dimension(:,:,:):: qnr
! --- CAPS ---
    real(r_kind),dimension(grd_anl%lat2,grd_anl%lon2,grd_anl%nsig),target::dum3
    real(r_kind),dimension(grd_anl%lat2,grd_anl%lon2),target::dum2

    associate( this => this ) ! eliminates warning for unused dummy argument needed for binding
    end associate
 
  !      create simple regular grid
          call gsi_gridcreate(grid_anl,grd_anl%lat2,grd_anl%lon2,grd_anl%nsig)
          call gsi_gridcreate(grid_ens,grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)
  
  !      create two internal bundles, one on analysis grid and one on ensemble grid
  
         call gsi_bundlecreate (suba,grid_anl,'ensemble work',istatus, &
                                   names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
         if(istatus/=0) then
            write(6,*)' in ens_spread_dualres_regional: trouble creating bundle_anl bundle'
            call stop2(999)
         endif
         call gsi_bundlecreate (sube,grid_ens,'ensemble work ens',istatus, &
                                   names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
         if(istatus/=0) then
            write(6,*)' ens_spread_dualres_regional: trouble creating bundle_ens bundle'
            call stop2(999)
         endif
  
    sp_norm=(one/float(n_ens))
  
    sube%values=zero
  !
  
    if(regional_ensemble_option == 1)then
       print *,'global ensemble'
       sig_norm_sq_inv=n_ens-one
  
       do n=1,n_ens
          do i=1,nelen
             sube%values(i)=sube%values(i) &
               +en_perts(n,1)%valuesr4(i)*en_perts(n,1)%valuesr4(i)
          end do
       end do
  
       do i=1,nelen
         sube%values(i) = sqrt(sp_norm*sig_norm_sq_inv*sube%values(i))
       end do
    else
       do n=1,n_ens
          do i=1,nelen
             sube%values(i)=sube%values(i) &
               +(en_perts(n,1)%valuesr4(i)-en_bar%values(i))*(en_perts(n,1)%valuesr4(i)-en_bar%values(i))
          end do
       end do
   
       do i=1,nelen
         sube%values(i) = sqrt(sp_norm*sube%values(i))
       end do
    end if
  
    if(grd_ens%latlon1n == grd_anl%latlon1n) then
       do i=1,nelen
          suba%values(i)=sube%values(i)
       end do
    else
       inner_vars=1
       num_fields=max(0,nc3d)*grd_ens%nsig+max(0,nc2d)
       allocate(vector(num_fields))
       vector=.false.
       do ic3=1,nc3d
          if(trim(cvars3d(ic3))=='sf'.or.trim(cvars3d(ic3))=='vp') then
             do k=1,grd_ens%nsig
                vector((ic3-1)*grd_ens%nsig+k)=uv_hyb_ens
             end do
          end if
       end do
       call general_sub2grid_create_info(se,inner_vars,grd_ens%nlat,grd_ens%nlon,grd_ens%nsig,num_fields, &
                                         regional,vector)
       call general_sub2grid_create_info(sa,inner_vars,grd_anl%nlat,grd_anl%nlon,grd_anl%nsig,num_fields, &
                                         regional,vector)
       deallocate(vector)
       call general_sube2suba(se,sa,p_e2a,sube%values,suba%values,regional)
    end if
  
    dum2=zero
    dum3=zero
    call gsi_bundlegetpointer(suba,'sf',st,istat)
    if(istat/=0) then
       write(6,*)' no sf pointer in ens_spread_dualres, point st at dum3 array'
       st => dum3
    end if
    call gsi_bundlegetpointer(suba,'vp',vp,istat)
    if(istat/=0) then
       write(6,*)' no vp pointer in ens_spread_dualres, point vp at dum3 array'
       vp => dum3
    end if
    call gsi_bundlegetpointer(suba,'t',tv,istat)
    if(istat/=0) then
       write(6,*)' no t pointer in ens_spread_dualres, point tv at dum3 array'
       tv => dum3
    end if
    call gsi_bundlegetpointer(suba,'q',rh,istat)
    if(istat/=0) then
       write(6,*)' no q pointer in ens_spread_dualres, point rh at dum3 array'
       rh => dum3
    end if
    call gsi_bundlegetpointer(suba,'oz',oz,istat)
    if(istat/=0) then
       write(6,*)' no oz pointer in ens_spread_dualres, point oz at dum3 array'
       oz => dum3
    end if
! --- CAPS ___
! cloud hydrometer mixing ratio
    call gsi_bundlegetpointer(suba,'qr',qr,istat)
    if(istat/=0) then
       write(6,*)' no qr pointer in ens_spread_dualres, point qr at dum3 array'
       qr => dum3
    end if
    call gsi_bundlegetpointer(suba,'qs',qs,istat)
    if(istat/=0) then
       write(6,*)' no qs pointer in ens_spread_dualres, point qs at dum3 array'
       qs => dum3
    end if
    call gsi_bundlegetpointer(suba,'qg',qg,istat)
    if(istat/=0) then
       write(6,*)' no qg pointer in ens_spread_dualres, point qg at dum3 array'
       qg => dum3
    end if
    call gsi_bundlegetpointer(suba,'qnr',qnr,istat)
    if(istat/=0) then
       write(6,*)' no qnr pointer in ens_spread_dualres, point qnr at dum3 array'
       qnr => dum3
    end if
! -- CAPS ---
    call gsi_bundlegetpointer(suba,'cw',cw,istat)
    if(istat/=0) then
       write(6,*)' no cw pointer in ens_spread_dualres, point cw at dum3 array'
       cw => dum3
    end if
    call gsi_bundlegetpointer(suba,'ps',ps,istat)
    if(istat/=0) then
       write(6,*)' no ps pointer in ens_spread_dualres, point ps at dum2 array'
       ps => dum2
    end if
  
    if(write_ens_sprd) call write_spread_dualres(st,vp,tv,rh,oz,cw,ps,mype)
    if(write_ens_sprd) call write_spread_dualres_qcld_regional(qr,qs,qg,qnr,mype) ! CAPS
  
    return
  end subroutine ens_spread_dualres_regional_wrf

subroutine write_spread_dualres_qcld_regional(a,b,c,d,mype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    write_spread_dualres_qcld_regional   write ensemble spread for diagnostics
!   prgmmr: cTong           org: CAPS/OU              date: 2019-04-22
!
! abstract: write ensemble spread (previously interpolated to analysis grid)
!             for hydrometeor variables for diagnostic purposes.
!
!
! program history log:
!  
!   input argument list:
!     a    -  spread variable 1
!     b    -  spread variable 2
!     c    -  spread variable 3
!     d    -  spread variable 4
!     mype -  current processor number
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$ end documentation block
  use kinds, only: r_kind,i_kind,r_single
  use hybrid_ensemble_parameters, only: grd_anl
  use constants, only: zero
  implicit none

  integer(i_kind),intent(in):: mype
  character(255):: grdfile

  real(r_kind),dimension(grd_anl%lat2,grd_anl%lon2,grd_anl%nsig),intent(in):: a,b,c,d
  real(r_kind),dimension(grd_anl%lat2,grd_anl%lon2,grd_anl%nsig,4):: g3in

  real(r_kind),dimension(grd_anl%nlat,grd_anl%nlon,grd_anl%nsig):: work8_3d

  real(r_single),dimension(grd_anl%nlon,grd_anl%nlat,grd_anl%nsig):: work4_3d

  integer(i_kind) ncfggg,iret,i,j,k,n,mem2d,mem3d,num3d

! Initial memory used by 2d and 3d grids
  mem2d = 4*grd_anl%nlat*grd_anl%nlon
  mem3d = 4*grd_anl%nlat*grd_anl%nlon*grd_anl%nsig
  num3d=4

! transfer 2d arrays to generic work aray
  do k=1,grd_anl%nsig
    do j=1,grd_anl%lon2
       do i=1,grd_anl%lat2
         g3in(i,j,k,1)=a(i,j,k)
         g3in(i,j,k,2)=b(i,j,k)
         g3in(i,j,k,3)=c(i,j,k)
         g3in(i,j,k,4)=d(i,j,k)
       end do
     end do
  end do

  if (mype==0) then
    grdfile='ens_spread_qcld_reg.grd'
    ncfggg=len_trim(grdfile)
    call baopenwt(22,grdfile(1:ncfggg),iret)
    write(6,*)'WRITE_SPREAD_DUALRES_QCLD_REGIONAL:  open 22 to ',trim(grdfile),' with iret=',iret
  endif

! Process 3d arrays
  do n=1,num3d
    work8_3d=zero
    do k=1,grd_anl%nsig
      call gather_stuff2(g3in(1,1,k,n),work8_3d(1,1,k),mype,0)
    end do
    if (mype==0) then
      do k=1,grd_anl%nsig
        do j=1,grd_anl%nlon
          do i=1,grd_anl%nlat
            work4_3d(j,i,k) =work8_3d(i,j,k)
          end do
        end do
      end do
      call wryte(22,mem3d,work4_3d)
      write(6,*)'WRITE_SPREAD_DUALRES_QCLD_REGIONAL FOR VARIABLE NUM ',n,' min/max:',minval(work4_3d),maxval(work4_3d)
    endif
  end do

! Close byte-addressable binary file for grads
  if (mype==0) then
     call baclose(22,iret)
     write(6,*)'WRITE_SPREAD_DUALRES_QCLD_REGIONAL:  close 22 with iret=',iret
  end if

  return
end subroutine write_spread_dualres_qcld_regional

  subroutine get_fv3_mass_ensperts_fv3(this,en_perts,nelen,ps_bar)

!$$$  subprogram documentation block
!
! subprogram:   get_fv3_ensperts
!   prgmmr: cTong           org: CAPS/OU                date: 2019-04-23
!
! abstract: read ensemble members, and construct ensemble perturbations, for use
!             with hybrid ensemble option.
!
! program history log:
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!
!$$$ end documentation block
   use kinds, only: r_kind,i_kind,r_single
   use constants,only: zero,one,one_tenth,ten,fv
   use mpimod, only: mpi_comm_world,ierror,mype
   use gridmod,  only: eta1_ll
   use hybrid_ensemble_parameters, only: n_ens,grd_ens,q_hyb_ens
   use hybrid_ensemble_parameters, only: ensemble_path
   use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
   use gsi_bundlemod, only: gsi_bundlecreate
   use gsi_bundlemod, only: gsi_grid
   use gsi_bundlemod, only: gsi_bundle
   use gsi_bundlemod, only: gsi_bundlegetpointer
   use gsi_bundlemod, only: gsi_bundledestroy
   use gsi_bundlemod, only: gsi_gridcreate
   use guess_grids,   only: ntguessig
   use caps_radaruse_mod, only: l_use_log_qx, l_use_log_qx_pval, l_use_dbz_caps, l_use_log_nt, cld_nt_updt
   use caps_radaruse_mod, only: l_cvpnr, cvpnr_pval

   implicit none
   class(get_wrf_mass_ensperts_class), intent(inout) :: this
   type(gsi_bundle),allocatable, intent(inout) :: en_perts(:,:)
   integer(i_kind), intent(in   ):: nelen
   real(r_single),dimension(:,:,:),allocatable:: ps_bar

   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: u,v,ts,tv,oz,rh,qst
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: qr,qs,qg
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: qnr
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: w
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: ql,qi
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: prsl

   real(r_single),pointer,dimension(:,:,:):: w3
   real(r_single),pointer,dimension(:,:):: w2
   real(r_kind),pointer,dimension(:,:,:):: x3
   real(r_kind),pointer,dimension(:,:):: x2
   type(gsi_bundle):: en_bar
   type(gsi_grid):: grid_ens
   real(r_kind):: bar_norm,sig_norm,kapr,kap1

   integer(i_kind):: i,j,k,n,istatus
   integer(i_kind):: ic2,ic3
   integer(i_kind):: its,ite, it

   character(255) filelists(2)
   character(255) dynmfile
   character(255) tracfile

! --- CAPS --->
     !add by chenll
   integer,parameter   :: USEZG = 1
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig,n_ens):: qr_4d_tmp,qs_4d_tmp,qg_4d_tmp

   if ( l_use_dbz_caps ) then ! JP: Set flag to run this for CAPS only
      do n=1,n_ens,1
         do k=1,grd_ens%nsig
            do i=1,grd_ens%lon2
               do j=1,grd_ens%lat2
                  qr_4d_tmp(j,i,k,n)=zero
                  qs_4d_tmp(j,i,k,n)=zero
                  qg_4d_tmp(j,i,k,n)=zero
               end do
            end do
         end do
      end do
   end if
! <--- CAPS ----




   call gsi_gridcreate(grid_ens,grd_ens%lat2,grd_ens%lon2,grd_ens%nsig)
   call gsi_bundlecreate(en_bar,grid_ens,'ensemble',istatus,names2d=cvars2d,names3d=cvars3d,bundle_kind=r_kind)
   if(istatus/=0) then
      write(6,*)' get_fv3_mass_ensperts_fv3: trouble creating en_bar bundle'
      call stop2(999)
   endif

   write(filelists(1),'("filelist02")')
   write(filelists(2),'("filelist00")')
   its=ntguessig
   ite=ntguessig

   do it=its,ite
      if (mype == 0) write(*,*) 'number of ensemble file=',n_ens
      open(10,file=trim(filelists(1)),form='formatted',err=30)
      open(20,file=trim(filelists(2)),form='formatted',err=40)

  !
  ! INITIALIZE ENSEMBLE MEAN ACCUMULATORS
      en_bar%values=zero

      do n=1,n_ens
         en_perts(n,it)%valuesr4 = zero
      enddo

  !
  ! LOOP OVER ENSEMBLE MEMBERS
      do n=1,n_ens
  !
  ! DEFINE INPUT FILE NAME
          read(10,'(a)',err=20,end=20)dynmfile
          read(20,'(a)',err=20,end=20)tracfile
          dynmfile=trim(ensemble_path) // trim(dynmfile)
          tracfile=trim(ensemble_path) // trim(tracfile)
  !
  ! READ ENEMBLE MEMBERS DATA
          if (mype == 0) write(6,'(a,a)') 'CALL READ_FV3_MASS_ENSPERTS FOR ENS DATA (dynamic): ',trim(dynmfile)
          if (mype == 0) write(6,'(a,a)') 'CALL READ_FV3_MASS_ENSPERTS FOR ENS DATA ( tracer): ',trim(tracfile)
!          call this%general_read_fv3_mass(dynmfile,tracfile,prsl,u,v,tv,rh,oz,qr,qs,qg,qnr,w,mype) 
          call this%general_read_fv3_mass(dynmfile,tracfile,prsl,u,v,tv,rh,oz,ql,qi,qr,qs,qg,qnr,w,mype)  

  ! SAVE ENSEMBLE MEMBER DATA IN COLUMN VECTOR
            do ic3=1,nc3d

               call gsi_bundlegetpointer(en_perts(n,it),trim(cvars3d(ic3)),w3,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' for ensemble member ',n
                  call stop2(999)
               end if
               call gsi_bundlegetpointer(en_bar,trim(cvars3d(ic3)),x3,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars3d(ic3)),' for en_bar'
                  call stop2(999)
               end if

               select case (trim(cvars3d(ic3)))

                  case('sf','SF')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = u(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+u(j,i,k)
                           end do
                        end do
                     end do

                  case('vp','VP')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = v(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+v(j,i,k)
                           end do
                        end do
                     end do

                  case('t','T')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = tv(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+tv(j,i,k)
                           end do
                        end do
                     end do

                  case('q','Q')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = rh(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+rh(j,i,k)
                           end do
                        end do
                     end do

                  case('oz','OZ')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = oz(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+oz(j,i,k)
                           end do
                        end do
                     end do          

                  case('ql','QL')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = ql(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+ql(j,i,k)
                           end do
                        end do
                     end do

                  case('qi','QI')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = qi(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+qi(j,i,k)
                           end do
                        end do
                     end do

                  case('qr','QR')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              if (l_use_log_qx) then
                                  if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'log transform for qr : from member-->',n
                                  if (qr(j,i,k) <= 1.0E-5_r_kind) then !Originally Gang used 5.0E-5
                                      qr(j,i,k) = 1.0E-5_r_kind         !Rong Kong
                                  end if
                                  if (l_use_log_qx_pval .gt. 0.0_r_kind ) then ! CVpq
                                     qr(j,i,k) =((qr(j,i,k)**l_use_log_qx_pval)-1)/l_use_log_qx_pval  !chenll
                                     qr_4d_tmp(j,i,k,n)=qr(j,i,k) !chenll
                                  else  ! CVlogq
                                     qr(j,i,k) = log(qr(j,i,k))
                                  end if
                              end if
                              w3(j,i,k) = qr(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+qr(j,i,k)
                           end do
                        end do
                     end do

                  case('qs','QS')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              if (l_use_log_qx) then
                                  if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'log transform for qs : from member-->',n
                                  if (qs(j,i,k) <= 1.0E-5_r_kind) then 
                                      qs(j,i,k) = 1.0E-5_r_kind 
                                  end if
                                  if (l_use_log_qx_pval .gt. 0.0_r_kind ) then ! CVpq
                                     qs(j,i,k)=((qs(j,i,k)**l_use_log_qx_pval)-1)/l_use_log_qx_pval !chenll
                                     qs_4d_tmp(j,i,k,n)=qs(j,i,k) !chenll
                                  else  ! CVlogq
                                     qs(j,i,k) = log(qs(j,i,k))
                                  end if
                              end if
                              w3(j,i,k) = qs(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+qs(j,i,k)
                           end do
                        end do
                     end do

                  case('qg','QG')

                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              if (l_use_log_qx) then
                                  if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'log transform for qg : from member-->',n
                                  if (qg(j,i,k) <= 1.0E-5_r_kind) then
                                     qg(j,i,k) = 1.0E-5_r_kind
                                  end if
                                  if (l_use_log_qx_pval .gt. 0.0_r_kind ) then ! CVpq
                                     qg(j,i,k)=((qg(j,i,k)**l_use_log_qx_pval)-1)/l_use_log_qx_pval !chenll
                                     qg_4d_tmp(j,i,k,n)=qg(j,i,k) !chenll
                                  else  ! CVlogq
                                     qg(j,i,k) = log(qg(j,i,k))
                                  end if
                              end if
                              w3(j,i,k) = qg(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+qg(j,i,k)
                           end do
                        end do
                     end do

                  case('qnr','QNR')

                     if ( cld_nt_updt .gt. 0 ) then
                       do k=1,grd_ens%nsig
                          do i=1,grd_ens%lon2
                             do j=1,grd_ens%lat2
                                 if (l_use_log_nt) then
                                   if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'updating qnr in hybrid analysis and log transform for qnr : from member-->',n
                                     if (qnr(j,i,k) < one) then
                                       qnr(j,i,k) = one
                                     end if
                                     qnr(j,i,k) = log(qnr(j,i,k))
                                 else
                                   if(mype==0 .and. i==10 .and. j==10 .and. k==10) write(6,*)'updating qnr in hybrid analysis : NO log-transform to member-->',n
                                 end if
                                   if (l_cvpnr) then
                                       if (qnr(j,i,k) < one) then
                                           qnr(j,i,k) = one
                                       end if
                                       qnr(j,i,k) =((qnr(j,i,k)**cvpnr_pval)-1)/cvpnr_pval
                                   endif

                                 w3(j,i,k) = qnr(j,i,k)
                                 x3(j,i,k)=x3(j,i,k)+qnr(j,i,k)
                             end do
                          end do
                       end do
                     else
                       write(6,*)'qnr is not analyzed in hybrid analysis : member-->',n
                     end if

                  case('w','W')
                     do k=1,grd_ens%nsig
                        do i=1,grd_ens%lon2
                           do j=1,grd_ens%lat2
                              w3(j,i,k) = w(j,i,k)
                              x3(j,i,k)=x3(j,i,k)+w(j,i,k)
                           end do
                        end do
                     end do 

               end select
            end do

            do ic2=1,nc2d

               call gsi_bundlegetpointer(en_perts(n,it),trim(cvars2d(ic2)),w2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for ensemble member ',n
                  call stop2(999)
               end if
               call gsi_bundlegetpointer(en_bar,trim(cvars2d(ic2)),x2,istatus)
               if(istatus/=0) then
                  write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for en_bar'
                  call stop2(999)
               end if

               select case (trim(cvars2d(ic2)))

                  case('ps','PS')

                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = prsl(j,i,1)
                           x2(j,i)=x2(j,i)+prsl(j,i,1)
                        end do
                     end do

                  case('sst','SST')
  ! IGNORE SST IN HYBRID for now

                     do i=1,grd_ens%lon2
                        do j=1,grd_ens%lat2
                           w2(j,i) = zero
                           x2(j,i)=zero
                        end do
                     end do

               end select
            end do
      enddo ! n ensemble mem loop
  !

!     if ( 1 == 0) then  ! CAPS
!     call calc_num_rms_caps(qr_4d_tmp,qs_4d_tmp,qg_4d_tmp,l_use_log_qx_pval)
!     end if

  ! CALCULATE ENSEMBLE MEAN
      bar_norm = one/float(n_ens)
      en_bar%values=en_bar%values*bar_norm

  ! Copy pbar to module array.  ps_bar may be needed for vertical localization
  ! in terms of scale heights/normalized p/p
      do ic2=1,nc2d

         if(trim(cvars2d(ic2)) == 'ps'.or.trim(cvars2d(ic2)) == 'PS') then

            call gsi_bundlegetpointer(en_bar,trim(cvars2d(ic2)),x2,istatus)
            if(istatus/=0) then
               write(6,*)' error retrieving pointer to ',trim(cvars2d(ic2)),' for en_bar to get ps_bar'
               call stop2(999)
            end if

            do i=1,grd_ens%lon2
               do j=1,grd_ens%lat2
                  ps_bar(j,i,1)=x2(j,i)
               end do
            end do
            exit
         end if
      end do

      call mpi_barrier(mpi_comm_world,ierror)
  !
  ! CALCULATE ENSEMBLE SPREAD
      call this%ens_spread_dualres_regional(mype,en_perts,nelen,en_bar)
      call mpi_barrier(mpi_comm_world,ierror)
  !
  ! CONVERT ENSEMBLE MEMBERS TO ENSEMBLE PERTURBATIONS
      sig_norm=sqrt(one/max(one,n_ens-one))

      do n=1,n_ens
         do i=1,nelen
            en_perts(n,it)%valuesr4(i)=(en_perts(n,it)%valuesr4(i)-en_bar%values(i))*sig_norm
         end do
      end do

   enddo ! it 4d loop
  !
   call gsi_bundledestroy(en_bar,istatus)
   if(istatus/=0) then
      write(6,*)' in get_fv3_mass_ensperts_fv3: trouble destroying en_bar bundle'
             call stop2(999)
          endif

   return

30 write(6,*) 'get_fv3_mass_ensperts_fv3: open filelist of dynamic failed '
   call stop2(555)
40 write(6,*) 'get_fv3_mass_ensperts_fv3: open filelist of tracer failed '
   call stop2(555)
20 write(6,*) 'get_fv3_mass_ensperts_fv3: read FV3 ens failed ',n
   call stop2(555)

  end subroutine get_fv3_mass_ensperts_fv3

  subroutine general_read_fv3_mass(this,dynmfile,tracfile,g_prsl,g_u,g_v,g_tv,g_rh,g_oz,  &
                                                        g_ql,g_qi,g_qr,g_qs,g_qg,g_qnr,g_w,mype)
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    general_read_fv3_mass  read fv3 model ensemble members
  !   prgmmr: cTong            org: CAPS/OU            date: 2019-04-24
  !
  ! abstract: read ensemble members from the fv3 model netcdf format (dynamic and tracer)
  !           for use with hybrid ensemble option.
  !
  ! program history log:

  !   input argument list:
  !
  !   output argument list:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm RS/6000 SP
  !
  !$$$ end documentation block
   
    use netcdf, only: nf90_nowrite,nf90_inquire
    use netcdf, only: nf90_open,nf90_close
    use netcdf, only: nf90_inq_varid,nf90_inquire_variable,nf90_get_var
    use netcdf, only: nf90_inq_dimid,nf90_inquire_dimension
    use kinds, only: r_kind,r_single,i_kind
    use gridmod, only: lat2,lon2,nsig,itotsub,ijn_s,nlon,nlat,eta1_ll
    use constants, only: zero,one,fv,zero_single,rd_over_cp_mass,one_tenth,h300
    use hybrid_ensemble_parameters, only: grd_ens,q_hyb_ens
    use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
    use guess_grids, only: ntguessig
    use netcdf_mod, only: nc_check
    use mod_fv3_lola, only: fv3_h_to_ll,nya,nxa,fv3uv2earth
    use general_commvars_mod, only: ltosi_s,ltosj_s

    implicit none
  !
  ! Declare passed variables
    class(get_wrf_mass_ensperts_class), intent(inout) :: this
    real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig),intent(out):: &
                                                  g_prsl,g_u,g_v,g_tv,g_rh,g_oz, &
                                                  g_qr,g_qs,g_qg,g_qnr,g_w,      &
                                                  g_ql,g_qi

    character(255),intent(in):: dynmfile
    character(255),intent(in):: tracfile
  !
  ! Declare local parameters
    real(r_kind),parameter:: r0_01 = 0.01_r_kind
    real(r_kind),parameter:: r10   = 10.0_r_kind
    real(r_kind),parameter:: r100  = 100.0_r_kind
  !
  ! Declare local variables
    character(len=128) :: name
    integer(i_kind),allocatable,dimension(:):: dim_id,dim,ijns,displss
    real(r_kind),allocatable,dimension(:,:,:)::tempU_3d,tempV_3d,temp_prsl,tempW_3d,temp_3d
    real(r_kind),allocatable,dimension(:,:,:)::prsl_3d,delp_3d,ts_3d,tv_3d,q_3d
    real(r_kind),allocatable,dimension(:):: work_u ,work_v ,work_tv,work_rh,work_p
    real(r_kind),allocatable,dimension(:):: work_qr,work_qs,work_qg,work_qnr
    real(r_kind),allocatable,dimension(:):: work_w
    real(r_kind),allocatable,dimension(:):: work_ql,work_qi
    real(r_kind),allocatable,dimension(:,:):: a
    real(r_kind),allocatable,dimension(:,:):: u,v
    real(r_kind),allocatable,dimension(:,:,:)::tsn,prsl,qst

    integer(i_kind):: n,ns,nx,ny,nz,i,j,k,d_max,file_id,var_id,ndim,mype,nzp1,ir,iret,len,ii,jj,kk
    integer(i_kind):: Time_id,s_n_id,w_e_id,b_t_id,s_n_stag_id,w_e_stag_id,b_t_stag_id
    integer(i_kind):: Time_len,s_n_len,w_e_len,b_t_len,s_n_stag_len,w_e_stag_len,b_t_stag_len
    integer(i_kind) iderivative
    integer(i_kind) it,ier,istatus
    integer(i_kind) ndimensions,nvariables,nattributes,unlimiteddimid

    real(r_kind):: deltasigma
    real(r_kind) psfc_this_dry,psfc_this

    logical ice

    character(len=24),parameter :: myname_ = 'general_read_fv3_mass'

    allocate(ijns(npe),displss(npe))
    do i=1,npe
       ijns(i)=ijn_s(i)*nsig
    enddo
    displss(1)=0
    do i=2,npe
       displss(i)=displss(i-1)+ ijns(i-1)
    enddo
  !
  ! OPEN ENSEMBLE MEMBER DATA FILE & EXAMINE DIMENSIONS
    if (mype==0) then !
      allocate(work_u   (itotsub*nsig))
      allocate(work_v   (itotsub*nsig))
      allocate(work_tv  (itotsub*nsig))
      allocate(work_rh  (itotsub*nsig))
      allocate(work_qr  (itotsub*nsig))
      allocate(work_qs  (itotsub*nsig))
      allocate(work_qg  (itotsub*nsig))
      allocate(work_qnr (itotsub*nsig))
      allocate(work_p   (itotsub*nsig))
      allocate(work_w   (itotsub*nsig))
      allocate(work_ql  (itotsub*nsig))
      allocate(work_qi  (itotsub*nsig))
     
      call nc_check( nf90_open(trim(dynmfile),nf90_nowrite,file_id),&
          myname_,'open '//trim(dynmfile) )
   !
   ! FV3 DYNAMIC FILE DIMENSIONS
      call nc_check( nf90_inq_dimid(file_id,'Time',Time_id),&
          myname_,'inq_dimid Time '//trim(dynmfile) )
      call nc_check( nf90_inq_dimid(file_id,'yaxis_2',s_n_id),&
          myname_,'inq_dimid south_north '//trim(dynmfile) )
      call nc_check( nf90_inq_dimid(file_id,'xaxis_1',w_e_id),&
          myname_,'inq_dimid west_east '//trim(dynmfile) )
      call nc_check( nf90_inq_dimid(file_id,'zaxis_1',b_t_id),&
          myname_,'inq_dimid bottom_top '//trim(dynmfile) )
      call nc_check( nf90_inq_dimid(file_id,'yaxis_1',s_n_stag_id),&
          myname_,'inq_dimid south_north_stag '//trim(dynmfile) )
      call nc_check( nf90_inq_dimid(file_id,'xaxis_2',w_e_stag_id),&
          myname_,'inq_dimid west_east_stag '//trim(dynmfile) )

      call nc_check( nf90_inquire_dimension(file_id,Time_id,len=Time_len),&
          myname_,'inquire_dimension Time '//trim(dynmfile) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_id,len=s_n_len),&
          myname_,'inquire_dimension south_north '//trim(dynmfile) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_id,len=w_e_len),&
          myname_,'inquire_dimension west_east '//trim(dynmfile) )
      call nc_check( nf90_inquire_dimension(file_id,b_t_id,len=b_t_len),&
          myname_,'inquire_dimension bottom_top '//trim(dynmfile) )
      call nc_check( nf90_inquire_dimension(file_id,s_n_stag_id,len=s_n_stag_len),&
          myname_,'inquire_dimension south_north_stag '//trim(dynmfile) )
      call nc_check( nf90_inquire_dimension(file_id,w_e_stag_id,len=w_e_stag_len),&
          myname_,'inquire_dimension west_east_stag '//trim(dynmfile) )

      nx=w_e_len
      ny=s_n_len
      nz=b_t_len
      if (nx /= grd_ens%nlon .or. ny /= grd_ens%nlat .or. nz /= grd_ens%nsig) then
       print *,'incorrect grid size in netcdf file'
       print *,'nx,ny,nz,nlon,nlat,nsig',nx,ny,nz,grd_ens%nlon,grd_ens%nlat,grd_ens%nsig
       call stop2(999)
      endif

      iret=nf90_inquire(file_id,ndimensions,nvariables,nattributes,unlimiteddimid)

      allocate(dim(ndimensions))
      allocate(a(nya,nxa))

      do k=1,ndimensions
         iret=nf90_inquire_dimension(file_id,k,name,len)
         dim(k)=len
      enddo
      
      allocate(u(dim(1),dim(4)))
      allocate(v(dim(1),dim(4)))
!
! READ FIELDS FROM dynamic file
!   Read u wind (m/s)
      call nc_check( nf90_inq_varid(file_id,'u',var_id),&
          myname_,'inq_varid U '//trim(dynmfile) )
   
      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable U '//trim(dynmfile) )
      allocate(dim_id(ndim))
      
      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable U '//trim(dynmfile) )
      allocate(tempU_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
      
      call nc_check( nf90_get_var(file_id,var_id,tempU_3d),&
          myname_,'get_var U '//trim(dynmfile) )

!   Read v wind (m/s)
      call nc_check( nf90_inq_varid(file_id,'v',var_id),&
          myname_,'inq_varid V '//trim(dynmfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable V '//trim(dynmfile) )
      if(allocated(dim_id    )) deallocate(dim_id    )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable V '//trim(dynmfile) )
      allocate(tempV_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,tempV_3d),&
          myname_,'get_var V '//trim(dynmfile) )

      nz=nsig
      nzp1=nz+1
      do i=1,nz
         ir=nzp1-i
         call fv3uv2earth(tempU_3d(:,:,i),tempV_3d(:,:,i),nx,ny,u,v)
         call fv3_h_to_ll(u,a,nx,ny,nxa,nya)
         kk=0
         do n=1,npe
            ns=displss(n)+(ir-1)*ijn_s(n)
            do j=1,ijn_s(n)
               ns=ns+1
               kk=kk+1
               ii=ltosi_s(kk)
               jj=ltosj_s(kk)
               work_u(ns)=a(ii,jj)
            end do
         end do
      enddo ! i
      print *,'min/max u',minval(work_u),maxval(work_u)
    
      do i=1,nz
         ir=nzp1-i
         call fv3uv2earth(tempU_3d(:,:,i),tempV_3d(:,:,i),nx,ny,u,v)
         call fv3_h_to_ll(v,a,nx,ny,nxa,nya)
         kk=0
         do n=1,npe
            ns=displss(n)+(ir-1)*ijn_s(n)
            do j=1,ijn_s(n)
               ns=ns+1
               kk=kk+1
               ii=ltosi_s(kk)
               jj=ltosj_s(kk)
               work_v(ns)=a(ii,jj)
            end do
         end do
      enddo ! i
      print *,'min/max v',minval(work_v),maxval(work_v)
      deallocate(tempU_3d,tempV_3d,dim_id)

!   Read w wind (m/s)
      call nc_check( nf90_inq_varid(file_id,'W',var_id),&
          myname_,'inq_varid W '//trim(dynmfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable W '//trim(dynmfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable W '//trim(dynmfile) )
      allocate(tempW_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,tempW_3d),&
          myname_,'get_var W '//trim(dynmfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(tempW_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_w(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max w',minval(work_w),maxval(work_w)
      deallocate(tempW_3d,dim_id)

!   Read tsen (K)
      call nc_check( nf90_inq_varid(file_id,'T',var_id),&
          myname_,'inq_varid T '//trim(dynmfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable T '//trim(dynmfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable T '//trim(dynmfile) )
      allocate(ts_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,ts_3d),&
          myname_,'get_var T '//trim(dynmfile) )
      deallocate(dim_id)

!   Read delp (pa), and transfer it to p 
      call nc_check( nf90_inq_varid(file_id,'delp',var_id),&
          myname_,'inq_varid DELP '//trim(dynmfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable DELP '//trim(dynmfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable DELP '//trim(dynmfile) )
      allocate(delp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,delp_3d),&
          myname_,'get_var DELP '//trim(dynmfile) )

      allocate(temp_prsl(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))+1))
      allocate(  prsl_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
      allocate(     prsl(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
      temp_prsl(:,:,1)=eta1_ll(nsig+1)
      do i=2,nsig+1
         temp_prsl(:,:,i)=delp_3d(:,:,i-1)*0.001_r_kind+temp_prsl(:,:,i-1)
         prsl_3d(:,:,i-1)=temp_prsl(:,:,i) ! Save 3D pressure field (kPa) for later RH calculation
         prsl(:,:,nsig+2-i)=temp_prsl(:,:,i) ! UPRIGHT prsl, needed by RH calculation
      enddo
      deallocate(delp_3d,temp_prsl)
      call nc_check( nf90_close(file_id), myname_,'close '//trim(dynmfile) ) 
      
      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(prsl_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_p(ns)=a(ii,jj)
             end do
          end do
          IF(i==nz) print *,'min/max ps',minval(a)*1000,maxval(a)*1000
      enddo ! i      
      deallocate(dim_id)
      deallocate(dim,a)
!
! READ FIELDS FROM tracer file
      call nc_check( nf90_open(trim(tracfile),nf90_nowrite,file_id),&
          myname_,'open '//trim(tracfile) )
      iret=nf90_inquire(file_id,ndimensions,nvariables,nattributes,unlimiteddimid)

      allocate(dim(ndimensions))
      allocate(a(nlat,nlon))

      do k=1,ndimensions
         iret=nf90_inquire_dimension(file_id,k,name,len)
         dim(k)=len
      enddo
!   Read q (specific humidity)
      call nc_check( nf90_inq_varid(file_id,'sphum',var_id),&
          myname_,'inq_varid SPHUM '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable SPHUM '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable SPHUM '//trim(tracfile) )
      allocate(q_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,q_3d),&
          myname_,'get_var SPHUM '//trim(tracfile) )

      ! Calculating tv with ts and q
      allocate( tv_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
      allocate(   tsn(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))
      do k=1,nsig
        do i=1,dim(dim_id(1))
          do j=1,dim(dim_id(2))
             tv_3d(i,j,k)=ts_3d(i,j,k)*(one+fv*q_3d(i,j,k))
             ! recompute sensible temp from virtual temp
             tsn(i,j,nsig+1-k)= tv_3d(i,j,k)/(one+fv*max(zero,q_3d(i,j,k))) ! UPRIGHT tsn, needed by RH calculation
          enddo
        enddo
      enddo

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(tv_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_tv(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max tv',minval(work_tv),maxval(work_tv)

      ! Calculating RH based on q, ts, and prsl if RH option
      if (.not.q_hyb_ens) then
        allocate(qst(dim(dim_id(1)),dim(dim_id(2)),nsig)) ! UPRIGHT qst, needed by RH calculation
        ice=.true.
        iderivative=0
        call genqsat(qst,tsn,prsl,dim(dim_id(1)),dim(dim_id(2)),nsig,ice,iderivative)
        do k=1,nsig
          do i=1,dim(dim_id(1))
            do j=1,dim(dim_id(2))
               q_3d(i,j,k) = q_3d(i,j,k)/qst(i,j,nsig+1-k)
            enddo
          enddo
        enddo
        deallocate(qst,prsl,tsn)
!        print *,'min/max rh',minval(q_3d),maxval(q_3d)
      else
!        print *,'min/max q ',minval(q_3d),maxval(q_3d)
      end if

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(q_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_rh(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      if (.not.q_hyb_ens) then
        print *,'min/max rh',minval(work_rh),maxval(work_rh)
      else
        print *,'min/max q ',minval(work_rh),maxval(work_rh)
      endif
      deallocate(dim_id)
      deallocate(ts_3d,tv_3d,prsl_3d,q_3d)

!   read ql (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'liq_wat',var_id),&
          myname_,'inq_varid LIQ_WAT '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable LIQ_WAT '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable LIQ_WAT '//trim(tracfile) )
      if(allocated(temp_3d   )) deallocate(temp_3d   )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var LIQ_WAT '//trim(tracfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(temp_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_ql(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max ql',minval(work_ql),maxval(work_ql)
      deallocate(temp_3d,dim_id)

!   read qi (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'ice_wat',var_id),&
          myname_,'inq_varid ICE_WAT '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable ICE_WAT '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable ICE_WAT '//trim(tracfile) )
      if(allocated(temp_3d   )) deallocate(temp_3d   )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var ICE_WAT '//trim(tracfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(temp_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_qi(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max qi',minval(work_qi),maxval(work_qi)
      deallocate(temp_3d,dim_id)

!   read qr (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'rainwat',var_id),&
          myname_,'inq_varid RAINWAT '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable RAINWAT '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable RAINWAT '//trim(tracfile) )
      if(allocated(temp_3d   )) deallocate(temp_3d   )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var RAINWAT '//trim(tracfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(temp_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_qr(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max qr',minval(work_qr),maxval(work_qr)
      deallocate(temp_3d,dim_id)

!   read qs (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'snowwat',var_id),&
          myname_,'inq_varid SNOWWAT '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable SNOWWAT '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable SNOWWAT '//trim(tracfile) )
      if(allocated(temp_3d   )) deallocate(temp_3d   )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var SNOWWAT '//trim(tracfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(temp_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_qs(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max qs',minval(work_qs),maxval(work_qs)
      deallocate(temp_3d,dim_id)

!   read qg (kg/kg)
      call nc_check( nf90_inq_varid(file_id,'graupel',var_id),&
          myname_,'inq_varid GRAUPEL '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable GRAUPEL '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable GRAUPEL '//trim(tracfile) )
      if(allocated(temp_3d   )) deallocate(temp_3d   )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var GRAUPEL '//trim(tracfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(temp_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_qg(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max qg',minval(work_qg),maxval(work_qg)
      deallocate(temp_3d,dim_id)

!   read qnr (#/kg)
      call nc_check( nf90_inq_varid(file_id,'rain_nc',var_id),&
          myname_,'inq_varid RAIN_NC '//trim(tracfile) )

      call nc_check( nf90_inquire_variable(file_id,var_id,ndims=ndim),&
          myname_,'inquire_variable RAIN_NC '//trim(tracfile) )
      allocate(dim_id(ndim))

      call nc_check( nf90_inquire_variable(file_id,var_id,dimids=dim_id),&
          myname_,'inquire_variable RAIN_NC '//trim(tracfile) )
      if(allocated(temp_3d   )) deallocate(temp_3d   )
      allocate(temp_3d(dim(dim_id(1)),dim(dim_id(2)),dim(dim_id(3))))

      call nc_check( nf90_get_var(file_id,var_id,temp_3d),&
          myname_,'get_var RAIN_NC '//trim(tracfile) )

      do i=1,nz
          ir=nzp1-i
          call fv3_h_to_ll(temp_3d(:,:,i),a,dim(dim_id(1)),dim(dim_id(2)),nlon,nlat)
          kk=0
          do n=1,npe
             ns=displss(n)+(ir-1)*ijn_s(n)
             do j=1,ijn_s(n)
                ns=ns+1
                kk=kk+1
                ii=ltosi_s(kk)
                jj=ltosj_s(kk)
                work_qnr(ns)=a(ii,jj)
             end do
          end do
      enddo ! i
      print *,'min/max qnr',minval(work_qnr),maxval(work_qnr)
      deallocate(temp_3d,dim_id)
      deallocate(dim,a)

    endif

    call mpi_scatterv(work_p,ijns,displss,mpi_rtype,&
            g_prsl,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)  
    call mpi_scatterv(work_u,ijns,displss,mpi_rtype,&
            g_u,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_v,ijns,displss,mpi_rtype,&
            g_v,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_tv,ijns,displss,mpi_rtype,&
            g_tv,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)   
    call mpi_scatterv(work_rh,ijns,displss,mpi_rtype,&
            g_rh,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_qr,ijns,displss,mpi_rtype,&
            g_qr,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_qs,ijns,displss,mpi_rtype,&
            g_qs,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_qg,ijns,displss,mpi_rtype,&
            g_qg,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_qnr,ijns,displss,mpi_rtype,&
            g_qnr,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_w,ijns,displss,mpi_rtype,&
            g_w,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_ql,ijns,displss,mpi_rtype,&
            g_ql,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)
    call mpi_scatterv(work_qi,ijns,displss,mpi_rtype,&
            g_qi,ijns(mype+1),mpi_rtype,0,mpi_comm_world,ierror)

    if (mype==0) deallocate(work_u , work_v,work_tv, work_p)
    if (mype==0) deallocate(work_rh,work_qr,work_qs,work_qg,work_qnr) 
    if (mype==0) deallocate(work_w)
    if (mype==0) deallocate(work_ql,work_qi)
    ! for now, don't do anything with oz
    g_oz = 0.

  end subroutine general_read_fv3_mass

! --- CAPS --->
  subroutine calc_num_rms_caps(qr_4d_temp,qs_4d_temp,qg_4d_temp,log_qx_pval)
   use hybrid_ensemble_parameters, only: n_ens,grd_ens
   use kinds, only: r_kind
   use guess_grids, only: ges_tsen  !chenll
   use guess_grids, only: nfldsig   !chenll
   use constants, only: zero

   real(r_kind),dimension(:,:,:,:),intent(in)                    :: qr_4d_temp,qs_4d_temp,qg_4d_temp
   real(r_kind),                   intent(in)                    :: log_qx_pval
   real(r_kind),dimension(grd_ens%lat2,grd_ens%lon2,grd_ens%nsig):: qr_mean_tmp,qs_mean_tmp,qg_mean_tmp
   real(r_kind),dimension(31)    :: rms_qr_tmp,rms_qs_tmp,rms_qg_tmp
   real(r_kind),dimension(31)    :: num_qr_tmp,num_qs_tmp,num_qg_tmp
   real(r_kind)                  :: ges_tsen_tmp
   real(r_kind)                  :: qval_tmp
   integer                       :: i,j,k,n

! initialization
   do k=1,grd_ens%nsig
      do i=1,grd_ens%lon2
         do j=1,grd_ens%lat2
            qr_mean_tmp(j,i,k)=zero
            qs_mean_tmp(j,i,k)=zero
            qg_mean_tmp(j,i,k)=zero
         end do
      end do
   end do
   do i=1,31
      rms_qr_tmp(i)=zero
      rms_qs_tmp(i)=zero
      rms_qg_tmp(i)=zero
      num_qr_tmp(i)=zero
      num_qs_tmp(i)=zero
      num_qg_tmp(i)=zero
   end do
   qval_tmp=(5.0E-4_r_kind**log_qx_pval-1.0_r_kind)/log_qx_pval
!
! add by chenll
      !calculate ensemble mean of qr qs qg
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           qr_mean_tmp(j,i,k)=qr_mean_tmp(j,i,k)+qr_4d_temp(j,i,k,n)
           qs_mean_tmp(j,i,k)=qs_mean_tmp(j,i,k)+qs_4d_temp(j,i,k,n)
           qg_mean_tmp(j,i,k)=qg_mean_tmp(j,i,k)+qg_4d_temp(j,i,k,n)
         end do
       end do
     end do
   end do
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
           qr_mean_tmp(j,i,k)=qr_mean_tmp(j,i,k)/float(n_ens)
           qs_mean_tmp(j,i,k)=qs_mean_tmp(j,i,k)/float(n_ens)
           qg_mean_tmp(j,i,k)=qg_mean_tmp(j,i,k)/float(n_ens)
           !write(6,*) 'qr_mean_tmp=',qr_mean_tmp(j,i,k)
       end do
     end do
   end do

   !calculate rms for qr 
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         !write(6,*) 'ges(i,j,k£©=', ges_tsen(j,i,k,nfldsig)
         do n=1,n_ens 
           !ges_tsen_tmp=ges_tsen(j,i,k,nfldsig)-273.15_r_kind  
           !if((ges_tsen_tmp .le. -62.0_r_kind)&        
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -62.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(1)=rms_qr_tmp(1)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2            
              num_qr_tmp(1)=num_qr_tmp(1)+1.0_r_kind              
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(1)=sqrt(rms_qr_tmp(1)/num_qr_tmp(1))


   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -62.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -54.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(2)=rms_qr_tmp(2)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(2)=num_qr_tmp(2)+1.0_r_kind
            end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(2)=sqrt(rms_qr_tmp(2)/num_qr_tmp(2))
   
   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -54.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -46.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(3)=rms_qr_tmp(3)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2              
              num_qr_tmp(3)=num_qr_tmp(3)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(3)=sqrt(rms_qr_tmp(3)/num_qr_tmp(3))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -46.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -38.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(4)=rms_qr_tmp(4)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(4)=num_qr_tmp(4)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(4)=sqrt(rms_qr_tmp(4)/num_qr_tmp(4))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -38.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -30.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(5)=rms_qr_tmp(5)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(5)=num_qr_tmp(5)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(5)=sqrt(rms_qr_tmp(5)/num_qr_tmp(5))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -30.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -26.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(6)=rms_qr_tmp(6)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(6)=num_qr_tmp(6)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(6)=sqrt(rms_qr_tmp(6)/num_qr_tmp(6))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -26.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -22.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(7)=rms_qr_tmp(7)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(7)=num_qr_tmp(7)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(7)=sqrt(rms_qr_tmp(7)/num_qr_tmp(7))

         
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -22.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -18.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(8)=rms_qr_tmp(8)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(8)=num_qr_tmp(8)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(8)=sqrt(rms_qr_tmp(8)/num_qr_tmp(8))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -18.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -14.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(9)=rms_qr_tmp(9)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(9)=num_qr_tmp(9)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(9)=sqrt(rms_qr_tmp(9)/num_qr_tmp(9))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -14.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -10.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(10)=rms_qr_tmp(10)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(10)=num_qr_tmp(10)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(10)=sqrt(rms_qr_tmp(10)/num_qr_tmp(10))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -10.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -8.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(11)=rms_qr_tmp(11)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(11)=num_qr_tmp(11)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(11)=sqrt(rms_qr_tmp(11)/num_qr_tmp(11))


      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -8.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -6.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(12)=rms_qr_tmp(12)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(12)=num_qr_tmp(12)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(12)=sqrt(rms_qr_tmp(12)/num_qr_tmp(12))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -6.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -4.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(13)=rms_qr_tmp(13)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(13)=num_qr_tmp(13)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(13)=sqrt(rms_qr_tmp(13)/num_qr_tmp(13))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -4.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -2.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(14)=rms_qr_tmp(14)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(14)=num_qr_tmp(14)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(14)=sqrt(rms_qr_tmp(14)/num_qr_tmp(14))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -2.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 0.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(15)=rms_qr_tmp(15)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(15)=num_qr_tmp(15)+1.0_r_kind
            end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(15)=sqrt(rms_qr_tmp(15)/num_qr_tmp(15))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 0.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 2.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(16)=rms_qr_tmp(16)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(16)=num_qr_tmp(16)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(16)=sqrt(rms_qr_tmp(16)/num_qr_tmp(16))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 2.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 4.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(17)=rms_qr_tmp(17)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(17)=num_qr_tmp(17)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(17)=sqrt(rms_qr_tmp(17)/num_qr_tmp(17))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 4.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 6.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(18)=rms_qr_tmp(18)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(18)=num_qr_tmp(18)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(18)=sqrt(rms_qr_tmp(18)/num_qr_tmp(18))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 6.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 8.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(19)=rms_qr_tmp(19)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(19)=num_qr_tmp(19)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(19)=sqrt(rms_qr_tmp(19)/num_qr_tmp(19))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 8.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 10.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(20)=rms_qr_tmp(20)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(20)=num_qr_tmp(20)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(20)=sqrt(rms_qr_tmp(20)/num_qr_tmp(20))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 10.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 12.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(21)=rms_qr_tmp(21)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2              
              num_qr_tmp(21)=num_qr_tmp(21)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(21)=sqrt(rms_qr_tmp(21)/num_qr_tmp(21))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 12.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 14.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(22)=rms_qr_tmp(22)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(22)=num_qr_tmp(22)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(22)=sqrt(rms_qr_tmp(22)/num_qr_tmp(22))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 14.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 16.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(23)=rms_qr_tmp(23)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(23)=num_qr_tmp(23)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(23)=sqrt(rms_qr_tmp(23)/num_qr_tmp(23))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 16.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 18.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(24)=rms_qr_tmp(24)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(24)=num_qr_tmp(24)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(24)=sqrt(rms_qr_tmp(24)/num_qr_tmp(24))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 18.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 20.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(25)=rms_qr_tmp(25)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(25)=num_qr_tmp(25)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(25)=sqrt(rms_qr_tmp(25)/num_qr_tmp(25))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 20.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 22.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(26)=rms_qr_tmp(26)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(26)=num_qr_tmp(26)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(26)=sqrt(rms_qr_tmp(26)/num_qr_tmp(26))

         
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 22.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 24.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(27)=rms_qr_tmp(27)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(27)=num_qr_tmp(27)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(27)=sqrt(rms_qr_tmp(27)/num_qr_tmp(27))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 24.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 26.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(28)=rms_qr_tmp(28)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(28)=num_qr_tmp(28)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(28)=sqrt(rms_qr_tmp(28)/num_qr_tmp(28))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 26.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 28.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(29)=rms_qr_tmp(29)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(29)=num_qr_tmp(29)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(29)=sqrt(rms_qr_tmp(29)/num_qr_tmp(29))

       
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 28.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 30.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(30)=rms_qr_tmp(30)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(30)=num_qr_tmp(30)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(30)=sqrt(rms_qr_tmp(30)/num_qr_tmp(30))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 30.0_r_kind)&
            .and. (qr_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qr_tmp(31)=rms_qr_tmp(31)+(qr_4d_temp(j,i,k,n)-qr_mean_tmp(j,i,k))**2
              num_qr_tmp(31)=num_qr_tmp(31)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qr_tmp(31)=sqrt(rms_qr_tmp(31)/num_qr_tmp(31))

   
   
   
!calculate rms for qs 
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         !write(6,*) 'ges(i,j,k£©=', ges_tsen(j,i,k,nfldsig)
         do n=1,n_ens 
           !ges_tsen_tmp=ges_tsen(j,i,k,nfldsig)-273.15_r_kind  
           !if((ges_tsen_tmp .le. -62.0_r_kind)&        
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -62.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(1)=rms_qs_tmp(1)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2            
              num_qs_tmp(1)=num_qs_tmp(1)+1.0_r_kind              
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(1)=sqrt(rms_qs_tmp(1)/num_qs_tmp(1))


   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -62.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -54.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(2)=rms_qs_tmp(2)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(2)=num_qs_tmp(2)+1.0_r_kind
            end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(2)=sqrt(rms_qs_tmp(2)/num_qs_tmp(2))
   
   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -54.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -46.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(3)=rms_qs_tmp(3)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2              
              num_qs_tmp(3)=num_qs_tmp(3)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(3)=sqrt(rms_qs_tmp(3)/num_qs_tmp(3))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -46.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -38.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(4)=rms_qs_tmp(4)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(4)=num_qs_tmp(4)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(4)=sqrt(rms_qs_tmp(4)/num_qs_tmp(4))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -38.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -30.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(5)=rms_qs_tmp(5)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(5)=num_qs_tmp(5)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(5)=sqrt(rms_qs_tmp(5)/num_qs_tmp(5))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -30.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -26.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(6)=rms_qs_tmp(6)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(6)=num_qs_tmp(6)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(6)=sqrt(rms_qs_tmp(6)/num_qs_tmp(6))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -26.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -22.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(7)=rms_qs_tmp(7)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(7)=num_qs_tmp(7)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(7)=sqrt(rms_qs_tmp(7)/num_qs_tmp(7))

         
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -22.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -18.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(8)=rms_qs_tmp(8)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(8)=num_qs_tmp(8)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(8)=sqrt(rms_qs_tmp(8)/num_qs_tmp(8))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -18.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -14.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(9)=rms_qs_tmp(9)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(9)=num_qs_tmp(9)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(9)=sqrt(rms_qs_tmp(9)/num_qs_tmp(9))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -14.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -10.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(10)=rms_qs_tmp(10)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(10)=num_qs_tmp(10)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(10)=sqrt(rms_qs_tmp(10)/num_qs_tmp(10))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -10.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -8.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(11)=rms_qs_tmp(11)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(11)=num_qs_tmp(11)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(11)=sqrt(rms_qs_tmp(11)/num_qs_tmp(11))


      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -8.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -6.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(12)=rms_qs_tmp(12)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(12)=num_qs_tmp(12)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(12)=sqrt(rms_qs_tmp(12)/num_qs_tmp(12))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -6.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -4.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(13)=rms_qs_tmp(13)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(13)=num_qs_tmp(13)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(13)=sqrt(rms_qs_tmp(13)/num_qs_tmp(13))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -4.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -2.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(14)=rms_qs_tmp(14)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(14)=num_qs_tmp(14)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(14)=sqrt(rms_qs_tmp(14)/num_qs_tmp(14))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -2.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 0.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(15)=rms_qs_tmp(15)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(15)=num_qs_tmp(15)+1.0_r_kind
            end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(15)=sqrt(rms_qs_tmp(15)/num_qs_tmp(15))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 0.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 2.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(16)=rms_qs_tmp(16)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(16)=num_qs_tmp(16)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(16)=sqrt(rms_qs_tmp(16)/num_qs_tmp(16))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 2.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 4.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(17)=rms_qs_tmp(17)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(17)=num_qs_tmp(17)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(17)=sqrt(rms_qs_tmp(17)/num_qs_tmp(17))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 4.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 6.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(18)=rms_qs_tmp(18)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(18)=num_qs_tmp(18)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(18)=sqrt(rms_qs_tmp(18)/num_qs_tmp(18))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 6.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 8.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(19)=rms_qs_tmp(19)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(19)=num_qs_tmp(19)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(19)=sqrt(rms_qs_tmp(19)/num_qs_tmp(19))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 8.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 10.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(20)=rms_qs_tmp(20)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(20)=num_qs_tmp(20)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(20)=sqrt(rms_qs_tmp(20)/num_qs_tmp(20))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 10.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 12.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(21)=rms_qs_tmp(21)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2              
              num_qs_tmp(21)=num_qs_tmp(21)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(21)=sqrt(rms_qs_tmp(21)/num_qs_tmp(21))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 12.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 14.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(22)=rms_qs_tmp(22)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(22)=num_qs_tmp(22)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(22)=sqrt(rms_qs_tmp(22)/num_qs_tmp(22))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 14.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 16.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(23)=rms_qs_tmp(23)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(23)=num_qs_tmp(23)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(23)=sqrt(rms_qs_tmp(23)/num_qs_tmp(23))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 16.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 18.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(24)=rms_qs_tmp(24)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(24)=num_qs_tmp(24)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(24)=sqrt(rms_qs_tmp(24)/num_qs_tmp(24))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 18.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 20.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(25)=rms_qs_tmp(25)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(25)=num_qs_tmp(25)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(25)=sqrt(rms_qs_tmp(25)/num_qs_tmp(25))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 20.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 22.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(26)=rms_qs_tmp(26)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(26)=num_qs_tmp(26)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(26)=sqrt(rms_qs_tmp(26)/num_qs_tmp(26))

         
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 22.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 24.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(27)=rms_qs_tmp(27)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(27)=num_qs_tmp(27)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(27)=sqrt(rms_qs_tmp(27)/num_qs_tmp(27))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 24.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 26.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(28)=rms_qs_tmp(28)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(28)=num_qs_tmp(28)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(28)=sqrt(rms_qs_tmp(28)/num_qs_tmp(28))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 26.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 28.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(29)=rms_qs_tmp(29)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(29)=num_qs_tmp(29)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(29)=sqrt(rms_qs_tmp(29)/num_qs_tmp(29))

       
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 28.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 30.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(30)=rms_qs_tmp(30)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(30)=num_qs_tmp(30)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(30)=sqrt(rms_qs_tmp(30)/num_qs_tmp(30))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 30.0_r_kind)&
            .and. (qs_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qs_tmp(31)=rms_qs_tmp(31)+(qs_4d_temp(j,i,k,n)-qs_mean_tmp(j,i,k))**2
              num_qs_tmp(31)=num_qs_tmp(31)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qs_tmp(31)=sqrt(rms_qs_tmp(31)/num_qs_tmp(31))




!calculate rms for qg 
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         !write(6,*) 'ges(i,j,k£©=', ges_tsen(j,i,k,nfldsig)
         do n=1,n_ens 
           !ges_tsen_tmp=ges_tsen(j,i,k,nfldsig)-273.15_r_kind  
           !if((ges_tsen_tmp .le. -62.0_r_kind)&        
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -62.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(1)=rms_qg_tmp(1)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2            
              num_qg_tmp(1)=num_qg_tmp(1)+1.0_r_kind              
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(1)=sqrt(rms_qg_tmp(1)/num_qg_tmp(1))


   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -62.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -54.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(2)=rms_qg_tmp(2)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(2)=num_qg_tmp(2)+1.0_r_kind
            end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(2)=sqrt(rms_qg_tmp(2)/num_qg_tmp(2))
   
   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -54.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -46.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(3)=rms_qg_tmp(3)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2              
              num_qg_tmp(3)=num_qg_tmp(3)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(3)=sqrt(rms_qg_tmp(3)/num_qg_tmp(3))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -46.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -38.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(4)=rms_qg_tmp(4)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(4)=num_qg_tmp(4)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(4)=sqrt(rms_qg_tmp(4)/num_qg_tmp(4))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -38.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -30.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(5)=rms_qg_tmp(5)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(5)=num_qg_tmp(5)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(5)=sqrt(rms_qg_tmp(5)/num_qg_tmp(5))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -30.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -26.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(6)=rms_qg_tmp(6)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(6)=num_qg_tmp(6)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(6)=sqrt(rms_qg_tmp(6)/num_qg_tmp(6))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -26.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -22.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(7)=rms_qg_tmp(7)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(7)=num_qg_tmp(7)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(7)=sqrt(rms_qg_tmp(7)/num_qg_tmp(7))

         
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -22.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -18.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(8)=rms_qg_tmp(8)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(8)=num_qg_tmp(8)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(8)=sqrt(rms_qg_tmp(8)/num_qg_tmp(8))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -18.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -14.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(9)=rms_qg_tmp(9)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(9)=num_qg_tmp(9)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(9)=sqrt(rms_qg_tmp(9)/num_qg_tmp(9))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -14.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -10.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(10)=rms_qg_tmp(10)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(10)=num_qg_tmp(10)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(10)=sqrt(rms_qg_tmp(10)/num_qg_tmp(10))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -10.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -8.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(11)=rms_qg_tmp(11)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(11)=num_qg_tmp(11)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(11)=sqrt(rms_qg_tmp(11)/num_qg_tmp(11))


      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -8.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -6.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(12)=rms_qg_tmp(12)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(12)=num_qg_tmp(12)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(12)=sqrt(rms_qg_tmp(12)/num_qg_tmp(12))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -6.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -4.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(13)=rms_qg_tmp(13)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(13)=num_qg_tmp(13)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(13)=sqrt(rms_qg_tmp(13)/num_qg_tmp(13))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -4.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. -2.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(14)=rms_qg_tmp(14)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(14)=num_qg_tmp(14)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(14)=sqrt(rms_qg_tmp(14)/num_qg_tmp(14))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. -2.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 0.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(15)=rms_qg_tmp(15)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(15)=num_qg_tmp(15)+1.0_r_kind
            end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(15)=sqrt(rms_qg_tmp(15)/num_qg_tmp(15))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 0.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 2.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(16)=rms_qg_tmp(16)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(16)=num_qg_tmp(16)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(16)=sqrt(rms_qg_tmp(16)/num_qg_tmp(16))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 2.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 4.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(17)=rms_qg_tmp(17)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(17)=num_qg_tmp(17)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(17)=sqrt(rms_qg_tmp(17)/num_qg_tmp(17))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 4.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 6.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(18)=rms_qg_tmp(18)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(18)=num_qg_tmp(18)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(18)=sqrt(rms_qg_tmp(18)/num_qg_tmp(18))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 6.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 8.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(19)=rms_qg_tmp(19)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(19)=num_qg_tmp(19)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(19)=sqrt(rms_qg_tmp(19)/num_qg_tmp(19))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 8.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 10.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(20)=rms_qg_tmp(20)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(20)=num_qg_tmp(20)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(20)=sqrt(rms_qg_tmp(20)/num_qg_tmp(20))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 10.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 12.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(21)=rms_qg_tmp(21)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2              
              num_qg_tmp(21)=num_qg_tmp(21)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(21)=sqrt(rms_qg_tmp(21)/num_qg_tmp(21))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 12.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 14.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(22)=rms_qg_tmp(22)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(22)=num_qg_tmp(22)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(22)=sqrt(rms_qg_tmp(22)/num_qg_tmp(22))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 14.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 16.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(23)=rms_qg_tmp(23)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(23)=num_qg_tmp(23)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(23)=sqrt(rms_qg_tmp(23)/num_qg_tmp(23))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 16.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 18.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(24)=rms_qg_tmp(24)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(24)=num_qg_tmp(24)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(24)=sqrt(rms_qg_tmp(24)/num_qg_tmp(24))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 18.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 20.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(25)=rms_qg_tmp(25)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(25)=num_qg_tmp(25)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(25)=sqrt(rms_qg_tmp(25)/num_qg_tmp(25))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 20.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 22.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(26)=rms_qg_tmp(26)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(26)=num_qg_tmp(26)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(26)=sqrt(rms_qg_tmp(26)/num_qg_tmp(26))

         
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 22.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 24.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(27)=rms_qg_tmp(27)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(27)=num_qg_tmp(27)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(27)=sqrt(rms_qg_tmp(27)/num_qg_tmp(27))

   
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 24.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 26.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(28)=rms_qg_tmp(28)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(28)=num_qg_tmp(28)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(28)=sqrt(rms_qg_tmp(28)/num_qg_tmp(28))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 26.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 28.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(29)=rms_qg_tmp(29)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(29)=num_qg_tmp(29)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(29)=sqrt(rms_qg_tmp(29)/num_qg_tmp(29))

       
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 28.0_r_kind)&
            .and. ((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .le. 30.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(30)=rms_qg_tmp(30)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(30)=num_qg_tmp(30)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(30)=sqrt(rms_qg_tmp(30)/num_qg_tmp(30))

      
   do k=1,grd_ens%nsig
     do i=1,grd_ens%lon2
       do j=1,grd_ens%lat2
         do n=1,n_ens
           if(((ges_tsen(j,i,k,nfldsig)-273.15_r_kind) .gt. 30.0_r_kind)&
            .and. (qg_4d_temp(j,i,k,n) .ge. qval_tmp)) then
              rms_qg_tmp(31)=rms_qg_tmp(31)+(qg_4d_temp(j,i,k,n)-qg_mean_tmp(j,i,k))**2
              num_qg_tmp(31)=num_qg_tmp(31)+1.0_r_kind
           end if
         end do
       end do
     end do
   end do
   rms_qg_tmp(31)=sqrt(rms_qg_tmp(31)/num_qg_tmp(31))
  
   !output
   write(6,*)' nownow output rms_qr'
   do i=1,31
     write(6,*) rms_qr_tmp(i) 
   end do
   write(6,*)' now output rms_qs'
   do i=1,31
     write(6,*) rms_qs_tmp(i) 
   end do   
   write(6,*)' now output rms_qg'
   do i=1,31
     write(6,*) rms_qg_tmp(i) 
   end do 
   write(6,*)' now output num_qr'  
   do i=1,31
     write(6,*) num_qr_tmp(i) 
   end do 
   write(6,*)' now output num_qs'   
   do i=1,31
     write(6,*) num_qs_tmp(i) 
   end do  
   write(6,*)' now output num_qg'     
   do i=1,31
     write(6,*) num_qg_tmp(i) 
   end do
  end subroutine calc_num_rms_caps
! <--- CAPS ---
end module get_wrf_mass_ensperts_mod
