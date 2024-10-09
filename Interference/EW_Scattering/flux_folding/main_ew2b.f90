program ew_2body
   use mc_module
   use dirac_matrices
   use mympi
   use mathtool
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=939.565d0,xmp=938.272d0,xmpi=139.d0
   integer*4 :: nw,nev,i,xA,i_intf,i_fg,np,ne,j,i_fsi,nwlk,nwfold,n_cos_bins,nenu
   integer*8, allocatable :: irn(:),irn0(:)   
   real*8 :: wmax,ee,thetalept,xpf,hw,xmavg,cos_thetamin,cos_thetamax
   real*8, allocatable :: sig(:,:),sig_err(:,:),w(:)
   real*8, allocatable :: sig_fold(:,:),wfold(:),fold(:)
   character*50 :: fname,intf_char,en_char,theta_char,int_char,dec_char
   real*8 :: ti,tf,TA,dwfold,ffold,hwfold    
   character*40 :: nk_fname, flux_fname
   !open(unit=7, file='xsec_ew.out')
   ! initialize mpi
   call init0()


   if (myrank().eq.0) then
      read(5,*) nev
      read(5,*) nwlk 
      read(5,*) cos_thetamin,cos_thetamax,n_cos_bins    
      read(5,*) wmax,nw
      read(5,*) i_intf
      read(5,*) xpf
      read(5,*) xA
      read(5,*) i_fg
      read(5,*) nk_fname
      read(5,*) flux_fname
      read(5,*) nenu
      read(5,*) np,ne
      read(5,*) i_fsi      

      write(int_char,'(i3)') int(thetalept)
      write(dec_char,'(i1)') int(mod(thetalept*10.0d0,10.0d0))

      theta_char=trim(int_char)//'p'//trim(dec_char)
      theta_char=adjustl(theta_char)
      theta_char=trim(theta_char)
      write(en_char,'(i4)') int(ee)
      en_char=adjustl(en_char)
      en_char=trim(en_char)

      fname='MB_C12_QMC_fluxfold_2b_08_09_lattice.out'
      fname=trim(fname)

      open(unit=7, file=fname)
      !open(unit=7, file='test.out')
   endif   


   call bcast(nev)
   call bcast(wmax)
   call bcast(nw)
   call bcast(cos_thetamin)
   call bcast(cos_thetamax)
   call bcast(n_cos_bins)
   call bcast(i_intf)
   call bcast(xpf)
   call bcast(xA)
   call bcast(i_fg)
   call bcast(nk_fname)
   call bcast(flux_fname)
   call bcast(nenu)
   call bcast(np)
   call bcast(ne)
   call bcast(i_fsi)
   
   ti=MPI_Wtime()



   allocate(irn0(nwlk))
   do i=1,nwlk
       irn0(i)=19+i
    enddo
    if (myrank().eq.0) then
       write (6,'(''number of cpus ='',t50,i10)') nproc()
       if (mod(nwlk,nproc()).ne.0) then
          write(6,*)'Error: nwalk must me a multiple of nproc'
          stop
       endif
    endif
    nwlk=nwlk/nproc()
    
    allocate(irn(nwlk))
    irn(:)=irn0(myrank()*nwlk+1:myrank()*nwlk+nwlk)

   xmavg = (xmn+xmp)/2.0d0
   allocate(w(nw),sig(2,nw),sig_err(2,nw))
   call dirac_matrices_in(xmd,xmavg,xmpi)

   !compute the cross section 
   call mc_init(i_intf,i_fg,i_fsi,irn,nev,nwlk,xpf, &
      &   cos_thetamin,cos_thetamax,n_cos_bins,xmpi,xmd,xmavg, &
      &   xA,np,ne,nk_fname,flux_fname,nenu)
   
   hw=wmax/dble(nw)
   do i=1,nw
      w(i)=(-0.5d0+dble(i))*hw
      call mc_eval(w(i),sig(:,i),sig_err(:,i))
      if (myrank().eq.0) then      
        write(6,*) 'xsec:', w(i), sig(1,i), sig(2,i)
        write(7,*) w(i), sig(1,i), sig(2,i)
      endif  
   enddo

   close(7)


   if(i_fsi.eq.1) then
       open(unit=21,file='Responses/FSI/folding.in',status='unknown',form='formatted')
       read(21,*)TA,hwfold,nwfold
       allocate(sig_fold(2,nw),wfold(nwfold),fold(nwfold))
       do i=1,nwfold
          read(21,*) wfold(i),fold(i)
       enddo
       close(21)
       write(6,*)'TA=',TA
       TA=1.0d0-2.0d0*sum(fold(:))*hwfold
       write(6,*)'norm folding',TA
       do i=1,nw 
          sig_fold(:,i)=0.0d0
          do j=1,nw
             dwfold=abs(w(i)-w(j))
             if (dwfold.le.wfold(nwfold)) then
                call interpolint(wfold,fold,nwfold,dwfold,ffold,3)
                sig_fold(:,i)=sig_fold(:,i)+sig(:,j)*ffold*hw!+rl(j)*ffold*hw
             endif
          enddo
          sig_fold(:,i)=sig_fold(:,i)+TA*sig(:,i)
       enddo
       do i=1,nw
          write(6,*) w(i),sig_fold(1,i),sig_fold(2,i)
          write(14,*) w(i),sig_fold(1,i),sig_fold(2,i)
       enddo
    endif
    close(14)


      tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()   

end program
