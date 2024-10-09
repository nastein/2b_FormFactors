program ew_2body
   use mc_module
   use dirac_matrices
   use mympi
   use mathtool
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=939.565d0,xmp=938.272d0,xmpi=139.d0
   integer*4 :: nw,nev,i,xA,i_intf,i_fg,np,ne,j,i_fsi,nwlk,nwfold,i_exc,i_dir
   integer*8, allocatable :: irn(:),irn0(:)   
   real*8 :: wmax,ee,thetalept,xpf,hw,xmavg
   real*8, allocatable :: sig(:,:),sig_err(:,:),w(:)
   real*8, allocatable :: sig_fold(:,:),wfold(:),fold(:)
   character*50 :: fname,intf_char,en_char,theta_char,int_char,dec_char
   real*8 :: ti,tf,TA,dwfold,ffold,hwfold    
   character*40 :: nk_fname
   !open(unit=7, file='xsec_ew.out')
   ! initialize mpi
   call init0()


   if (myrank().eq.0) then
      read(5,*) nev
      read(5,*) nwlk
      read(5,*) ee,thetalept      
      read(5,*) wmax,nw
      read(5,*) i_intf
      read(5,*) i_exc
      read(5,*) i_dir
      read(5,*) xpf
      read(5,*) xA
      read(5,*) i_fg
      read(5,*) nk_fname
      read(5,*) np,ne   

      write(int_char,'(i3)') int(thetalept)
      write(dec_char,'(i1)') int(mod(thetalept*10.0d0,10.0d0))

      theta_char=trim(int_char)//'p'//trim(dec_char)
      theta_char=adjustl(theta_char)
      theta_char=trim(theta_char)
      write(en_char,'(i4)') int(ee)
      en_char=adjustl(en_char)
      en_char=trim(en_char)

      !fname='test.out'
      !fname='C12_EW_12b_'//trim(en_char)//'_'//trim(theta_char)//'_nodir.out'
      fname='test.out'
      fname=trim(fname)

      open(unit=7, file=fname)
      !open(unit=7, file='test.out')
   endif   


   call bcast(nev)
   call bcast(wmax)
   call bcast(nw)
   call bcast(ee)
   call bcast(thetalept)
   call bcast(i_intf)
   call bcast(i_exc)
   call bcast(i_dir)
   call bcast(xpf)
   call bcast(xA)
   call bcast(i_fg)
   call bcast(nk_fname)
   call bcast(np)
   call bcast(ne)
   
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

    thetalept=thetalept/180.0d0*pi

   xmavg = (xmn+xmp)/2.0d0
   allocate(w(nw),sig(2,nw),sig_err(2,nw))
   call dirac_matrices_in(xmd,xmavg,xmpi,0.0d0,105.658357d0)

   if (myrank().eq.0) then
      write(6,*)'Starting up'
   endif

   !compute the cross section 
   call mc_init(i_intf,i_exc,i_dir,i_fg,0,irn,nev,nwlk,xpf,thetalept,xmpi,xmd,xmavg,xA,np,ne,nk_fname)
   


   hw=(wmax)/dble(nw)
   do i=1,nw
      w(i)= dble(i)*hw
      call mc_eval(ee,w(i),sig(:,i),sig_err(:,i))
      if (myrank().eq.0) then      
        write(6,*) 'xsec:', w(i), sig(1,i), sig(2,i)
        write(7,*) w(i), sig(1,i), sig(2,i)
      endif  
   enddo

   close(7)

      tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()   

end program
