program ew_2body
   use mc_module
   use dirac_matrices
   use mympi
   use mathtool
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1236.0d0,xmn=939.565d0,xmp=938.272d0,xmpi=139.d0
   integer*4 :: nw,nev,i,xA,i_intf,i_fg,np,ne,j,i_fsi,nwlk,nwfold
   integer*8, allocatable :: irn(:),irn0(:)   
   real*8 :: wmax,qval,xpf,hw,xmavg
   real*8, allocatable :: resp(:,:,:),resp_err(:,:,:),sig(:,:),sig_err(:,:),w(:)
   real*8, allocatable :: resp_fold(:,:,:),wfold(:),fold(:)
   character*50 :: fname,fname2,intf_char,en_char,theta_char,int_char,dec_char
   real*8 :: ti,tf,TA,dwfold,ffold,hwfold    
   character*40 :: nk_fname

   ! initialize mpi
   call init0()


   if (myrank().eq.0) then
      read(5,*) nev
      read(5,*) nwlk
      read(5,*) qval     
      read(5,*) wmax,nw
      read(5,*) i_intf
      read(5,*) xpf
      read(5,*) xA
      read(5,*) i_fg
      read(5,*) nk_fname
      read(5,*) np,ne
      read(5,*) i_fsi      

   
      !write(en_char,'(i4)') int(qval)
      !en_char=adjustl(en_char)
      !en_char=trim(en_char)

      fname='2b_q500_CC.out'
      fname=trim(fname)

      fname2='2b_q500_CC_FSI.out'
      fname2=trim(fname2)

      open(unit=7, file=fname)
      if(i_fsi.eq.1) then
         open(unit=14, file=fname2)
      endif

   endif   


   call bcast(nev)
   call bcast(wmax)
   call bcast(nw)
   call bcast(qval)
   call bcast(i_intf)
   call bcast(xpf)
   call bcast(xA)
   call bcast(i_fg)
   call bcast(nk_fname)
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
   allocate(w(nw),resp(2,5,nw),resp_err(2,5,nw))
   call dirac_matrices_in(xmd,xmavg,xmpi)

   !compute the responses
   call mc_init(i_intf,i_fg,i_fsi,irn,nev,nwlk,xpf,qval,xmpi,xmd,xmavg,xA,np,ne,nk_fname)
   
   write(6,*) 'R_cc, R_cl, R_ll, R_t, R_tp'
   !Computing response tensor elements instead of functions
   !write(6,*) 'R_00, R_0Z, R_ZZ, R_XX, R_XY'

   hw=wmax/dble(nw)
   do i=1,nw
      w(i)=dble(i)*hw
      call mc_eval(w(i),resp(:,:,i),resp_err(:,:,i))
      if (myrank().eq.0) then      
        write(6,*) 'xsec:', w(i), resp(i_intf+1,1,i), resp(i_intf+1,2,i), resp(i_intf+1,3,i), resp(i_intf+1,4,i), resp(i_intf+1,5,i)
        !write(7,*) w(i), resp(i_intf+1,1,i), resp(i_intf+1,2,i), resp(i_intf+1,3,i), resp(i_intf+1,4,i), resp(i_intf+1,5,i)
        !flush(7)
      endif  
   enddo

   !close(7)


   if(i_fsi.eq.1) then
       open(unit=21,file='FSI/folding.in',status='unknown',form='formatted')
       read(21,*)TA,hwfold,nwfold
       allocate(resp_fold(2,5,nw),wfold(nwfold),fold(nwfold))
       do i=1,nwfold
          read(21,*) wfold(i),fold(i)
       enddo
       close(21)
       write(6,*)'TA=',TA
       TA=1.0d0-2.0d0*sum(fold(:))*hwfold
       write(6,*)'norm folding',TA
       do i=1,nw 
          resp_fold(:,:,i)=0.0d0
          do j=1,nw
             dwfold=abs(w(i)-w(j))
             if (dwfold.le.wfold(nwfold)) then
                call interpolint(wfold,fold,nwfold,dwfold,ffold,3)
                resp_fold(:,:,i)=resp_fold(:,:,i)+resp(:,:,j)*ffold*hw!+rl(j)*ffold*hw
             endif
          enddo
          resp_fold(:,:,i)=resp_fold(:,:,i)+TA*resp(:,:,i)
       enddo
       do i=1,nw
          !write(6,*) w(i),sig_fold(1,i),sig_fold(2,i)
          write(14,*) w(i),resp_fold(i_intf+1,1,i), resp_fold(i_intf+1,2,i) &
          & , resp_fold(i_intf+1,3,i), resp_fold(i_intf+1,4,i), resp_fold(i_intf+1,5,i)
       enddo
    endif
    close(14)


      tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()   

end program
