program em_2body
   use mc_module
   use dirac_matrices
   use mympi
   use mathtool
   
   implicit none
   real*8, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0
   real*8, parameter :: xmd=1232.25d0,xmn=938.91875d0,xmpi=139.570179d0
   integer*4 :: nw,nev,i,xA,i_fg,np,ne,j,i_fsi,nwlk,nwfold,i_intf,i_exc,i_dir
   integer*8, allocatable :: irn(:),irn0(:)   
   real*8 :: wmax,qval,xpf,hw
   real*8, allocatable :: resp(:,:,:),resp_err(:,:,:),w(:)
   real*8, allocatable :: resp_fold(:,:,:),wfold(:),fold(:)
   character*50 :: fname,intf_char,en_char,qval_char,int_char,dec_char
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
      read(5,*) i_exc
      read(5,*) i_dir
      read(5,*) xpf
      read(5,*) xA
      read(5,*) i_fg
      read(5,*) nk_fname
      read(5,*) np,ne
      read(5,*) i_fsi    

      !fname='bodek_results/C12_FSI_EM_1b_QMC_q570.out'
      fname='test_FG_delta_exactR_q500_Ca.out'
      fname=trim(fname)
      open(unit=7, file=fname)
   endif   

   call bcast(nev)
   call bcast(nwlk)
   call bcast(wmax)
   call bcast(nw)
   call bcast(qval)
   call bcast(i_intf)
   call bcast(i_exc)
   call bcast(i_dir)
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

   allocate(w(nw),resp(2,5,nw),resp_err(2,5,nw))
   call dirac_matrices_in(xmd,xmn,xmpi)

   !Initialize currents and spinors, i_intf determines 1 and 2b intf calculation
   call mc_init(i_intf,i_exc,i_dir,i_fg,i_fsi,irn,nev,nwlk,xpf,qval,xmpi,xmd,xmn,xA,np,ne,nk_fname)


   
   write(6,*) 'R_L, R_T'
   hw=wmax/dble(nw)
   do i=1,nw
      w(i)=dble(i)*hw
      call mc_eval(w(i),resp(:,:,i),resp_err(:,:,i))
      if (myrank().eq.0) then      
      !! I want only the electromagnetic
        write(6,*) 'wval and ph sp', w(i), resp(i_intf+1,1,i),resp(i_intf+1,4,i)
        write(7,*)w(i), resp(i_intf+1,1,i),resp(i_intf+1,4,i)
        !flush(7)
      endif  
   enddo
   !close(7)




   if(i_fsi.eq.1) then
       open(unit=21,file='../../FSI/folding.in',status='unknown',form='formatted')
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
          write(7,*) w(i),resp_fold(i_intf+1,1,i),resp_fold(i_intf+1,4,i)
       enddo
    endif
    close(7)


      tf=MPI_Wtime()
   if (myrank().eq.0) then
      write(6,*)'Elapsed time is',tf-ti
   endif
   call done()   

end program
