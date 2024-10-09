module mc_module
   implicit none
   integer*4, private, save :: nev,xA,i_fg,i_fsi,np,ne,nwlk,npot,np_del,i_intf
   integer*4, private, save :: i_exc, i_dir
   integer*4, private, parameter :: neq=10000,nvoid=100
   real*8, private, save ::  xpf,xpmax
   real*8, private, save:: xmpi,xmd,xmn,norm,qval
   real*8, private, parameter :: pi=acos(-1.0d0),hbarc=197.327053d0,xmp=938.0d0,ppmax=1.0d0*1.e3
   real*8, private,parameter :: alpha=1.0d0/137.0d0
   real*8, private, allocatable :: pv(:),dp(:,:),ep(:),Pke(:,:),nk(:)
   real*8, private, allocatable :: kin(:),pot(:),pdel(:),pot_del(:)
   integer*8, private, allocatable, save :: irn(:)
contains

subroutine mc_init(i_intf_in,i_exc_in,i_dir_in, &
   &  i_fg_in,i_fsi_in,irn_in,nev_in,nwlk_in, &
   &  xpf_in,qval_in,xmpi_in,xmd_in,xmn_in,xA_in, np_in,ne_in,nk_fname_in)

  use mathtool
   implicit none
   integer*8 :: irn_in(nwlk_in)
   integer*4 :: nev_in,nwlk_in,xA_in,i_fg_in,np_in,i,j,ne_in,np0,ne0,ien,i_intf_in
   integer*4 :: ipot,i_fsi_in,i_exc_in,i_dir_in
   real*8 :: xpf_in,xmpi_in,xmd_in,xmn_in,mlept_in,hp,he,qval_in,dummy
   real*8, allocatable :: pv0(:),dp0(:,:),ep0(:)
   character*40 :: nk_fname_in
   
   nev=nev_in
   nwlk=nwlk_in
   i_intf=i_intf_in
   i_dir=i_dir_in
   i_exc=i_exc_in
   xpf=xpf_in
   xmpi=xmpi_in
   xmd=xmd_in
   xmn=xmn_in
   qval=qval_in
   xA=xA_in
   i_fg=i_fg_in
   np0=np_in
   ne0=ne_in
   i_fsi=i_fsi_in

   allocate(irn(nwlk))
   irn(:)=irn_in(:)

   if(i_fg.ne.1) then

      if(1==1) then     
         open(unit=10,file='Pke_c12_sep.out',status='unknown',form='formatted')
         read(10,*) np, nE
         allocate(pv(np),Pke(np,ne),ep(nE),nk(np))

         !read the full SF if we are computing 1b cross section
         if(i_intf.ne.1) then
            do i=1,ne
               do j=1,np
                  read(10,*) pv(j), ep(i), Pke(j,i)
               enddo
               read(10,*)
            enddo

         !read only the MF SF if we are computing the 2b cross section
         else
            do i=1,ne
               do j=1,np
                  read(10,*) pv(j), ep(i),dummy, Pke(j,i)  !..the third column of the file has the full SF, the fourth is only MF, the fifth is the BG
               enddo
               read(10,*)
            enddo
         endif

         pv = pv*hbarc
         Pke=Pke/hbarc**3/(2.0d0*pi)**3
         close(10)
         he=(ep(2)-ep(1))
      endif

      if(1==0) then
         open(unit=10,file='pke12_tot.data',status='unknown',form='formatted')
         read(10,*) nE, np
         allocate(pv(np),Pke(np,ne),ep(nE),nk(np))

         do j=1,np
            read(10,*) pv(j)
            read(10,'(4(f6.1,2x,e10.3))')(ep(i),Pke(j,i),i=1,nE)
         enddo
         close(10)

         hp=pv(2)-pv(1)!pmax/dble(nbox)
         he=ep(2)-ep(1)
      endif
  endif
   
   if(i_fg.eq.1) then
      np=2*np0
      ne=1
      allocate(pv(np),ep(ne),Pke(np,ne),nk(np))
      hp=xpf/dble(np)
      he=1.0d0
      do i=1,np
         pv(i)=dble(i-0.5d0)*hp
         Pke(i,1)=1.0d0/(4.0d0*pi*xpf**3/3.0d0)
      enddo
   endif


      norm=0.0d0
      do i=1,np
         norm=norm+sum(Pke(i,:))*pv(i)**2*4.0d0*pi*(pv(2)-pv(1))*he 
      enddo
      write(6,*) 'norm',norm
      ! this needs to be updated
      Pke=Pke/norm
      do i=1,np
        nk(i)=sum(Pke(i,:))*he 
      enddo



      open(10, file='rho_0p5.dat')
      read(10,*) np_del
      allocate(pdel(np_del),pot_del(np_del))
      do i=1,np_del
         read(10,*) pdel(i),pot_del(i)
      enddo


     if(i_fsi.eq.1)then
       open(8, file='../../FSI/realOP_12C_EDAI.dat')
       read(8,*) npot
       allocate(kin(npot),pot(npot))
       do ipot=1,npot
          read(8,*)kin(ipot),pot(ipot)
          pot(ipot)=pot(ipot)
          !....kin and pot are in MeV
       enddo
    endif

      
   return
end subroutine   


subroutine mc_eval(w,resp_avg_tot,resp_err_tot)
  use mathtool
  use mympi
   implicit none
   integer*4 :: i,ie1,ie2,ne1,ne2,j
   real*8 :: w,emax,ee
   real*8 :: resp_o(2,5,nwlk)
   real*8 :: resp_avg(2,5),resp_err(2,5)
   real*8 :: resp_avg_tot(2,5),resp_err_tot(2,5)

   real*8 :: sig
   real*8 :: wmax,q2,q2max
   integer*4 :: i_acc,i_avg,i_acc_tot,i_avg_tot
   integer*4 :: ip1_o(nwlk),ie1_o(nwlk),ip2_o(nwlk),ie2_o(nwlk)
   integer*4 :: ip1_n(nwlk),ip2_n(nwlk),ie1_n(nwlk),ie2_n(nwlk)
   real*8 ::  g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk) 

   resp_avg=0.0d0
   resp_err=0.0d0
   i_acc=0
   i_avg=0
   g_o=0.0d0
   if (i_fg.eq.1) then
      xpmax=pv(np)
      emax=1.0d0
   else
      xpmax=pv(np)!-pv(1)
      emax=ep(ne)-ep(1)
   endif
   
   !
   do i=1,nwlk
      call setrn(irn(i))
      do while(g_o(i).le.0.0d0)
         ip1_o(i)=1+int(np*ran())
         ie1_o(i)=1+int(ne*ran())
         call g_eval(pv(ip1_o(i)),PkE(ip1_o(i),ie1_o(i)),xpmax,emax,norm,g_o(i))
      enddo
      call getrn(irn(i))
   enddo

   
      
   do i=1,nev
      do j=1,nwlk
         call setrn(irn(j))
         ip1_n(j)=nint(ip1_o(j)+0.05d0*np*(-1.0d0+2.0d0*ran()))
         ie1_n(j)=nint(ie1_o(j)+0.05d0*ne*(-1.0d0+2.0d0*ran()))
         if (ip1_n(j).le.np.and.ip1_n(j).ge.1.and.ie1_n(j).le.nE.and.ie1_n(j).ge.1) then
            call g_eval(pv(ip1_n(j)),PkE(ip1_n(j),ie1_n(j)), &
                 &   xpmax,emax,norm,g_n(j))
         else
            g_n(j)=0.0d0
         endif

         if (g_n(j)/g_o(j).ge.ran()) then
            ip1_o(j)=ip1_n(j)
            ie1_o(j)=ie1_n(j)
            g_o(j)=g_n(j)
            i_acc=i_acc+1
         endif
         if (i.ge.neq.and.mod(i,nvoid).eq.0) then
            if(ip1_o(j).gt.np) cycle
            call f_eval(pv(ip1_o(j)),ip1_o(j),ie1_o(j), w,resp_o(:,:,j))    
            resp_o(:,:,j)=resp_o(:,:,j)/g_o(j)
            !write(6,*)'resp_o = ', resp_o(:,:,j)
            resp_avg=resp_avg+resp_o(:,:,j)
            resp_err=resp_err+resp_o(:,:,j)**2
            i_avg=i_avg+1
         endif
         call getrn(irn(j))
      enddo
   enddo
   call addall(resp_avg,resp_avg_tot) 
   call addall(resp_err,resp_err_tot)
   call addall(i_avg,i_avg_tot) 
   call addall(i_acc,i_acc_tot) 
   if (myrank().eq.0) then
      resp_avg_tot=resp_avg_tot/dble(i_avg_tot)
      !write(6,*)'resp_avg_tot = ', resp_avg_tot
      resp_err_tot=resp_err_tot/dble(i_avg_tot)
      resp_err_tot=sqrt((resp_err_tot-resp_avg_tot**2)/dble(i_avg_tot-1))
      write(6,*)'acceptance=',dble(i_acc_tot)/dble(nev*nwlk*nproc())
   endif
 
   return
end subroutine   

subroutine f_eval(p1,ip1,ie1,w,resp)
   use mathtool
  implicit none
  integer*4 :: ip1,ie1
  real*8 :: ee,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,eef,w,q2
  real*8 :: cos_theta,jac_c,tan2
  real*8 :: v_ll,v_t
  real*8 :: tnl2,sig0,delta,rho,rhop
  real*8 :: resp(2,5)


  !p2=ran()*xpf !  to be used to mediate with a RFG on particle 2
  p2=ran()*xpmax  ! to be used to mediate with a SF on particle 2

  ctp2=-1.0d0+2.0d0*ran()
  phip2=2.0d0*pi*ran()
  phip1=2.0d0*pi*ran()

  call int_eval(p2,ctp2,phip2,p1,phip1,ip1,ie1,w,resp)
  
  return

end subroutine f_eval

subroutine int_eval(p2,ctp2,phip2,p1,phip1,ip1,ie1,w,r_now)
   use dirac_matrices         
   use mathtool
   implicit none
   real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
   real*8, parameter :: fstar=2.13d0,eps=1.0d0
   integer*4 :: ie1,ie2,ip1,ip2,i,j
   real*8 :: w,wt,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,stpp1,stp1,stp2
   real*8 :: at,bt,vt,par1,par2,pp1,den,jac,arg
   real*8 :: q2,q2e,qt2,rho,ca5,cv3,cv4,cv5,gep,gen,gmp,gmn,fpiem
   real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4)
   real*8 :: k2e_4(4),k1e_4(4), E_onshell
   real*8 :: r_cc_pi,r_cl_pi,r_ll_pi,r_t_pi,r_tp_pi
   real*8 :: r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del   
   real*8 :: r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int  
   real*8 :: dp1,dp2,delta_w
   real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2,tkin_pf,u_pq
   real*8 :: onebody(5)
   real*8 :: had_del_p_dir(4,4),had_pi_p_dir(4,4)
   real*8 :: had_del_n_dir(4,4),had_pi_n_dir(4,4)
   real*8 :: had_del_p_exc(4,4),had_pi_p_exc(4,4)
   real*8 :: had_del_n_exc(4,4),had_pi_n_exc(4,4)
   real*8 :: total_tens(4,4)
   real*8 :: dir(5),exc(5),r_now(2,5),res_p(4,4),res_n(4,4)
   real*8 :: f_1p,f_1n,f_2p,f_2n,f_A,f_P


   
   r_now=0.0d0
 
   q_4(1)=w
   q_4(2:3)=0.0d0
   q_4(4)=qval
   p1_4(1)=sqrt(p1**2+xmn**2)

   if(i_fg.eq.1) then
      q_4(1)= w - 20.0d0
   else
      q_4(1)=w-p1_4(1)-ep(ie1)+xmn  
   endif


    u_pq=0.0d0
    if(i_fsi.eq.1)then
         tkin_pf=sqrt(qval**2+xmn**2)-xmn                   
         u_pq=0.0d0
         if((tkin_pf.lt.kin(npot)).and.(tkin_pf.gt.kin(1))) call interpolint(kin,pot,npot,tkin_pf,u_pq,3)
    endif
   ctp1=((q_4(1) + p1_4(1)-u_pq)**2-p1**2-qval**2-xmn**2)/(2.0d0*p1*qval)
   q2=q_4(1)**2-qval**2

   if(abs(ctp1).gt.1.0d0) then   
      r_now=0.0d0
      return
   endif  
   

   stp1=sqrt(1.0d0-ctp1**2)
   stp2=sqrt(1.0d0-ctp2**2)


   p1_4(2)=p1*stp1*cos(phip1)
   p1_4(3)=p1*stp1*sin(phip1)
   p1_4(4)=p1*ctp1
   p2_4(1)=sqrt(p2**2+xmn**2)
   p2_4(2)=p2*stp2*cos(phip2)
   p2_4(3)=p2*stp2*sin(phip2)
   p2_4(4)=p2*ctp2

 !......define constants and ff
   !ffgnd= fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)
   cv3=fstar/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)

   cv4= -1.15/(1.0d0-q2/lsq)**2/(1.0d0-q2/4.0d0/lsq)*sqrt(3.0d0/2.0d0)

   cv5 = 0.48/(1.0d0-q2/lsq)**2/(1.0d0-q2/0.776/lsq)*sqrt(3.0d0/2.0d0)


   ca5=0.0d0!1.20d0/(1.0d0-q2/xma2)**2/(1.0d0-q2/3.0d0/xma2)*sqrt(3.0d0/2.0d0)
!   ffgnd=gep/sqrt(1.0d0-q2/(xmn+xmd)**2)/sqrt(1.0d0-q2/l3)
   !cv3 = 0.0d0 
   cv4 = 0.0d0  
   cv5 = 0.0d0
   !ca5 = 0.0d0
   rho=xpf**3/(1.5d0*pi**2)

   !Note for EM we do not use F_A or F_P, just here for consistency
   call nform(-q2/hbarc**2,f_1p,f_1n,f_2p,f_2n,gep,gen,gmp,gmn,f_A,f_P)

      
!....at this point we can define pp1_4
   pp1_4(:)=p1_4(:)+q_4(:)
!....Pauli blocking
   pp1=sqrt(sum(pp1_4(2:4)**2))
   if(pp1.lt.xpf) then   
      r_now=0.0d0
      return
   endif  


   !...define pp2
   pp2_4(:)=p2_4(:)


   !MANUALLY FIX KINEMATICS FOR TESTING
   p1_4 = (/892.28663457549249d0,-77.036892302107034d0,-26.857649369116537d0,67.942565322702237d0/)
   p2_4 = (/1002.007566138d0, 96.6849111519814d0, -193.77595264909456d0, 274.8702302143906d0/)
   pp1_4 = (/1090.4949377893222d0, 281.3919769883637d0, -74.60793248850116d0, 470.7861894114025d0/)
   pp2_4 = p2_4 
   q_4 =  (/198.20830321382988d0,358.42886929047074d0,-47.750283119384619d0 ,402.84362408870038d0/)
   w = q_4(1)
   E_onshell = sqrt(sum(p1_4(2:4)**2) + xmn**2)
   q_4(1) = w + p1_4(1) - E_onshell 
   p1_4(1) = E_onshell


   cv3 = 1.2997163009564785d0
   cv4 = 0.0d0
   cv5 = 0.0d0
   ca5 = 0.0d0
   fpiem = 0.47263563216390736d0

   !for proton knockout
   f_1p = 0.59137717615671226d0
   f_2p = 0.89838354961982203d0

   f_1n = -.017340697711910772d0
   f_2n = -0.99011610632251468d0

!...define pion momenta
   k1_4(:)=pp1_4(:)-p1_4(:)
   k2_4(:)=q_4(:)-k1_4(:)

   k1e_4(:)=pp2_4(:)-p1_4(:) 
   k2e_4(:)=pp1_4(:)-p2_4(:)

   had_pi_p_dir=0.0d0
   had_del_p_dir=0.0d0
   had_pi_n_dir=0.0d0
   had_del_n_dir=0.00d0

   had_pi_p_exc=0.0d0
   had_del_p_exc=0.00d0
   had_pi_n_exc=0.0d0
   had_del_n_exc=0.0d0
   
!.......currents
   
   if(i_intf.eq.0) then
      call current_init(w,p1_4,p2_4,pp1_4,pp2_4,q_4,k1_4,k2_4,1)
      call define_spinors()
      call det_Ja(f_1p,f_2p)
      call hadr_tens(res_p)
      call det_Ja(f_1n,f_2n)
      call hadr_tens(res_n)
   endif

   if(i_intf.eq.1) then
      !Direct contribution first
      if(i_dir.eq.1) then
         call current_init(w,p1_4,p2_4,pp1_4,pp2_4,q_4,k1_4,k2_4,1)
         call define_spinors()

         !proton knockout  
         call det_Ja(f_1p,f_2p) 
         call det_Jpi(fpiem)
         call det_JaJb_JcJd(cv3,cv4,cv5,ca5,np_del,pdel,pot_del)
         call det_J1Jdel_dir(had_del_p_dir,1)
         call det_J1Jpi_dir(had_pi_p_dir,1)
      
         !neutron knockout 
         !call det_Ja(f_1n,f_2n)
         !call det_J1Jdel_dir(had_del_n_dir,-1)
         !call det_J1Jpi_dir(had_pi_n_dir,-1)
      endif

      !Exchange contribution
      if(i_exc.eq.1) then
         call current_init(w,p1_4,p2_4,pp1_4,pp2_4,q_4,k1e_4,k2e_4,2)
         call define_spinors()

         !proton knockout  
         call det_Ja(f_1p,f_2p)
         call det_Jpi(fpiem)
         call det_JaJb_JcJd(cv3,cv4,cv5,ca5,np_del,pdel,pot_del)
         call det_J1Jdel_exc(had_del_p_exc,1)
         call det_J1Jpi_exc(had_pi_p_exc,1)

         !neutron knockout 
         !call det_Ja(f_1n,f_2n)
         !call det_J1Jdel_exc(had_del_n_exc,-1)
         !call det_J1Jpi_exc(had_pi_n_exc,-1)
      endif
   endif


   onebody(1) = res_p(1,1)+res_n(1,1)
   onebody(2) = 0.0d0
   onebody(3) = 0.0d0
   onebody(4) = res_p(2,2) + res_p(3,3) + res_n(2,2) + res_n(3,3)
   onebody(5) = 0.0d0
 
   exc(:)=0.0d0 
   dir(:)=0.0d0
   
   exc(1)=had_del_n_exc(1,1)+had_pi_n_exc(1,1)+had_del_p_exc(1,1)+had_pi_p_exc(1,1)
   exc(2)=0.0d0
   exc(3)=0.0d0
   exc(4)=had_del_n_exc(2,2)+had_del_n_exc(3,3)+ &
      &  had_pi_n_exc(2,2)+had_pi_n_exc(3,3)+ &
      &  had_del_p_exc(2,2)+had_del_p_exc(3,3)+ &
      &  had_pi_p_exc(2,2)+had_pi_p_exc(3,3) 
   exc(5)=0.0d0

   dir(1)=had_del_n_dir(1,1)+had_pi_n_dir(1,1)+had_del_p_dir(1,1)+had_pi_p_dir(1,1)!
   dir(2)=0.0d0
   dir(3)=0.0d0
   dir(4)=had_del_n_dir(2,2)+had_del_n_dir(3,3)+ &
      &  had_pi_n_dir(2,2)+had_pi_n_dir(3,3)+ &
      &  had_del_p_dir(2,2)+had_del_p_dir(3,3)+ &
      &  had_pi_p_dir(2,2)+had_pi_p_dir(3,3)  
   dir(5)=0.0d0

   !commenting out spectral function and normalization
    dp1=(4.0d0*pi*xpf**3/3.0d0)*(1.0d0/(dble(xA)/2.0d0))

  !    dp2=1 !to be used to mediate with a RFG on particle 2

    if(i_fg.eq.1) then
       dp2=1.0d0
    else
       !commenting out momentum distribution and normalization
       !all interpolint(pv,nk,np,p2,dp2,3)
       dp2=(4.0d0*pi*xpf**3/3.0d0)*(1.0d0/(dble(xA)/2.0d0))
    endif
!+
   !had_del_n_dir + had_del_p_dir  &
      !&  - had_del_n_exc - had_del_p_exc
   total_tens =  had_pi_p_dir + had_pi_n_dir - had_pi_p_exc - had_pi_n_exc + &
      &        had_del_n_dir + had_del_p_dir - had_del_n_exc - had_del_p_exc
   total_tens = (dble(xA)/rho)*total_tens*dp1*dp2/(2.0d0*pi)**6
   do i=1,4
      do j=1,4
         print*,total_tens(i,j)
      enddo
   enddo
      
     return
end subroutine   




subroutine g_eval(pj1,pke1,pmax,emax,norm,g)
    implicit none
    real*8, parameter :: pi=acos(-1.0d0)
    real*8 :: pj1,pj2,pke1,g,pmax,emax,norm
    g=(4.0d0*pi)*pj1**2*pke1
    !g=g/norm
    !g=g/norm
    return
    end subroutine


end module        
