module mc_module
   implicit none
   integer*4, private, save :: nev,xA,i_fg,i_fsi,np,ne,nwlk,npot,np_del,i_intf,i_exc,i_dir
   integer*4, private, parameter :: neq=10000,nvoid=100
   real*8, private, save ::  xpf,xpmax
   real*8, private, save:: xmpi,xmd,xmn,norm,thetalept
   complex*16, private, parameter :: ci    = (0.0d0,1.0d0), czero = (0.0d0,0.0d0)
   real*8, private, parameter :: pi=acos(-1.0d0),xme=0.0d0,xmmu=105.658357d0, ppmax=1.0d0*1.e3
   real*8, private,parameter :: hbarc=197.327053d0, G_F = 1.1664e-11*hbarc**2
   real*8, private,parameter :: cb =0.9741699d0,sb=sqrt(1.0d0-cb**2)
   real*8, private, allocatable :: pv(:),dp(:,:),ep(:),Pke(:,:),nk(:)
   real*8, private, allocatable :: kin(:),pot(:),pdel(:),pot_del(:)
   integer*8, private, allocatable, save :: irn(:)
contains

subroutine mc_init(i_intf_in,i_exc_in,i_dir_in, &
   &  i_fg_in,i_fsi_in,irn_in,nev_in,nwlk_in,xpf_in,&
   &  thetalept_in,xmpi_in,xmd_in,xmn_in,xA_in, &
   &  np_in,ne_in,nk_fname_in)

  use mathtool
   implicit none
   integer*8 :: irn_in(nwlk_in)
   integer*4 :: nev_in,nwlk_in,xA_in,i_fg_in,np_in,i,j,ne_in,np0,ne0,ien,i_intf_in
   integer*4 :: ipot,i_fsi_in,i_exc_in,i_dir_in
   real*8 :: xpf_in,xmpi_in,xmd_in,xmn_in,mlept_in,hp,he,thetalept_in,dummy
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
   thetalept=thetalept_in
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
       open(8, file='Responses/FSI/realOP_12C_EDAI.dat')
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


subroutine mc_eval(ee,w,sig_avg_tot,sig_err_tot)
  use mathtool
  use mympi
   implicit none
   integer*4 :: i,ie1,ie2,ne1,ne2,j
   real*8 :: w,emax,ee
   real*8 :: sig_o(2,nwlk)
   real*8 :: sig_avg(2),sig_err(2)
   real*8 :: sig_avg_tot(2),sig_err_tot(2)

   real*8 :: qval_in,sig
   real*8 :: wmax,q2,q2max
   integer*4 :: i_acc,i_avg,i_acc_tot,i_avg_tot
   integer*4 :: ip1_o(nwlk),ie1_o(nwlk),ip2_o(nwlk),ie2_o(nwlk)
   integer*4 :: ip1_n(nwlk),ip2_n(nwlk),ie1_n(nwlk),ie2_n(nwlk)
   real*8 ::  g_o(nwlk),g_n(nwlk),f_o(nwlk),f_n(nwlk) 

   sig_avg=0.0d0
   sig_err=0.0d0
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
            call f_eval(ee,pv(ip1_o(j)),ip1_o(j),ie1_o(j), w,sig_o(:,j))    
            sig_o(:,j)=sig_o(:,j)/g_o(j)
            sig_avg=sig_avg+sig_o(:,j)
            sig_err=sig_err+sig_o(:,j)**2
            i_avg=i_avg+1
         endif
         call getrn(irn(j))
      enddo
   enddo
   call addall(sig_avg,sig_avg_tot) 
   call addall(sig_err,sig_err_tot)
   call addall(i_avg,i_avg_tot) 
   call addall(i_acc,i_acc_tot) 
   if (myrank().eq.0) then
      sig_avg_tot=sig_avg_tot/dble(i_avg_tot)
      sig_err_tot=sig_err_tot/dble(i_avg_tot)
      sig_err_tot=sqrt((sig_err_tot-sig_avg_tot**2)/dble(i_avg_tot-1))
      write(6,*)'acceptance=',dble(i_acc_tot)/dble(nev*nwlk*nproc())
   endif
 
   return
end subroutine   

subroutine f_eval(ee,p1,ip1,ie1,w,sig)
  use mathtool
  use dirac_matrices, only: define_lept_spinors, lepton_current_init, lept_tens, contract !Don't want to accidently use saved variables
  implicit none
  integer*4 :: ip1,ie1
  real*8 :: ee,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,w,q2
  real*8 :: cos_theta,jac_c,qval,kdotq,elept,xklept
  real*8 :: tnl2,sig0,delta,rho,rhop
  real*8 :: qp(4), q(4), kprobe_4(4), klept_4(4), sig(2), re_contraction(2)
  complex*16 :: lepton_tens(4,4), r_now(2,4,4), contraction(2)

  !  p2=ran()*xpf   to be used to mediate with a RFG on particle 2
  p2=ran()*xpmax

  ctp2=-1.0d0+2.0d0*ran()
  phip2=2.0d0*pi*ran()
  phip1=2.0d0*pi*ran()
  sig=0.0d0
  elept = ee - w
  if(elept.lt.0.0d0) stop
  xklept = sqrt(elept**2 - xmmu**2)
  cos_theta=cos(thetalept)

  sig0 = (G_F*cb/hbarc)**2 /(2.0d0*pi)*xklept*elept

  q2=2.0d0*ee*elept - 2.0d0*ee*xklept*cos_theta - xmmu**2

  if(q2.lt.0.0d0) stop
  qval=sqrt(q2+w**2)
  kdotq = xmmu**2/2.0d0 + ee*w - (w**2 - qval**2)/2.0d0

  qp(1)=ee+elept
  qp(4)=2.0d0*kdotq/qval - qval
  qp(3)=0.0d0 
  qp(2)=2.0d0*sqrt(ee**2 - kdotq**2/qval**2)

  q(1)=w
  q(2)=0.0d0 
  q(3)=0.0d0 
  q(4)=qval

  kprobe_4=(q+qp)/2.0d0 
  klept_4=(qp-q)/2.0d0

  call lepton_current_init(kprobe_4, klept_4)
  call define_lept_spinors()
  call lept_tens(lepton_tens)

  call int_eval(p2,ctp2,phip2,p1,phip1,ip1,ie1,w,qval,r_now)

  !Contract the lepton and 1b hadron/12b hadron tensors separately
  contraction(1) = contract(r_now(1,:,:),lepton_tens)
  contraction(2) = contract(r_now(2,:,:),lepton_tens)


  re_contraction(1) = contraction(1)
  !r12munu Lmunu + r12^* munu Lmunu
  re_contraction(2) = contraction(2) + conjg(contraction(2))

  sig(:)=sig0*(re_contraction(:))*1.e15 !1e-15 fm^2
  return

end subroutine f_eval

subroutine int_eval(p2,ctp2,phip2,p1,phip1,ip1,ie1,w,qval,r_now)
   use dirac_matrices         
   use mathtool
   implicit none
   real*8, parameter :: lsq=0.71*1.e6,l3=3.5d0*1.e6,xma2=1.1025d0*1.e6
   real*8, parameter :: fstar=2.13d0,eps=1.0d0
   integer*4 :: ie1,ie2,ip1,ip2
   real*8 :: w,ctpp1,p2,ctp2,phip2,p1,ctp1,phip1,stpp1,stp1,stp2
   real*8 :: at,bt,vt,par1,par2,pp1,den,jac,arg
   real*8 :: q2,q2e,rho,ca5,cv3,cv4,cv5,gep,gen,gmp,gmn,qval
   real*8 :: p1_4(4),p2_4(4),pp1_4(4),pp2_4(4),k2_4(4),k1_4(4),q_4(4),pp_4(4),q_4e(4)
   real*8 :: k2e_4(4),k1e_4(4)
   real*8 :: r_cc_pi,r_cl_pi,r_ll_pi,r_t_pi,r_tp_pi
   real*8 :: r_cc_del,r_cl_del,r_ll_del,r_t_del,r_tp_del   
   real*8 :: r_cc_int,r_cl_int,r_ll_int,r_t_int,r_tp_int  
   real*8 :: dp1,dp2,delta_w
   real*8 :: tkin_pp1,tkin_pp2, u_pp1,u_pp2,tkin_pf,u_pq
   complex*16 :: had_del_dir(4,4),had_pi_dir(4,4)
   complex*16 :: had_del_exc(4,4),had_pi_exc(4,4),had_intf(4,4)
   complex*16 :: onebody(4,4)
   complex*16 :: r_now(2,4,4)
   real*8 :: f_1p,f_1n,f_2p,f_2n,f_A,f_P,ff1v,ff2v


   
   r_now=0.0d0
 
   q_4(1)=w
   q_4(2:3)=0.0d0
   q_4(4)=qval
   q_4e = q_4
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
   q2e = q_4e(1)**2 - qval**2

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

   ca5=1.20d0/(1.0d0-q2/xma2)**2/(1.0d0-q2/3.0d0/xma2)*sqrt(3.0d0/2.0d0)
!   ffgnd=gep/sqrt(1.0d0-q2/(xmn+xmd)**2)/sqrt(1.0d0-q2/l3)

   !cv3 = 0.0d0 
   !cv4 = 0.0d0  
   !cv5 = 0.0d0
   !ca5 = 0.0d0
   rho=xpf**3/(1.5d0*pi**2)



   call nform(-q2/hbarc**2,f_1p,f_1n,f_2p,f_2n,gep,gen,gmp,gmn,f_A,f_P)
   ff1v=f_1p-f_1n
   ff2v=f_2p-f_2n
   !Setting pseudoscalar to 0 for now
   f_P = 0.0d0
   print*,ff1v,ff2v,f_A
      
!....at this point we can define pp1_4
   pp1_4(:)=p1_4(:)+q_4(:)
!....Pauli blocking
   pp1=sqrt(sum(pp1_4(2:4)**2))
   
   !!!!!!!!!!!!!!!!!!!!!!
   !Commenting out pauli blocking for a moment
   if(pp1.lt.xpf) then   
      r_now=0.0d0
      return
   endif  

   !...define pp2
   pp2_4(:)=p2_4(:)
   
!...define pion momenta
   k1_4(:)=pp1_4(:)-p1_4(:)
   k2_4(:)=q_4(:)-k1_4(:)
   
   k1e_4(:)=pp2_4(:)-p1_4(:) 
   k2e_4(:)= pp1_4(:)-p2_4(:)


   had_pi_dir=0.0d0
   had_del_dir=0.0d0
   had_pi_exc=0.0d0
   had_del_exc=0.0d0
   
!.......currents
   !direct
   call current_init(w,p1_4,p2_4,pp1_4,pp2_4,q_4,k1_4,k2_4,1)
   call define_spinors()
   call det_Ja(ff1v,ff2v,f_A,f_P)
   
   if(i_intf.eq.0) then
      call hadr_tens(onebody)
   endif

   if(i_intf.eq.1) then 
      if(i_dir.eq.1) then   
         call det_Jpi(gep-gen)
         call det_JaJb_JcJd(cv3,cv4,cv5,ca5,np_del,pdel,pot_del)
         call det_J1Jdel_dir(had_del_dir,-1)
         call det_J1Jpi_dir(had_pi_dir,-1)
      endif
      if(i_exc.eq.1) then
         call current_init(w,p1_4,p2_4,pp1_4,pp2_4,q_4,k1e_4,k2e_4,2)
         call define_spinors()
         call det_Ja(ff1v,ff2v,f_A,f_P)
         call det_Jpi(gep-gen)
         call det_JaJb_JcJd(cv3,cv4,cv5,ca5,np_del,pdel,pot_del)
         call det_J1Jdel_exc(had_del_exc,-1)
         call det_J1Jpi_exc(had_pi_exc,-1) 
      endif
   endif

   !Dir - Exchange
   had_intf = had_del_dir + had_pi_dir - had_del_exc - had_pi_exc

    dp1=PkE(ip1,ie1)*(4.0d0*pi*xpf**3/3.0d0)*(norm/(dble(xA)/2.0d0))

  !    dp2=1 !to be used to mediate with a RFG on particle 2

    if(i_fg.eq.1) then
       dp2=1.0d0
    else
       call interpolint(pv,nk,np,p2,dp2,3)
       dp2=dp2*(4.0d0*pi*xpf**3/3.0d0)*(norm/(dble(xA)/2.0d0))
    endif


      !ONEBODY
      if(i_intf.ne.1) then
         r_now(1,:,:) = (dble(xA)/rho)*dp1*p1**2*(2.0d0*pi)*onebody(:,:)*pp1_4(1)  &  
         &     /(p1*qval)/(2.0d0*pi)**3!/(norm/(dble(xA)/2.0d0)) !*dp2*p2**2*(4.0d0*pi)/(4.0d0*pi*xpf**3/3.0d0)  
         ! note that since I am using the MF SF I can't compute the 1b contribution in a correct way

         r_now(2,:,:) = 0.0d0

      !TWOBODY
      else
         r_now(1,:,:) = 0.0d0

         r_now(2,:,:) = (dble(xA)/rho)*dp1*p1**2*(2.0d0*pi)*(had_intf(:,:))*pp1_4(1)  &
         &     /(p1*qval)/(2.0d0*pi)**3 *dp2*p2**2*(4.0d0*pi*xpmax)/(2.0d0*pi)**3 ! use xpf  to mediate with a RFG on particle 2
      endif
      
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
