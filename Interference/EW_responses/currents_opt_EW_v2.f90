    module dirac_matrices
    implicit none
    integer*4, private, save :: i_fl
    complex*16, private, parameter :: czero = (0.0d0,0.0d0)
    complex*16, private, parameter :: cone  = (1.0d0,0.0d0)
    complex*16, private, parameter :: ci    = (0.0d0,1.0d0)
    real*8, private, parameter :: pi=acos(-1.0d0)    
    real*8, private, parameter :: fgnd=5.0d0,fpind=0.54d0
    real*8, private, parameter :: fstar=2.13d0, xmrho=775.8d0,ga=1.26d0,fpinn2=0.081*4.0d0*pi!1.0094d0! 2.14/2.13 from JUAN, !=0.08*4.0d0*pi ARTURO
    real*8, private, parameter :: lpi=1300.0d0,lpind=1150.0d0
    real*8, private, save :: mqe, qval
    complex*16, private, save :: sig(3,2,2),id(2,2),id4(4,4),up(2),down(2)
    complex*16, private, save :: up1(2,4),up2(2,4),upp1(2,4),upp2(2,4), &
            &   ubarp1(2,4),ubarp2(2,4),ubarpp1(2,4),ubarpp2(2,4)
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4), sigma_munu(4,4,4,4)
    complex*16, private, save :: p1_sl(4,4),p2_sl(4,4),pp1_sl(4,4),pp2_sl(4,4), &
         &   k1_sl(4,4),k2_sl(4,4),q_sl(4,4), &
         &   Pi_k1(4,4),Pi_k2(4,4),Pi_k1e(4,4),Pi_k2e(4,4)
    real*8, private, save ::  p1(4),p2(4),pp1(4),pp2(4),q(4),k1(4),k2(4)
    complex*16, private, save :: J_a_mu(4,4,4),J_b_mu(4,4,4),J_c_mu(4,4,4),J_d_mu(4,4,4)
    complex*16, private, save :: Je_a_mu(4,4,4),Je_b_mu(4,4,4),Je_c_mu(4,4,4),Je_d_mu(4,4,4)
    complex*16, private, save :: J_pif(4,4,4),J_sea1(4,4,4),J_sea2(4,4,4),J_pl1(4,4,4),J_pl2(4,4,4)
    complex*16, private, save :: J_1(4,4,4)    
    real*8, private,save :: xmd,xmn,xmpi,w
contains

subroutine dirac_matrices_in(xmd_in,xmn_in,xmpi_in)
    implicit none
    integer*4 :: i,j
    real*8 :: xmd_in,xmn_in,xmpi_in
    xmd=xmd_in
    xmn=xmn_in
    xmpi=xmpi_in
    sig(:,:,:)=czero
    id(:,:)=czero
    id(1,1)=cone;id(2,2)=cone
    sig(1,1,2)=cone;sig(1,2,1)=cone
    sig(2,1,2)=-ci;sig(2,2,1)=ci
    sig(3,1,1)=cone;sig(3,2,2)=-cone
    gamma_mu=czero
    gamma_mu(1:2,1:2,1)=id;gamma_mu(3:4,3:4,1)=-id
    id4=czero    
    id4(1:2,1:2)=id;id4(3:4,3:4)=id
    do i=2,4
      gamma_mu(1:2,3:4,i)=sig(i-1,:,:)
      gamma_mu(3:4,1:2,i)=-sig(i-1,:,:)
    enddo
    gamma_mu(1:2,3:4,5)=id
    gamma_mu(3:4,1:2,5)=id
    g_munu=czero
    g_munu(1,1)=cone;g_munu(2,2)=-cone;g_munu(3,3)=-cone;g_munu(4,4)=-cone
    up(1)=cone;up(2)=czero
    down(1)=czero;down(2)=cone
    do i=1,4
       do j=1,4
          sigma_munu(:,:,i,j)=ci*0.5d0*(matmul(gamma_mu(:,:,i),gamma_mu(:,:,j)) &
               &     -matmul(gamma_mu(:,:,j),gamma_mu(:,:,i)))
       enddo
    enddo
end subroutine 


subroutine define_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigp1(2,2),sigp2(2,2),sigpp1(2,2),sigpp2(2,2)
    real*8 :: cp1,cp2,cpp1,cpp2
    sigp1=czero
    sigp2=czero
    sigpp1=czero
    sigpp2=czero
    !.....initialize quadrispinors
    up1=czero
    up2=czero
    upp1=czero
    upp2=czero
!.......initialize normalization factors
    cp1=sqrt((p1(1)+xmn)/(2.0d0*p1(1)))
    cp2=sqrt((p2(1)+xmn)/(2.0d0*p2(1)))
    cpp1=sqrt((pp1(1)+xmn)/(2.0d0*pp1(1)))
    cpp2=sqrt((pp2(1)+xmn)/(2.0d0*pp2(1)))
!.....define sigma*p
    do i=1,3
      sigp1=sigp1+sig(i,:,:)*p1(i+1)
      sigp2=sigp2+sig(i,:,:)*p2(i+1)
      sigpp1=sigpp1+sig(i,:,:)*pp1(i+1)
      sigpp2=sigpp2+sig(i,:,:)*pp2(i+1)
    enddo
!.....build quadri-spinors    
    up1(1,1:2)=up(:)
    up1(1,3:4)=matmul(sigp1(:,:),up(:))/(p1(1)+xmn)
    up1(2,1:2)=down(:)
    up1(2,3:4)=matmul(sigp1(:,:),down(:))/(p1(1)+xmn)
    up1(:,:)=cp1*up1(:,:)
!
    up2(1,1:2)=up(:)
    up2(1,3:4)=matmul(sigp2(:,:),up(:))/(p2(1)+xmn)
    up2(2,1:2)=down(:)
    up2(2,3:4)=matmul(sigp2(:,:),down(:))/(p2(1)+xmn)
    up2(:,:)=cp2*up2(:,:)
!
    upp1(1,1:2)=up(:)
    upp1(1,3:4)=matmul(sigpp1(:,:),up(:))/(pp1(1)+xmn)
    upp1(2,1:2)=down(:)
    upp1(2,3:4)=matmul(sigpp1(:,:),down(:))/(pp1(1)+xmn)
    upp1(:,:)=cpp1*upp1(:,:)
!
    upp2(1,1:2)=up(:)
    upp2(1,3:4)=matmul(sigpp2(:,:),up(:))/(pp2(1)+xmn)
    upp2(2,1:2)=down(:)
    upp2(2,3:4)=matmul(sigpp2(:,:),down(:))/(pp2(1)+xmn)
    upp2(:,:)=cpp2*upp2(:,:)
!
    ubarp1(1,1:2)=up(:)
    ubarp1(1,3:4)=-matmul(up(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(2,1:2)=down(:)
    ubarp1(2,3:4)=-matmul(down(:),sigp1(:,:))/(p1(1)+xmn)
    ubarp1(:,:)=cp1*ubarp1(:,:)
!
    ubarp2(1,1:2)=up(:)
    ubarp2(1,3:4)=-matmul(up(:),sigp2(:,:))/(p2(1)+xmn)
    ubarp2(2,1:2)=down(:)
    ubarp2(2,3:4)=-matmul(down(:),sigp2(:,:))/(p2(1)+xmn)
    ubarp2(:,:)=cp2*ubarp2(:,:)
!
    ubarpp1(1,1:2)=up(:)
    ubarpp1(1,3:4)=-matmul(up(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(2,1:2)=down(:)
    ubarpp1(2,3:4)=-matmul(down(:),sigpp1(:,:))/(pp1(1)+xmn)
    ubarpp1(:,:)=cpp1*ubarpp1(:,:)
!
    ubarpp2(1,1:2)=up(:)
    ubarpp2(1,3:4)=-matmul(up(:),sigpp2(:,:))/(pp2(1)+xmn)
    ubarpp2(2,1:2)=down(:)
    ubarpp2(2,3:4)=-matmul(down(:),sigpp2(:,:))/(pp2(1)+xmn)
    ubarpp2(:,:)=cpp2*ubarpp2(:,:)
    return
end subroutine




subroutine current_init(win,p1_in,p2_in,pp1_in,pp2_in,q_in,k1_in,k2_in,i_fl_in)
    implicit none
    integer*4 :: i,i_fl_in
    real*8 ::  p1_in(4),p2_in(4),pp1_in(4),pp2_in(4),q_in(4),k1_in(4),k2_in(4)
    real*8 :: win
    i_fl=i_fl_in
    w=win
    
    p1=p1_in
    p2=p2_in
    pp1=pp1_in
    pp2=pp2_in
    k1=k1_in
    k2=k2_in
    q=q_in
!
    p1_sl=czero
    p2_sl=czero
    pp1_sl=czero
    pp2_sl=czero
    k1_sl=czero
    k2_sl=czero
    q_sl=czero
       
    do i=1,4
       p1_sl=p1_sl+g_munu(i,i)*gamma_mu(:,:,i)*p1(i)  
       p2_sl=p2_sl+g_munu(i,i)*gamma_mu(:,:,i)*p2(i)
       pp1_sl=pp1_sl+g_munu(i,i)*gamma_mu(:,:,i)*pp1(i)  
       pp2_sl=pp2_sl+g_munu(i,i)*gamma_mu(:,:,i)*pp2(i)
       k1_sl=k1_sl+g_munu(i,i)*gamma_mu(:,:,i)*k1(i)  
       k2_sl=k2_sl+g_munu(i,i)*gamma_mu(:,:,i)*k2(i)
       q_sl=q_sl+g_munu(i,i)*gamma_mu(:,:,i)*q(i)  
    enddo
    if(i_fl.eq.1) then
       Pi_k1(:,:)=matmul(gamma_mu(:,:,5),k1_sl(:,:))/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)
       Pi_k2(:,:)=matmul(gamma_mu(:,:,5),k2_sl(:,:))/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)
    elseif(i_fl.eq.2) then
       Pi_k1e(:,:)=matmul(gamma_mu(:,:,5),k1_sl(:,:))/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)
       Pi_k2e(:,:)=matmul(gamma_mu(:,:,5),k2_sl(:,:))/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)
       !Pi_k1e(:,:)=matmul(gamma_mu(:,:,5),k1_sl(:,:))/(sum(k1(2:4)**2)+xmpi**2)
       !Pi_k2e(:,:)=matmul(gamma_mu(:,:,5),k2_sl(:,:))/(sum(k2(2:4)**2)+xmpi**2)

    endif
    return
end subroutine





subroutine det_Ja(f1v,f2v,ffa,ffp)
  implicit none
  integer*4 :: mu,nu
  real*8 :: f1v,f2v,ffa,ffp
  complex*16 :: J_1_V(4,4,4),J_1_A(4,4,4)
  
  J_1_V=czero
  J_1_A=czero

  do mu=1,4
     J_1_V(:,:,mu)=czero
     J_1_A(:,:,mu)=czero
     do nu=1,4
        J_1_V(:,:,mu)=J_1_V(:,:,mu)+ci*f2v*sigma_munu(:,:,mu,nu)&
             &     *g_munu(nu,nu)*q(nu)/2.0d0/xmn
     enddo
     J_1_V(:,:,mu)=J_1_V(:,:,mu)+f1v*gamma_mu(:,:,mu)

     !! Now add axial pieces of current 
     J_1_A(:,:,mu)=J_1_A(:,:,mu)+ffa*matmul(gamma_mu(:,:,mu),gamma_mu(:,:,5))

     J_1_A(:,:,mu)=J_1_A(:,:,mu)+ffp*gamma_mu(:,:,5)*q(mu)/xmn

  enddo

  !! Applying current conservation first to the vector current before adding axial current
  J_1_V(:,:,4) = (w/q(4))*J_1_V(:,:,1)

  J_1 = J_1_V + J_1_A


  return
end subroutine det_Ja

subroutine hadr_tens(res)
   implicit none
   integer*4 :: i1,f1,i,j
   complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4)
   complex*16:: res(4,4)

   

   do i1=1,2
      do f1=1,2
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
            J_mu_dag(f1,i1,i)=conjg(J_mu(f1,i1,i))
         enddo
      enddo
   enddo
   
   res=0.0d0
   do i1=1,2
      do f1=1,2
         do i=1,4
            do j=1,4
               res(i,j)=res(i,j)+J_mu_dag(f1,i1,i)*J_mu(f1,i1,j)
               !write(6,*) i,j, res(i,j)
             enddo  
         enddo
      enddo
   enddo

   !res(1,4)=(w/q(4))*res(1,1)
   !res(4,1)=(w/q(4))*res(1,1)
   !res(4,4)=(w/q(4))**2*res(1,1)

  
  return
end subroutine hadr_tens



subroutine det_JaJb_JcJd(cv3,cv4,cv5,ca5,np_del,pdel,pot_del)
    use mathtool
    implicit none
    integer*4 :: i,j,mu,np_del
    real*8 :: pa(4),pb(4),pc(4),pd(4),width,fpik1,fpik2,fpindk2,fpindk1
    real*8 :: ga,gb,gc,gd,cv3,cv4,cv5,ca5,pa2,pb2,pc2,pd2,pdel(np_del),pot_del(np_del)
    real*8 :: pot_pa,pot_pb,pot_pc,pot_pd
    complex*16 :: pa_sl(4,4),pb_sl(4,4),pc_sl(4,4),pd_sl(4,4)
    complex*16 :: xmd_a,xmd_b,xmd_c,xmd_d
    complex*16 :: j_a_1(4,4,4),j_a_2(4,4,4,4),RSa(4,4,4,4),RSb(4,4,4,4),j_b_1(4,4,4,4),j_b_2(4,4,4)
    complex*16 :: j_c_1(4,4,4),j_c_2(4,4,4,4),RSc(4,4,4,4),RSd(4,4,4,4),j_d_1(4,4,4,4),j_d_2(4,4,4)
    complex*16 :: J_a(4,4,4),J_b(4,4,4),J_c(4,4,4),J_d(4,4,4)
  
    if(i_fl.eq.1)then
        pa(:)=p1(:)+q(:)
        pb(:)=pp1(:)-q(:)
        pc(:)=p2(:)+q(:)
        pd(:)=pp2(:)-q(:)
    elseif(i_fl.eq.2) then
        pa(:)=p1(:)+q(:)
        pb(:)=pp2(:)-q(:)
        pc(:)=p2(:)+q(:)
        pd(:)=pp1(:)-q(:)
    endif

    pa2=pa(1)**2-sum(pa(2:4)**2)
    pb2=pb(1)**2-sum(pb(2:4)**2)
    pc2=pc(1)**2-sum(pc(2:4)**2)
    pd2=pd(1)**2-sum(pd(2:4)**2)
    
    call interpolint(pdel,pot_del,np_del,sqrt(pa2),pot_pa,1)
    call interpolint(pdel,pot_del,np_del,sqrt(pb2),pot_pb,1)
    call interpolint(pdel,pot_del,np_del,sqrt(pc2),pot_pc,1)
    call interpolint(pdel,pot_del,np_del,sqrt(pd2),pot_pd,1)
    
    pa_sl=czero
    pb_sl=czero
    pc_sl=czero
    pd_sl=czero
    call delta_se(pa(1)**2-sum(pa(2:4)**2),ga,pot_pa)
    xmd_a=xmd-0.5d0*ci*ga
    call delta_se(pb(1)**2-sum(pb(2:4)**2),gb,pot_pb)
    xmd_b=xmd-0.5d0*ci*gb
    call delta_se(pc(1)**2-sum(pc(2:4)**2),gc,pot_pc)
    xmd_c=xmd-0.5d0*ci*gc
    call delta_se(pd(1)**2-sum(pd(2:4)**2),gd,pot_pd)
    xmd_d=xmd-0.5d0*ci*gd
    fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
    fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
    fpindk1=lpind**2/(lpind**2-k1(1)**2+sum(k1(2:4)**2))
    fpindk2=lpind**2/(lpind**2-k2(1)**2+sum(k2(2:4)**2))
    do i=1,4
       pa_sl=pa_sl+g_munu(i,i)*gamma_mu(:,:,i)*pa(i)
       pb_sl=pb_sl+g_munu(i,i)*gamma_mu(:,:,i)*pb(i)
       pc_sl=pc_sl+g_munu(i,i)*gamma_mu(:,:,i)*pc(i)
       pd_sl=pd_sl+g_munu(i,i)*gamma_mu(:,:,i)*pd(i)
    enddo
! costruisco i primi due termini della corrente a 2 corpi corrispondenti ai diagrammi a,b,c e d questo e' un passaggio intermedio,
! l'espressione finale di tali correnti e' data da j_a_mu, j_b_mu, j_c_mu, j_d_mu
    do i=1,4
      j_a_1(:,:,i)=k2(i)*id4(:,:)
      j_b_2(:,:,i)=k2(i)*id4(:,:)
      j_c_1(:,:,i)=k1(i)*id4(:,:)
      j_d_2(:,:,i)=k1(i)*id4(:,:)
      !!!...I AM USING THE FULL 
      do j=1,4
         RSa(:,:,i,j)=matmul(pa_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pa(i)*pa(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pa(j)-gamma_mu(:,:,j)*pa(i))/3.0d0/xmd) &
    &    *(1.0d0/(pa(1)**2-sum(pa(2:4)**2)-xmd_a**2))
!   &     *(pa2-xmd**2)/((pa2-xmd**2)**2+xmd**2*ga**2)

         RSb(:,:,i,j)=matmul(pb_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pb(i)*pb(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pb(j)-gamma_mu(:,:,j)*pb(i))/3.0d0/xmd) &
    &    *(1.0d0/(pb(1)**2-sum(pb(2:4)**2)-xmd_b**2))
!   &     *(pb2-xmd**2)/((pb2-xmd**2)**2+xmd**2*gb**2)
         RSc(:,:,i,j)=matmul(pc_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pc(i)*pc(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pc(j)-gamma_mu(:,:,j)*pc(i))/3.0d0/xmd) &
    &    *(1.0d0/(pc(1)**2-sum(pc(2:4)**2)-xmd_c**2))
!   &     *(pc2-xmd**2)/((pc2-xmd**2)**2+xmd**2*gc**2)
         RSd(:,:,i,j)=matmul(pd_sl(:,:)+xmd*id4(:,:),g_munu(i,j)*id4(:,:)-matmul(gamma_mu(:,:,i),gamma_mu(:,:,j))/3.0d0- &
    &    2.0d0*pd(i)*pd(j)/3.0d0/xmd**2*id4(:,:)-(gamma_mu(:,:,i)*pd(j)-gamma_mu(:,:,j)*pd(i))/3.0d0/xmd) &
    &    *(1.0d0/(pd(1)**2-sum(pd(2:4)**2)-xmd_d**2))
!   &     *(pd2-xmd**2)/((pd2-xmd**2)**2+xmd**2*gd**2)
         J_a_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
         J_b_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
         J_d_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)

     !    J_a_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
     !&        (q(1)*pa(1) -q(3)*pa(3))-q(i)*pa(j))

     !    J_b_1(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
     !&        (q(1)*pb(1) -q(3)*pb(3))-q(i)*pb(j))

    !     J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
    ! &        (q(1)*pc(1)-q(3)*pc(3))-q(i)*pc(j))

    !     J_d_1(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+(cv4/xmn+cv5/xmn)*(g_munu(i,j)* &
    ! &        (q(1)*pd(1)-q(3)*pd(3))-q(i)*pd(j))
         

!         J_b_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)
!         J_c_2(:,:,i,j)=cv3*matmul(g_munu(i,j)*q_sl(:,:)-q(i)*gamma_mu(:,:,j),gamma_mu(:,:,5))+ca5*xmn*g_munu(i,j)*id4(:,:)
!         J_d_1(:,:,i,j)=cv3*matmul(gamma_mu(:,:,5),g_munu(j,i)*q_sl(:,:)-q(j)*gamma_mu(:,:,i))+ca5*xmn*g_munu(j,i)*id4(:,:)



      enddo
    enddo
! costruisco Jmua, Jmub
   do mu=1,4
      J_a(:,:,mu)=czero
      J_b(:,:,mu)=czero
      J_c(:,:,mu)=czero
      J_d(:,:,mu)=czero
      do i=1,4
         do j=1,4
            J_a(:,:,mu)=J_a(:,:,mu)+matmul(J_a_1(:,:,i)*g_munu(i,i),matmul(RSa(:,:,i,j),g_munu(j,j)*J_a_2(:,:,j,mu)))
            J_b(:,:,mu)=J_b(:,:,mu)+matmul(J_b_1(:,:,mu,i)*g_munu(i,i),matmul(RSb(:,:,i,j),g_munu(j,j)*J_b_2(:,:,j))) 
            J_c(:,:,mu)=J_c(:,:,mu)+matmul(J_c_1(:,:,i)*g_munu(i,i),matmul(RSc(:,:,i,j),g_munu(j,j)*J_c_2(:,:,j,mu)))
            J_d(:,:,mu)=J_d(:,:,mu)+matmul(J_d_1(:,:,mu,i)*g_munu(i,i),matmul(RSd(:,:,i,j),g_munu(j,j)*J_d_2(:,:,j))) 
         enddo
      enddo
   enddo
   if(i_fl.eq.1) then
      J_a_mu=J_a*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      J_b_mu=J_b*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      J_c_mu=J_c*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
      J_d_mu=J_d*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
   elseif(i_fl.eq.2)then
      Je_a_mu=J_a*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      Je_b_mu=J_b*fpik2*fpindk2*sqrt(fpinn2)*fstar/xmpi**2/xmn
      Je_c_mu=J_c*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
      Je_d_mu=J_d*fpik1*fpindk1*sqrt(fpinn2)*fstar/xmpi**2/xmn
   endif

end subroutine

subroutine det_Jpi(gep)
   implicit none
   integer*4 :: mu
   real*8 :: fpik1,fpik2,gep,frho1,frho2,fact
   fpik1=(lpi**2-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2))
   fpik2=(lpi**2-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2))
   frho1=1.0d0/(1.0d0 - (k1(1)**2-sum(k1(2:4)**2))/xmrho**2) 
   frho2=1.0d0/(1.0d0 - (k2(1)**2-sum(k2(2:4)**2))/xmrho**2)


   !...this factor is needed to fulfill current conservation, see A3 Dekker
   fact=(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)*(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & *(1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2) &
        & - 1.0d0/(k1(1)**2-sum(k1(2:4)**2)-xmpi**2)/(lpi**2-k1(1)**2+sum(k1(2:4)**2)) &
        & - 1.0d0/(k2(1)**2-sum(k2(2:4)**2)-xmpi**2)/(lpi**2-k2(1)**2+sum(k2(2:4)**2)))
   do mu=1,4
      J_pif(:,:,mu)=gep*(k1(mu)-k2(mu))*Pi_k1e(:,:)!*fact
      J_sea1(:,:,mu)=-gep*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))-frho1/ga*gamma_mu(:,:,mu)!/fpik2**2
      J_sea2(:,:,mu)=gep*matmul(gamma_mu(:,:,5),gamma_mu(:,:,mu))+frho2/ga*gamma_mu(:,:,mu)!/fpik1**2
      J_pl1(:,:,mu)=frho1/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
      J_pl2(:,:,mu)=-frho2/ga*q(mu)*q_sl(:,:)/(q(1)**2-q(4)**2-xmpi**2)
   enddo
  J_pif=J_pif*fpik1*fpik2*fpinn2/xmpi**2 
  J_sea1=J_sea1*fpik2**2*fpinn2/xmpi**2
  J_sea2=J_sea2*fpik1**2*fpinn2/xmpi**2
  J_pl1=J_pl1*fpinn2/xmpi**2*fpik2**2
  J_pl2=J_pl2*fpinn2/xmpi**2*fpik1**2
  return
end subroutine 







   subroutine det_J1Jdel_exc(res)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: Je_12c(2,2,4),Je_12d(2,2,4)   
   complex*16 :: Je_12c_dag(2,2,4),Je_12d_dag(2,2,4)

   complex*16 :: Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)   
   complex*16 :: j1ja(4,4),j1jb(4,4), j1jc(4,4),j1jd(4,4),res(4,4)
   real*8 :: res_re(4,4)

   complex*16 :: J_mu(2,2,4)   
   do i1=1,2
      do f1=1,2
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
         enddo
      enddo
   enddo



   do i1=1,2
      do f1=1,2
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         Je_1(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k1e(:,:),up1(i1,:)))
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1))


         do i=1,4
            !J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
            
            Je_12a(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_a_mu(:,:,i),up1(i1,:)))
            Je_12a_dag(f1,i1,i)=conjg(Je_12a(f1,i1,i))
            !
            Je_12b(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_b_mu(:,:,i),up1(i1,:)))
            Je_12b_dag(f1,i1,i)=conjg(Je_12b(f1,i1,i))


            Je_12c(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(Je_c_mu(:,:,i),up2(i1,:)))
            Je_12c_dag(f1,i1,i)=conjg(Je_12c(f1,i1,i))

            Je_12d(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(Je_d_mu(:,:,i),up2(i1,:)))
            Je_12d_dag(f1,i1,i)=conjg(Je_12d(f1,i1,i))


         enddo
      enddo
   enddo   
   j1ja=czero
   j1jb=czero
   j1jc=czero
   j1jd=czero
   do i=1,4
      do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  j1jb(i,j)=j1jb(i,j)+Je_12b_dag(i2,i1,i)*J_mu(f1,i1,j)*Je_2_dag(f1,i2)
                  j1jc(i,j)=j1jc(i,j)+Je_12c_dag(f1,i2,i)*J_mu(f1,i1,j)*Je_1_dag(i2,i1)
               enddo
            enddo
         enddo
      enddo
   enddo
   res=j1jb(:,:)*(8.0d0/3.0d0) +j1jc(:,:)*(8.0d0/3.0d0) 
     
   !res_re = res + conjg(res)

  ! res_re(1,4)=(w/q(4))*res_re(1,1)
  ! res_re(4,1)=(w/q(4))*res_re(1,1)
  ! res_re(4,4)=(w/q(4))**2*res_re(1,1)


   return
 end subroutine det_J1Jdel_exc



subroutine det_JdelJdel_exc(res_re)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: Je_12c(2,2,4),Je_12d(2,2,4)   
   complex*16 :: Je_12c_dag(2,2,4),Je_12d_dag(2,2,4)

   complex*16 :: Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)   
   complex*16 :: jbjb(4,4), jcjc(4,4),jbjc(4,4),res(4,4)
   real*8 :: res_re(4,4)

   complex*16 :: J_mu(2,2,4)   
   do i1=1,2
      do f1=1,2
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
         enddo
      enddo
   enddo



   do i1=1,2
      do f1=1,2
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         Je_1(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k1e(:,:),up1(i1,:)))
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1))


         do i=1,4
            !J_12a(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_a_mu(:,:,i),up1(i1,:)))
            
            Je_12a(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_a_mu(:,:,i),up1(i1,:)))
            Je_12a_dag(f1,i1,i)=conjg(Je_12a(f1,i1,i))
            !
            Je_12b(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(Je_b_mu(:,:,i),up1(i1,:)))
            Je_12b_dag(f1,i1,i)=conjg(Je_12b(f1,i1,i))


            Je_12c(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(Je_c_mu(:,:,i),up2(i1,:)))
            Je_12c_dag(f1,i1,i)=conjg(Je_12c(f1,i1,i))

            Je_12d(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(Je_d_mu(:,:,i),up2(i1,:)))
            Je_12d_dag(f1,i1,i)=conjg(Je_12d(f1,i1,i))


         enddo
      enddo
   enddo   
   jbjc=czero
   jbjb=czero
   jcjc=czero
   do i=1,4
      do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  jbjb(i,j)=jbjb(i,j)+Je_12b_dag(i2,i1,i)*Je_12b(i2,i1,j)*Je_2_dag(f1,i2)*Je_2(f1,i2)

                  jcjc(i,j)=jcjc(i,j)+Je_12c_dag(f1,i2,i)*Je_12c(f1,i2,j)*Je_1(i2,i1)*Je_1_dag(i2,i1)

                  jbjc(i,j) = jbjc(i,j) + Je_12b_dag(i2,i1,i)*Je_12c(f1,i2,j)*Je_2_dag(f1,i2)*Je_1(i2,i1)

                  
                  !j2ja(i,j)=j2ja(i,j)+Je_12a_dag(f1,i2,i)*J_mu(f1,i1,j)*Je_2_dag(i2,i1)
                  !j2jb(i,j)=j2jb(i,j)+Je_21b_dag(f1,i2,i)*J_mu(f1,i1,j)*Je_22_dag(i2,i1)
               enddo
            enddo
         enddo
      enddo
   enddo
   res=jbjb(:,:)*(4.0d0/9.0d0) +jcjc(:,:)*4.0d0/9.0d0 +jbjc(:,:)*8.0d0/9.0d0

     
   res_re = res + conjg(res)

  ! res_re(1,4)=(w/q(4))*res_re(1,1)
  ! res_re(4,1)=(w/q(4))*res_re(1,1)
  ! res_re(4,4)=(w/q(4))**2*res_re(1,1)


   return
 end subroutine det_JdelJdel_exc




    subroutine det_J1Jpi_exc(res)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)
   
   complex*16 :: Je_s2(2,2,4),Je_s2_dag(2,2,4)
   complex*16 :: Je_f(2,2,4),Je_f_dag(2,2,4),Je_s1(2,2,4),Je_s1_dag(2,2,4)
   complex*16 :: Je_p1(2,2,4),Je_p2(2,2,4), Je_p1_dag(2,2,4), Je_p2_dag(2,2,4)
   complex*16 :: j1jf(4,4),j1js(4,4),j1jp(4,4),res(4,4)   
   real*8 :: res_re(4,4)
   
   complex*16 :: J_mu(2,2,4)   
   
   
   do i1=1,2
      do f1=1,2
         Je_1(f1,i1)=sum(ubarpp2(f1,:)*matmul(Pi_k1e(:,:),up1(i1,:)))
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1))
         Je_2(f1,i1)=sum(ubarpp1(f1,:)*matmul(Pi_k2e(:,:),up2(i1,:)))
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         do i=1,4
            J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
            Je_f(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_pif(:,:,i),up1(i1,:)))
            Je_f_dag(f1,i1,i)=conjg(Je_f(f1,i1,i))
            Je_s1(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_sea1(:,:,i),up1(i1,:)))
            Je_s1_dag(f1,i1,i)=conjg(Je_s1(f1,i1,i))
            Je_s2(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_sea2(:,:,i),up2(i1,:)))
            Je_s2_dag(f1,i1,i)=conjg(Je_s2(f1,i1,i))
            Je_p1(f1,i1,i)=sum(ubarpp2(f1,:)*matmul(J_pl1(:,:,i),up1(i1,:)))
            Je_p1_dag(f1,i1,i)=conjg(Je_p1(f1,i1,i))
            Je_p2(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_pl2(:,:,i),up2(i1,:)))
            Je_p2_dag(f1,i1,i)=conjg(Je_p2(f1,i1,i))
         enddo
      enddo
   enddo   
   j1jf=czero
   j1js=czero
   j1jp=czero
   do i=1,4
      do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  j1jf(i,j)=j1jf(i,j)+Je_f_dag(i2,i1,i)*J_mu(f1,i1,j)*Je_2_dag(f1,i2)
                  j1js(i,j)=j1js(i,j)+Je_s1_dag(i2,i1,i)*J_mu(f1,i1,j)*Je_2_dag(f1,i2) &
                          &  +Je_1_dag(i2,i1)*Je_s2_dag(f1,i2,i)*J_mu(f1,i1,j)
                  j1jp(i,j)=j1jp(i,j)+Je_p1_dag(i2,i1,i)*J_mu(f1,i1,j)*Je_2_dag(f1,i2) &
                          &  +Je_1_dag(i2,i1)*Je_p2_dag(f1,i2,i)*J_mu(f1,i1,j)
               enddo
            enddo
         enddo
      enddo
   enddo
   j1jf(:,:)=j1jf(:,:)*(4.0d0)
   j1js(:,:)=j1js(:,:)*(4.0d0)
   j1jp(:,:)=j1jp(:,:)*(4.0d0)

   res=j1jf(:,:) +j1js(:,:) +j1jp(:,:)
   
     
   !res_re = res + conjg(res)

!   res_re(1,4)=(w/q(4))*res_re(1,1)
!   res_re(4,1)=(w/q(4))*res_re(1,1)
!   res_re(4,4)=(w/q(4))**2*res_re(1,1)

   return
 end subroutine det_J1Jpi_exc











   subroutine det_J1Jpi_exc_nr(gd,gev,res_re)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   real*8 :: r_cc,r_cl,r_ll,r_t,r_tp,gd,gev
   complex*16 :: J_12a(2,2,4),J_12b(2,2,4)
   complex*16 :: Je_12a(2,2,4),Je_12b(2,2,4)   
   complex*16 :: Je_12a_dag(2,2,4),Je_12b_dag(2,2,4)
   complex*16 :: Je_2(2,2),Je_2_dag(2,2)
   complex*16 :: Je_1(2,2),Je_1_dag(2,2)
   
   complex*16 :: Je_s2(2,2,4),Je_s2_dag(2,2,4)
   complex*16 :: Je_f(2,2,4),Je_f_dag(2,2,4),Je_s1(2,2,4),Je_s1_dag(2,2,4)
   complex*16 :: j1jf(4,4),j1js(4,4),res(4,4)   
   real*8 :: res_re(4,4)
   real*8 :: fpik1, fpik2
   complex*16 :: J_mu(2,2,4)   
   

   fpik1=(lpi**2-xmpi**2)/(lpi**2+sum(k1(2:4)**2))
   fpik2=(lpi**2-xmpi**2)/(lpi**2+sum(k2(2:4)**2))

   



   do i1=1,2
      do f1=1,2
         Je_2(f1,i1)=sum(sig(1:3,f1,i1)*k2(2:4))/(-sum(k2(2:4)**2)-xmpi**2)*fpik2**2*fpinn2/xmpi**2
         Je_2_dag(f1,i1)=conjg(Je_2(f1,i1))
         Je_1(f1,i1)=sum(sig(1:3,f1,i1)*k1(2:4))/(-sum(k1(2:4)**2)-xmpi**2)*fpik1**2*fpinn2/xmpi**2
         Je_1_dag(f1,i1)=conjg(Je_1(f1,i1)) 
         J_mu(f1,i1,2)=sig(2,f1,i1)*q(4)*ci*gd/2.0d0/xmn
         J_mu(f1,i1,3)=-sig(1,f1,i1)*q(4)*ci*gd/2.0d0/xmn
         do i=2,4
         !   J_mu(f1,i1,i)=sum(ubarpp1(f1,:)*matmul(J_1(:,:,i),up1(i1,:)))
            Je_s1(f1,i1,i)=gev*sig(i-1,f1,i1)            
            Je_s1_dag(f1,i1,i)=conjg(Je_s1(f1,i1,i))
            Je_s2(f1,i1,i)=gev*sig(i-1,f1,i1)
            Je_s2_dag(f1,i1,i)=conjg(Je_s2(f1,i1,i))
         enddo
      enddo
   enddo   
   j1jf=czero
   j1js=czero
   do i=2,4
      !do j=1,4
         do i1=1,2
            do f1=1,2
               do i2=1,2
                  j=i
                  j1js(i,j)=j1js(i,j)+Je_s1_dag(i2,i1,i)*J_mu(f1,i1,j)*Je_2_dag(f1,i2) &
                          &  - Je_1_dag(i2,i1)*Je_s2_dag(f1,i2,i)*J_mu(f1,i1,j)
                  enddo
            enddo
         enddo
      !enddo
   enddo
   j1jf(:,:)=j1jf(:,:)*(2.0d0)
   j1js(:,:)=j1js(:,:)*(2.0d0)

   res=j1jf(:,:) +j1js(:,:)
   !res= +j1js(:,:)
   
     
   res_re = res + conjg(res)

!   res_re(1,4)=(w/q(4))*res_re(1,1)
!   res_re(4,1)=(w/q(4))*res_re(1,1)
!   res_re(4,4)=(w/q(4))**2*res_re(1,1)

   return
 end subroutine det_J1Jpi_exc_nr




   subroutine det_J1J1_nr(gd,res_re)
   implicit none
   integer*4 :: i1,f1,i2,f2,i,j
   complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4),res(4,4)
   real*8 :: res_re(4,4),gd
   
   
   J_mu=czero
   do i1=1,2
      do f1=1,2
         J_mu(f1,i1,2)=sig(2,f1,i1)*q(4)*ci*gd/2.0d0/xmn
         J_mu(f1,i1,3)=-sig(1,f1,i1)*q(4)*ci*gd/2.0d0/xmn
      enddo
   enddo   

   J_mu_dag(:,:,:)= conjg(J_mu(:,:,:))


  res=0.0d0
   do i1=1,2
      do f1=1,2
         do i=2,4
            do j=2,4
               res(i,j)=res(i,j)+J_mu_dag(f1,i1,i)*J_mu(f1,i1,j)
               !write(6,*) i,j, res(i,j)
             enddo  
         enddo
      enddo
   enddo

   res_re = res !+ conjg(res)
   

 !  res(1,4)=(w/q(4))*res(1,1)
 !  res(4,1)=(w/q(4))*res(1,1)
 !  res(4,4)=(w/q(4))**2*res(1,1)

  
  return
end subroutine det_J1J1_nr


  
        

subroutine delta_se(pd2,width,pot)
   implicit none
   real*8 :: pd2,width,kpi,ekpi,eknuc,r2a,rfa,pot
   width=0.0d0
   if (pd2.ge.(xmpi+xmn)**2)then
      kpi=dsqrt(1.0d0/4.0d0/pd2*(pd2-(xmn+xmpi)**2)*(pd2-(xmn-xmpi)**2))
      ekpi=dsqrt(kpi**2+xmpi**2)
      eknuc=dsqrt(kpi**2+xmn**2)
      r2a= (eknuc-ekpi)**2 -4.0d0*kpi**2
      rfa=sqrt(0.95d0*xmn**2/(0.95d0*xmn**2-r2a))
!.....definizione Gamma
      width=(4.0d0*fpind)**2/12.0d0/pi/xmpi**2*kpi**3/sqrt(pd2)*(xmn+eknuc)  !&
!   &  *rfa**2!*(lpind**2/(lpind**2-xmpi**2))**2
!      width=120.0d0

!      width=0.38/(3.0d0*xmpi**2)*kpi**3/sqrt(pd2)*(xmn+eknuc)  &
!   &  *rfa**2-pot*2.0d0!*(lpind**2/(lpind**2-xmpi**2))**2


      return
   endif
!.....Eq.3.26 Phys.Rev. C 49 2650

   return
end subroutine





end module dirac_matrices





