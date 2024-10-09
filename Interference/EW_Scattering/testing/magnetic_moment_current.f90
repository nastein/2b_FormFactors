module mag_mom_current
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
    complex*16, private, save :: uk1(2,4),ukp1(2,4), &
            &   ubark1(2,4),ubarkp1(2,4)
    complex*16, private, save :: t1(2,2),t2(2,2),t1p(2,2),t2p(2,2)
    complex*16, private, save :: gamma_mu(4,4,5),g_munu(4,4), sigma_munu(4,4,4,4)
    complex*16, private, save :: p1_sl(4,4),p2_sl(4,4),pp1_sl(4,4),pp2_sl(4,4), &
         &   k1_sl(4,4),k2_sl(4,4),q_sl(4,4), &
         &   Pi_k1(4,4),Pi_k2(4,4),Pi_k1e(4,4),Pi_k2e(4,4)
    real*8, private, save ::  p1(4),p2(4),pp1(4),pp2(4),q(4),k1(4),k2(4),k(4),kp(4)
    complex*16, private, save :: J_a_mu(4,4,4),J_b_mu(4,4,4),J_c_mu(4,4,4),J_d_mu(4,4,4)
    complex*16, private, save :: Je_a_mu(4,4,4),Je_b_mu(4,4,4),Je_c_mu(4,4,4),Je_d_mu(4,4,4)
    complex*16, private, save :: J_pif(4,4,4),J_sea1(4,4,4),J_sea2(4,4,4),J_pl1(4,4,4),J_pl2(4,4,4)
    complex*16, private, save :: J_1(4,4,4)    
    real*8, private,save :: xmd,xmn,xmpi,w,xmlept,xmprobe 
contains

subroutine dirac_matrices_in(xmprobe_in,xmlept_in)
    implicit none
    integer*4 :: i,j
    real*8 :: xmprobe_in,xmlept_in
    xmlept=xmlept_in
    xmprobe=xmprobe_in
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


subroutine define_lept_spinors()
    implicit none
    integer*4 :: i
    complex*16 :: sigk1(2,2),sigkp1(2,2)
    real*8 :: ck1,ckp1
    sigk1=czero
    sigkp1=czero
    !.....initialize quadrispinors
    uk1=czero
    ukp1=czero
!.......initialize normalization factors
    ck1=sqrt((k(1)+xmprobe)/(2.0d0*k(1)))
    ckp1=sqrt((kp(1)+xmlept)/(2.0d0*kp(1)))
!.....define sigma*p
    do i=1,3
      sigk1=sigk1+sig(i,:,:)*k(i+1)
      sigkp1=sigkp1+sig(i,:,:)*kp(i+1)
    enddo
!.....build quadri-spinors    
    uk1(1,1:2)=up(:)
    uk1(1,3:4)=matmul(sigk1(:,:),up(:))/(k(1)+xmprobe)
    uk1(2,1:2)=down(:)
    uk1(2,3:4)=matmul(sigk1(:,:),down(:))/(k(1)+xmprobe)
    uk1(:,:)=ck1*uk1(:,:)
!
    ukp1(1,1:2)=up(:)
    ukp1(1,3:4)=matmul(sigkp1(:,:),up(:))/(kp(1)+xmlept)
    ukp1(2,1:2)=down(:)
    ukp1(2,3:4)=matmul(sigkp1(:,:),down(:))/(kp(1)+xmlept)
    ukp1(:,:)=ckp1*ukp1(:,:)

!
    ubark1(1,1:2)=up(:)
    ubark1(1,3:4)=-matmul(up(:),sigk1(:,:))/(k(1)+xmprobe)
    ubark1(2,1:2)=down(:)
    ubark1(2,3:4)=-matmul(down(:),sigk1(:,:))/(k(1)+xmprobe)
    ubark1(:,:)=ck1*ubark1(:,:)

    ubarkp1(1,1:2)=up(:)
    ubarkp1(1,3:4)=-matmul(up(:),sigkp1(:,:))/(kp(1)+xmlept)
    ubarkp1(2,1:2)=down(:)
    ubarkp1(2,3:4)=-matmul(down(:),sigkp1(:,:))/(kp(1)+xmlept)
    ubarkp1(:,:)=ckp1*ubarkp1(:,:)

    return
end subroutine

subroutine lepton_current_init(k_in,kp_in,q_in)
    implicit none
    real*8 :: k_in(4),kp_in(4),q_in(4)

    k=k_in  
    kp=kp_in
    q=q_in

    return
end subroutine lepton_current_init


subroutine curr(tens)
    implicit none

    integer*4 :: i1,f1,i,j
    complex*16 :: qslash(4,4)
    complex*16 :: J_mu(2,2,4),J_mu_dag(2,2,4)
    complex*16 :: tens(4,4)

    !write(6,*)'q = ', q
    qslash = momslash(q)
    !write(6,*)'qslash = ',qslash

    do i1=1,2
      do f1=1,2
         do i=1,4
            
            J_mu(f1,i1,i)=sum(ubarkp1(f1,:)*matmul(gamma_mu(:,:,i),matmul(qslash,matmul(id4(:,:)-gamma_mu(:,:,5),uk1(i1,:)))))

            J_mu(f1,i1,i)=J_mu(f1,i1,i)-sum(ubarkp1(f1,:)* &
                & matmul(qslash,matmul(gamma_mu(:,:,i),matmul(id4(:,:)-gamma_mu(:,:,5),uk1(i1,:)))))

            J_mu_dag(f1,i1,i)=conjg(J_mu(f1,i1,i))
         enddo
      enddo
   enddo
   
   tens=0.0d0
   do i1=1,2
      do f1=1,2
         do i=1,4
            do j=1,4
               tens(i,j)=tens(i,j)+J_mu_dag(f1,i1,i)*J_mu(f1,i1,j)
            enddo   
         enddo
      enddo
   enddo

   return
end subroutine curr

function momslash(mom)
    implicit none
    integer*4 :: i
    real*8 :: mom(4)
    complex*16 :: momslash(4,4)
    momslash = czero
    do i = 1,4
        momslash = momslash+g_munu(i,i)*gamma_mu(:,:,i)*mom(i) 
        !write(6,*)'momslash interior = ', momslash
    enddo
    return
end function momslash


end module mag_mom_current