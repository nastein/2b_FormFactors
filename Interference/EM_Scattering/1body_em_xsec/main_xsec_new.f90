program response_ia
  use mathtool
  use dirac_matrices
    implicit none
    real*8, parameter :: mp=938.272046d0,mn=939.56563d0,mu=931.494061d0,hbarc=197.327053d0
    real*8, parameter :: xme=0.511d0
    real*8, parameter :: small=1e-12

    integer*4 :: np0,nh0,nA,nZ,np,nc,nw,i,j,k,nskip,iform,ie,iemax,nbox,fg
    real*8, allocatable :: p0(:),dp0(:),p(:),dp(:),cost_d(:)
    real*8, allocatable :: w(:),wt(:),sigma(:),Q2(:),pke(:,:),xe(:),he
    real*8 :: ee,thetalept,coste,eef,sig_p,sig_n, cost_te,norm,mqef
    real*8 :: pmax,hp,wmax,hw,pi,rdummy,mspec,mnuc,mqe,qval,qfm,eps,sig,arg,mstar
    real*8 :: cost,cost_rel,pf,dpf,kf,TA,dwfold,ffold,hwfold,hc,delta_w,espec,epf,ep,Q2lim
    character * 40 :: nk_fname
    character*50 :: fname,en_char,theta_char,int_char,dec_char
    logical :: explicit_d

! inputs

    read(5,*) ee,thetalept
    read(5,*) nZ, nA
    read(5,*) fg
    read(5,*) iform
    read(5,*) kF
    read(5,*) nw
    read(5,*) wmax
!........how to write the output file NOT INTERESTING
    write(int_char,'(i3)') int(thetalept)
    write(dec_char,'(i1)') int(mod(thetalept*10.0d0,10.0d0))

   theta_char=trim(int_char)//'p'//trim(dec_char)
    theta_char=adjustl(theta_char)
    theta_char=trim(theta_char)
    write(en_char,'(i4)') int(ee)
    en_char=adjustl(en_char)
    en_char=trim(en_char)

    fname='C12_CBF_em_'//trim(en_char)//'_'//trim(theta_char)//'.out'
    fname=trim(fname)

! write normalized responses on file
    open(unit=14,file=fname)
    !open(unit = 27,file='normalized_pke16.data')
!!!!!!!!!!!!!!!     

    nc=100
    !Genie imposed Q2 EM limit
    !Q2lim = 20000.0d0

! initialize useful constants
    mqe=0.5d0*(mp+mn)
    mqef=mqe/hbarc
    mnuc=dble(nA)*mu
    pi=acos(-1.0d0)
    qfm=qval/hbarc
    thetalept=thetalept*pi/180.0d0
    coste= cos(thetalept)
    call dirac_matrices_in(mqef)

    !.... read the spectral function
    !...if I want to read the full table...
    if(fg.ne.1)then
       if(1==1) then
           open(unit=4,file='pke12_tot.data',status='unknown',form='formatted') 
           read(4,*) iemax, nbox
           write(6,*)'iemax = ', iemax, ', nbox = ', nbox
           allocate(p(nbox),pke(iemax,nbox),dp(nbox),xe(iemax))
           do j=1,nbox
              read(4,*) p(j)
              read(4,'(4(f6.1,2x,e10.3))')(xe(i),pke(i,j),i=1,iemax)
           enddo
           close(4)

           !do j=1,nbox
            !write(27,*)p(j)
            !write(27,'(4(f6.1,2x,e10.3))')(xe(i),pke(i,j)/2.0d0,i=1,iemax)
           !enddo
           !close(27)



           pmax=p(nbox)
           write(6,*)'pmax = ', pmax
           write(6,*)'emax = ', xe(iemax)
           hp=p(2)-p(1)!pmax/dble(nbox)
           he=xe(2)-xe(1)
           pke=pke*he
           pke=pke/dble(nZ)
        endif
        if(1==0) then
           open(unit=4,file='Pke_c12_sep.out',status='unknown',form='formatted')
           read(4,*) nbox, iemax
           allocate(p(nbox),pke(iemax,nbox),xe(iemax),dp(nbox))
           do i=1,iemax
              do j=1,nbox
                 read(4,*) p(j), xe(i), pke(i,j)
              enddo
              read(4,*)
           enddo
           p = p*hbarc
           pke=pke/hbarc**3/(2.0d0*pi)**3
           close(4)
           he=xe(2)-xe(1) 
           hp=p(2)-p(1) 
           pke=pke*he 
        endif    
    else
       !...if I want to do a Global Fermi Gas
       iemax=1
       nbox=100
       allocate(p(nbox),pke(iemax,nbox),dp(nbox))
       p(nbox)=kF
       hp=p(nbox)/dble(nbox)
       do i=1,nbox
          p(i)=dble(i)*hp
          pke(1,i)=1.0d0
       enddo
    endif
!......I build the momentum distribution and check its normalization
    do j=1,nbox
       dp(j)=sum(pke(:,j))
    enddo
    norm=sum(p(:)**2*dp(:))*4.0d0*pi*hp
    pke=pke/norm
    write(6,*) 'norm', norm
    norm=0.0d0
    do j=1,nbox
       dp(j)=sum(pke(:,j))
    enddo
    norm=sum(p(:)**2*dp(:))*4.0d0*pi*hp
    write(6,*)'n(k) norm=',sum(p(:)**2*dp(:))*4.0d0*pi*hp!


! construct omega grid and the form factors
    allocate(w(nw),wt(nw),Q2(nw))
    hw=wmax/dble(nw)
    do i=1,nw
        w(i)=(dble(i)-0.5d0)*hw
        eef = ee - w(i)
        Q2(i) = 2.0d0*ee*eef*(1.0d0 - coste)
        write(6,*)'w = ', w(i)
    enddo

! compute the response functions in the impulse approximation
    allocate(sigma(nw))
    sigma(:)=0.0d0
    !...loop on the energy transfer
    do i=1,nw
       write(6,*)'i=',i
       qval=sqrt( Q2(i)+ w(i)**2 )
       !write(6,*)'qval = ', qval
       !...loop on the energy of the SF
       do ie=1,iemax
          !...loop on the momentum of the initial nucleon
          do j=1,nbox
             ep=sqrt(p(j)**2+mqe**2)
             !.....introduce the binding
             if(fg.eq.1)then
                wt(i)=w(i)-20.0d0
             else
                wt(i)=w(i)-xe(ie)+mqe-ep
             endif
             
             !.....I solve the energy conserving delta function analytically
             ! and extract the cosine between p and q
             cost_te=((wt(i)+ep)**2-p(j)**2-qval**2-mqe**2)/(2.0d0*p(j)*qval)
             if(abs(cost_te).gt.1.0d0) cycle

             !...momentum of the final lepton
             pf=sqrt(p(j)**2+qval**2+2.0d0*qval*p(j)*cost_te)
             
             epf=sqrt(mqe**2+pf**2)
             !...Next line is for Pauli blocking
             !if(pf.ge.kf) then
                !...iform decides which form factor we want to use
                call cc1(qval/hbarc,w(i),wt(i),p(j)/hbarc,pf/hbarc,ee,thetalept,iform,sig)
                sigma(i)=sigma(i)+p(j)**2*pke(ie,j)*(dble(nZ)*sig)*epf/(p(j)*qval)
             !endif
          enddo
       enddo
       sigma(i)=sigma(i)*hp*1.e9
       write(6,*) w(i),sigma(i)
       write(14,*) w(i),sigma(i)
       flush(14)
    enddo
    close(14)
 
    end program
    
