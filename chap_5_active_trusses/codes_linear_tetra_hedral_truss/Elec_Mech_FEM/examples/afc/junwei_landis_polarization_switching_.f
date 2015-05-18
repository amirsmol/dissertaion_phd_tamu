      module polarization_switching_

        real(8),parameter::tolrtao=1e-4
        real(8),parameter::tolrphi=1e-4
        real(8),parameter::tol01=1e-30

        real(8),parameter::YONG=61.0e09;
        real(8),parameter::POISS=0.31
        real(8),parameter::DIELC=2.77e-8

        real(8),parameter::DPZ33=5.93e-10;
        real(8),parameter::DPZ31=-2.74e-10;
        real(8),parameter::DPZ15=7.41e-10;


        real(8),parameter::P0=0.26;
        real(8),parameter::E0=8.0e5
        real(8),parameter::HE0 =7.14e5

        real(8),parameter::AME=1.4

      contains

      subroutine su3dsp(almda,tao,taorm,sigma)
c*******************************************************************************
c     given:  taorm, tao(1~6), sigma(7~9)
c     update: taorm, tao(7~9), sigma(1~6)
c*******************************************************************************
      implicit double precision (a-h,o-z)
      parameter (itmxmum=500)
      dimension
     1 tao(9), !input
     2 taorm(9), !input, output
     2 sigma(9), !output
     3 ptaorm(9),pr(3),gmtx(9,9),gpmtx(9,9),ginv(9,9),sigmabr(9),
     4 umtx(9,9),amtx(9,9),sigmabk(9),sigmaht(9),hmtx(9,9),anvct(9),
     5 amvct(9),vmtx(9,9),wmtx(9,9),rtao(9),bmtx(9,9),ajmtx(9,9),
     6 ajinv(9,9),deta(9,9),dtaorm(9),un(3) !local
!      common /matconst/ yong,poiss,dielc,dpz33,dpz 31,dpz15,p0,e0,he0,ame

      common /io/ it,in

c     initialization
      do i=1,9
      do j=1,9
      if(i.eq.j)then
              deta(i,j)=1.0d0
      else
              deta(i,j)=0.0d0
      endif
      end do
      end do

c     save taorm of the previous converged step
      do 10 i=1,9
      ptaorm(i)=taorm(i)
10    continue

c     ****** predictor step: compute the tiral state ******
      do i=1,3
      pr(i)=taorm(i+6)
      enddo
c     solve for sigma(1~6),tao(7~9)
      call gmat(pr,
     1          ginv,gpmtx,gmtx)
c
      do i=1,6
        sigma(i)=0.0d0
      do i1=1,6
        sigma(i)=sigma(i)+gpmtx(i,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        sigma(i)=sigma(i)+gpmtx(i,i1+6)*sigma(i1+6)
      end do
      end do
c
      do i=1,3
        tao(i+6)=taorm(i+6)
      do i1=1,6
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1+6)*sigma(i1+6)
      end do
      end do
c     sigmaht, phi
      call dgmat(pr,sigma,
     1           sigmabr,umtx,amtx)
      call flhdelc(pr,sigma,sigmabr,
     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
c     check consistency condition for tiral step
      if(phi.lt.0.0d0) then
              return
      endif

c     ****** corrector step: solve the residual equations ******
c     set guess vale of pr when entering the corrector step for the first time
      prmnt=0.0d0
      do i=1,3
          prmnt=prmnt+pr(i)*pr(i)
      end do
      prmnt=dsqrt(prmnt)
      if(prmnt.lt.tol01)then
          beta=1.0d0/dielc
          epcnst=beta*he0/(beta+he0)
          stsmnt=0.0d0
          do i=1,3
              stsmnt=stsmnt+sigma(i+6)*sigma(i+6)
          end do
          stsmnt=dsqrt(stsmnt)
          dstsmnt=stsmnt-e0
          do i=1,3
              un(i)=sigma(i+6)/stsmnt
          end do
          edisp=dstsmnt/epcnst
          edispe=dstsmnt/beta
          prmnt=edisp-edispe
          do i=1,3
              pr(i)=prmnt*un(i)
          end do
          do i=1,3
              taorm(i+6)=pr(i)
          end do
      end if
      call dgmat(pr,sigma,
     1           sigmabr,umtx,amtx)
      call flhdelc(pr,sigma,sigmabr,
     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
c
c___________________________________________________
c___________________________________________________
      i100=1
c___________________________________________________
c___________________________________________________
      iter=0
500   iter=iter+1
c      write(*,*)iter,taorm(9),tao(9)-taorm(9),sigma(9)
c      pause
c     initialization
      do i=1,9
      rtao(i)=0.0d0
      dtaorm(i)=0.0d0
      do j=1,9
      ajinv(i,j)=0.0d0
      end do
      end do
c     compute residuals
      do 30 i=1,9
      rtao(i)=rtao(i)-taorm(i)+ptaorm(i)+almda*anvct(i)
30    continue
      errphi=dabs(phi)/e0/e0
c     check for convergence
      anmstnrm=0.0d0
      anmprm=0.0d0
      do 40 i=1,6
      anmstnrm=anmstnrm+ptaorm(i)*ptaorm(i)
40    continue
      do 50 i=1,3
      anmprm=anmprm+ptaorm(i+6)*ptaorm(i+6)
50    continue
      anorm1=0.0d0
      anorm2=0.0d0
      do 60 i=1,6
      anorm1=anorm1+rtao(i)*rtao(i)
60    continue
      if (anmstnrm.gt.tol01)then
      anorm1=anorm1/anmstnrm
      endif
      do 70 i=1,3
      anorm2=anorm2+rtao(i+6)*rtao(i+6)
70    continue
      if(anmprm.gt.tol01)then
      anorm2=anorm2/anmprm
      endif
      errtao=anorm1+anorm2
      errtao=dsqrt(errtao)
c___________________________________________________
c___________________________________________________
      if(iter.eq.i100)then
          write(25,2010)errtao,errphi
          i100=i100+1
      end if
c___________________________________________________
c___________________________________________________
      if(errtao.lt.tolrtao .and.
     1   errphi.lt.tolrphi)then
              write(25,1030)iter-1,errtao,errphi
              return
      endif
      if(iter.gt.itmxmum)then
              write(25,1010)itmxmum
              write(25,1020)errtao,errphi
              write(25,*)sigma(9)
              stop
              return
      endif
c     compute inverted reduced jacobian
      call bmat(2,amtx,umtx,gmtx,hmtx,
     1          bmtx)
      do 80 i=1,9
      do 80 j=1,9
      ajinv(i,j)=ajinv(i,j)+deta(i,j)-almda*wmtx(i,j)
      do 80 ii=1,9
      ajinv(i,j)=ajinv(i,j)+almda*vmtx(i,ii)*bmtx(ii,j)
80    continue
      call mtxinv(9,ajinv,ajmtx)
c     calculate the increment to consistent parameter
      anumer=0.0d0
      adenom=0.0d0
      anumer=anumer+phi
      do 90 i=1,9
      do 90 j=1,9
      anumer=anumer+amvct(i)*ajmtx(i,j)*rtao(j)
      adenom=adenom-amvct(i)*ajmtx(i,j)*anvct(j)
      do 90 ii=1,9
      anumer=anumer-anvct(i)*bmtx(i,j)*ajmtx(j,ii)*rtao(ii)
      adenom=adenom+anvct(i)*bmtx(i,j)*ajmtx(j,ii)*anvct(ii)
90    continue
      dlmda=anumer/adenom
c     update variables
      almda=almda+dlmda
      do i=1,9
      do j=1,9
          dtaorm(i)=dtaorm(i)+ajmtx(i,j)*(rtao(j)+anvct(j)*dlmda)
      end do
          taorm(i)=taorm(i)+dtaorm(i)
      end do
      do i=1,3
      pr(i)=taorm(i+6)
      end do
c     update sigma(1~6),tao(7~9)
      call gmat(pr,
     1          ginv,gpmtx,gmtx)
c
      do i=1,6
        sigma(i)=0.0d0
      do i1=1,6
        sigma(i)=sigma(i)+gpmtx(i,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        sigma(i)=sigma(i)+gpmtx(i,i1+6)*sigma(i1+6)
      end do
      end do
c
      do i=1,3
        tao(i+6)=taorm(i+6)
      do i1=1,6
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1+6)*sigma(i1+6)
      end do
      end do
c     update sigmaht,phi
      call dgmat(pr,sigma,
     1           sigmabr,umtx,amtx)
      call flhdelc(pr,sigma,sigmabr,
     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
c
      goto 500

1010  format(/,8x,'***** convergence criterion is not satisfied
     1 after muximum #:',i9,'   of iteration allowed is reached*****')
1020  format(/,8x,'***** error after itmxmum iteration is:  errtao:',
     1 e12.4,'  errphi:',e12.4)
1030  format(/,8x,'***** converged!!!   error after',i9,
     1'   iterations is,  errtao:',e12.4,'  errphi:',e12.4)
2010  format(3x,2e12.4)

      end subroutine su3dsp









      subroutine bmat(itype,amtx,umtx,gmtx,hmtx,
     1                bmtx)
c    *****************************************************
c     itype=1: vector potential, for state update
c     itype=2: scalor potential, for state update
c     itype=3: stress input    , for state update
c     itype=4:                 , for consistent matrix
c    *****************************************************
      implicit double precision (a-h,o-z)
      dimension
     1 amtx(9,9),umtx(9,9),gmtx(9,9),hmtx(9,9), ! input
     2 bmtx(9,9), ! output
     3 deta(9,9),fmtx(9,9),omtx(9,9),g22mtx(3,3),g22inv(3,3) ! local
      data
     1    r0   ,rp5  ,r1   ,r2   ,r3   ,i3/
     2    0.0d0,0.5d0,1.0d0,2.0d0,3.0d0,3 /

c******************initialization********************
      do 10 i=1,9
      do 10 j=1,9
      bmtx(i,j)=r0
      fmtx(i,j)=r0
      omtx(i,j)=r0
10    continue

c*********** parameters deta(i,j)**********************
      do i=1,9
      do j=1,9
       if(i.eq.j)then
        deta(i,j)=r1
       else
        deta(i,j)=r0
       endif
      enddo
      enddo

c************** bmtx *****************
      do 30 i=1,9
      do 30 j=1,9
      bmtx(i,j)=bmtx(i,j)+hmtx(i,j)-umtx(i,j)
c
      if(itype.eq.1 .or. itype.eq.4)then
c
      do 20 ii=1,9
      do 20 jj=1,9
      bmtx(i,j)=bmtx(i,j)+(deta(i,ii)+amtx(ii,i))*gmtx(ii,jj)
     1          *(deta(jj,j)+amtx(jj,j))
20    continue
c
      elseif(itype.eq.2)then
c
      do 21 i1=1,9
      do 21 i2=1,9
        fmtx(i1,i2)=r0
      do 21 i3=1,9
        fmtx(i1,i2)=fmtx(i1,i2)+gmtx(i1,i3)*(deta(i3,i2)+amtx(i3,i2))
21    continue
      do 22 i1=1,3
      do 22 i2=1,3
        g22mtx(i1,i2)=gmtx(i1+6,i2+6)
22    continue
      call mtxinv(3,g22mtx,g22inv)
      do 23 i1=1,6
      do 23 i2=1,6
        omtx(i1,i2)=fmtx(i1,i2)
      do 23 i3=1,3
      do 23 i4=1,3
        omtx(i1,i2)=omtx(i1,i2)-gmtx(i1,i3+6)*g22inv(i3,i4)*
     1                                        fmtx(i4+6,i2)
23    continue
      do 24 i1=1,6
      do 24 i2=1,3
        omtx(i1,i2+6)=fmtx(i1,i2+6)
      do 24 i3=1,3
      do 24 i4=1,3
        omtx(i1,i2+6)=omtx(i1,i2+6)-gmtx(i1,i3+6)*g22inv(i3,i4)*
     1                                           fmtx(i4+6,i2+6)
24    continue
      do 25 i1=1,9
        bmtx(i,j)=bmtx(i,j)+(deta(i,i1)+amtx(i1,i))*omtx(i1,j)
25    continue
c
      else
c     do nothing
      end if
30    continue

      return
      end subroutine bmat




      subroutine dgmat(pr,sigma,
     1                 sigmabr,umtx,amtx)
      implicit double precision (a-h,o-z)
      dimension
     1  pr(3),sigma(9),   !input
     2  sigmabr(9),umtx(9,9),amtx(9,9),   !!! output
     3  dginv(9,9,9),ddginv(9,9,9,9),un(3),deta(3,3),alf(3,3) !local
!      common /matconst/ yong,poiss,dielc,dpz33,dpz31,dpz15,p0,e0,he0,ame
      common /errtol/ tolrtao,tolrphi,toldlmda,toldtaorm,tol01
      data
     1    r0   ,rp5  ,r1   ,r2   ,r3   ,i3/
     2    0.0d0,0.5d0,1.0d0,2.0d0,3.0d0,3 /

c******************initialization********************
      do 10 i=1,9
      sigmabr(i)=r0
      do 10 j=1,9
      umtx(i,j)=r0
      amtx(i,j)=r0
      do 10 ii=1,9
      dginv(i,j,ii)=r0
      do 10 jj=1,9
      ddginv(i,j,ii,jj)=r0
10    continue
      do i=1,3
      un(i)=r0
      do j=1,3
      alf(i,j)=r0
      if(i.eq.j)then
      deta(i,j)=r1
      else
      deta(i,j)=r0
      endif
      enddo
      enddo

c*********** parameters **********************
      prmnt=r0
      do 20 i=1,3
      prmnt=prmnt+pr(i)*pr(i)
20    continue
      if(prmnt.gt.tol01)then
      prmnt=sqrt(prmnt)
      do 30 i=1,3
      un(i)=pr(i)/prmnt
30    continue
      do 40 i=1,3
      do 40 j=1,3
      alf(i,j)=deta(i,j)-un(i)*un(j)
40    continue
      end if

c ******************** dginv ****************
      if(prmnt.gt.tol01)then
      do 80 i=1,3
      do 80 j=1,3
      do 80 k=1,6
      if(k.le.3)then
      ki=k
      kj=k
      elseif(k.eq.4)then
      ki=1
      kj=2
      elseif(k.eq.5)then
      ki=1
      kj=3
      else
      ki=2
      kj=3
      end if
      dginv(i+6,k,j+6)
     *  =dpz33/p0*(
     1   un(i)*un(ki)*un(kj)*un(j)+
     2   un(ki)*un(kj)*alf(i,j)+
     3   un(i)*alf(ki,j)*un(kj)+
     4   un(i)*un(ki)*alf(kj,j))
     *  +dpz31/p0*(
     1   alf(ki,kj)*alf(i,j)+
     2   un(i)*un(j)*alf(ki,kj)-
     3   un(i)*un(ki)*alf(kj,j)-
     4   un(i)*alf(ki,j)*un(kj))
     *  +dpz15/p0/r2*(
     1   alf(ki,j)*alf(kj,i)+
     2   alf(ki,i)*alf(kj,j)+
     3   un(ki)*un(j)*alf(kj,i)+
     4   un(j)*alf(ki,i)*un(kj)-
     5   un(i)*un(ki)*alf(kj,j)-
     6   un(i)*alf(ki,j)*un(kj)-
     7   un(ki)*un(kj)*alf(i,j)*r2)
80    continue
      do 90 i=1,6
      do 90 j=1,3
      do 90 ii=1,3
        dginv(i,j+6,ii+6)=dginv(j+6,i,ii+6)
90    continue
      end if

c**************** ddginv ****************
      if(prmnt.gt.tol01)then
      do 100 i=1,3
      do 100 ii=1,3
      do 100 jj=1,3
      do 100 k=1,6
      if(k.le.3)then
      ki=k
      kj=k
      elseif(k.eq.4)then
      ki=1
      kj=2
      elseif(k.eq.5)then
      ki=1
      kj=3
      else
      ki=2
      kj=3
      end if
      ddginv(i+6,k,ii+6,jj+6)
     *  =dpz33/p0/prmnt*(
     1   un(ki)*alf(i,ii)*alf(kj,jj)+
     2   un(ki)*alf(kj,ii)*alf(i,jj)+
     3   alf(ki,ii)*un(kj)*alf(i,jj)+
     4   alf(ki,jj)*un(kj)*alf(i,ii)+
     5   un(i)*alf(ki,ii)*alf(kj,jj)+
     6   un(i)*alf(ki,jj)*alf(kj,ii)-
     7   r2*un(i)*un(ki)*un(kj)*alf(ii,jj))
     *  +dpz31/p0/prmnt*(
     1   r2*un(i)*un(ki)*un(kj)*alf(ii,jj)-
     2   un(ki)*alf(i,ii)*alf(kj,jj)-
     3   un(ki)*alf(kj,ii)*alf(i,jj)-
     4   alf(ki,ii)*un(kj)*alf(i,jj)-
     5   alf(ki,jj)*un(kj)*alf(i,ii)-
     6   un(i)*alf(ki,ii)*alf(kj,jj)-
     7   un(i)*alf(ki,jj)*alf(kj,ii))
     *  +dpz15/p0/prmnt*(
     1   r2*un(i)*un(ki)*un(kj)*alf(ii,jj)-
     2   un(ki)*alf(i,ii)*alf(kj,jj)-
     3   un(ki)*alf(kj,ii)*alf(i,jj)-
     4   alf(ki,ii)*un(kj)*alf(i,jj)-
     5   alf(ki,jj)*un(kj)*alf(i,ii)-
     6   un(i)*alf(ki,ii)*alf(kj,jj)-
     7   un(i)*alf(ki,jj)*alf(kj,ii))
100   continue
      do 110 i=1,6
      do 110 j=1,3
      do 110 ii=1,3
      do 110 jj=1,3
      ddginv(i,j+6,ii+6,jj+6)=ddginv(j+6,i,ii+6,jj+6)
110   continue
      end if

c    ********sigmabr, amtx and umtx****************
      do 120 i=1,9
      do 120 j=1,9
      do 120 ii=1,9
      sigmabr(i)=sigmabr(i)-r1/r2*sigma(j)*sigma(ii)*dginv(j,ii,i)
      amtx(i,j)=amtx(i,j)+sigma(ii)*dginv(ii,i,j)
      do 120 jj=1,9
      umtx(i,j)=umtx(i,j)+r1/r2*sigma(ii)*sigma(jj)*ddginv(ii,jj,i,j)
120   continue

      return
      end subroutine dgmat







      subroutine flhdelc(pr,sigma,sigmabr,
     1                   sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
      implicit double precision (a-h,o-z)
      dimension
     1  pr(3),sigma(9),sigmabr(9), ! input
     2  sigmabk(9),sigmaht(9),hmtx(9,9),anvct(9),amvct(9),vmtx(9,9),
     *  wmtx(9,9), !!! output
     3  un(3),deta(3,3),alf(3,3) ! local
!      common /matconst/ yong,poiss,dielc,dpz33,dpz31,dpz15,p0,e0,he0,ame
      common /errtol/ tolrtao,tolrphi,toldlmda,toldtaorm,tol01
      data
     1    r0   ,rp5  ,r1   ,r2   ,r3   ,i3/
     2    0.0d0,0.5d0,1.0d0,2.0d0,3.0d0,3 /

c       ***initialization********************
      do 10 i=1,9
        sigmabk(i)=r0
        sigmaht(i)=r0
        anvct(i)=r0
        amvct(i)=r0
      do 10 j=1,9
        hmtx(i,j)=r0
        vmtx(i,j)=r0
        wmtx(i,j)=r0
10    continue
      do 20 i=1,3
       un(i)=r0
      do 20 j=1,3
       alf(i,j)=r0
       if(i.eq.j)then
        deta(i,j)=r1
       else
        deta(i,j)=r0
       endif
20    continue

c  *********parameters************
      prmnt=r0
      do 30 i=1,3
      prmnt=prmnt+pr(i)**2
30    continue
      if(prmnt.gt.tol01)then
      prmnt=sqrt(prmnt)
      do 40 i=1,3
      un(i)=pr(i)/prmnt
40    continue
      do 50 i=1,3
      do 50 j=1,3
      alf(i,j)=deta(i,j)-un(i)*un(j)
50    continue
      end if

c     ******sigmabk and hardening modulus**********
      if(prmnt.gt.tol01)then
      do 60 i=1,3
       sigmabk(i+6)=he0*p0/(ame-r1)*((r1-prmnt/p0)**(r1-ame)-r1)*un(i)
      do 60 j=1,3
       hmtx(i+6,j+6)=he0*(r1-prmnt/p0)**(-ame)*un(i)*un(j)+p0/prmnt*he0
     1 /(ame-r1)*((r1-prmnt/p0)**(r1-ame)-r1)*alf(i,j)
60    continue
c***********************************************************
c      else
c       hmtx(9,9)=he0*(r1-prmnt/p0)**(-ame)
c***********************************************************
      end if

c     ******sigmaht*********************
      do 70 i=1,3
       sigmaht(i+6)=sigma(i+6)-sigmabr(i+6)-sigmabk(i+6)
70    continue

c     ******anvct and vmtx, amvct=r0, wmtx=r0**********
      do 80 i=1,3
       anvct(i+6)=r2*sigmaht(i+6)
      do 80 j=1,3
       vmtx(i+6,j+6)=r2*deta(i,j)
80    continue

c     ****** phi **********
      phi=r0
      do 90 i=1,3
       phi=phi+sigmaht(i+6)*sigmaht(i+6)
90    continue
      phi=phi-e0*e0
c
      return
      end subroutine flhdelc











      subroutine gmat(pr,
     1                ginv,gpmtx,gmtx)
      implicit double precision (a-h,o-z)
      dimension
     1  pr(3),   !input
     2  ginv(9,9),gpmtx(9,9),gmtx(9,9), !!! output
     3  un(3),deta(3,3),alf(3,3),temp1(3,3),temp2(3,3)  !local
!      common /matconst/ yong,poiss,dielc,dpz33,dpz31,dpz15,p0,e0,he0,ame
      common /errtol/ tolrtao,tolrphi,toldlmda,toldtaorm,tol01
      data
     1    r0   ,rp5  ,r1   ,r2   ,r3   /
     2    0.0d0,0.5d0,1.0d0,2.0d0,3.0d0/

c******************initialization********************
      do 10 i=1,9
      do 10 j=1,9
      ginv(i,j)=r0
      gpmtx(i,j)=r0
      gmtx(i,j)=r0
10    continue
      do i=1,3
      un(i)=r0
      do j=1,3
      alf(i,j)=r0
      if(i.eq.j)then
      deta(i,j)=r1
      else
      deta(i,j)=r0
      endif
      enddo
      enddo

c*********** parameters **********************
      prmnt=r0
      do 20 i=1,3
      prmnt=prmnt+pr(i)*pr(i)
20    continue
      if(prmnt.gt.tol01)then
      prmnt=dsqrt(prmnt)
      do 30 i=1,3
      un(i)=pr(i)/prmnt
30    continue
      do 40 i=1,3
      do 40 j=1,3
      alf(i,j)=deta(i,j)-un(i)*un(j)
40    continue
      end if
      almd=yong*poiss/(r1+poiss)/(r1-r2*poiss)
      a2sh=yong/(r1+poiss)

c
c*************** ginv **************************************
c
c        *** elastic compliance****************
      ginv(1,1)=r1/yong
      ginv(2,2)=r1/yong
      ginv(3,3)=r1/yong
      ginv(1,2)=-poiss/yong
      ginv(2,1)=ginv(1,2)
      ginv(1,3)=-poiss/yong
      ginv(3,1)=ginv(1,3)
      ginv(2,3)=-poiss/yong
      ginv(3,2)=ginv(2,3)
      ginv(4,4)=r2/a2sh
      ginv(5,5)=r2/a2sh
      ginv(6,6)=r2/a2sh
c        *** piezo constants****************
      if(prmnt.gt.tol01)then
      do 50 i=1,3
      do 50 k=1,6
      if(k.le.3)then
      ki=k
      kj=k
      elseif(k.eq.4)then
      ki=1
      kj=2
      elseif(k.eq.5)then
      ki=1
      kj=3
      else
      ki=2
      kj=3
      end if
      ginv(i+6,k)
     *   =prmnt/p0*(
     1    dpz33*un(i)*un(ki)*un(kj)+
     2    dpz31*un(i)*alf(ki,kj)+
     3    r1/r2*dpz15*(un(ki)*alf(kj,i)+un(kj)*alf(ki,i)))
50    continue
      do 60 i=1,6
      do 60 j=1,3
      ginv(i,j+6)=ginv(j+6,i)
60    continue
      end if
c        *** dielectric constants****************
      do 70 i=1,3
      do 70 j=1,3
      ginv(i+6,j+6)=dielc*deta(i,j)
70    continue

c
c*************** gpmtx **************************************
c
c        *** elastic stiffness****************
      gpmtx(1,1)=almd+a2sh
      gpmtx(2,2)=almd+a2sh
      gpmtx(3,3)=almd+a2sh
      gpmtx(1,2)=almd
      gpmtx(2,1)=almd
      gpmtx(1,3)=almd
      gpmtx(3,1)=almd
      gpmtx(2,3)=almd
      gpmtx(3,2)=almd
      gpmtx(4,4)=a2sh/r2
      gpmtx(5,5)=a2sh/r2
      gpmtx(6,6)=a2sh/r2
c        *** piezo constants****************
      do i=1,3
       do j=1,6
        do ii=1,6
         gpmtx(i+6,j)=gpmtx(i+6,j)+ginv(i+6,ii)*gpmtx(ii,j)
        end do
       end do
      end do
      do i=1,6
       do j=1,3
        gpmtx(i,j+6)=-gpmtx(j+6,i)
       end do
      end do
c        *** dielectric constants****************
      do i=1,3
       do j=1,3
        gpmtx(i+6,j+6)=gpmtx(i+6,j+6)+ginv(i+6,j+6)
        do ii=1,6
         gpmtx(i+6,j+6)=gpmtx(i+6,j+6)+ginv(i+6,ii)*gpmtx(ii,j+6)
        end do
       end do
      end do

c
c*************** gmtx **************************************
c
c        *** dielectric constants****************
c
c         ***inverse of the dielectric matrix in gpmtx**********
      do i=1,3
       do j=1,3
        temp1(i,j)=gpmtx(i+6,j+6)
       end do
      end do
      call mtxinv(3,temp1,temp2)
c         ********************************************************
      do i=1,3
       do j=1,3
        gmtx(i+6,j+6)=temp2(i,j)
       end do
      end do
c        *** piezo constants****************
      do i=1,3
       do j=1,6
        do ii=1,3
         gmtx(i+6,j)=gmtx(i+6,j)-gmtx(i+6,ii+6)*gpmtx(ii+6,j)
        end do
       end do
      end do
      do i=1,6
       do j=1,3
        gmtx(i,j+6)=gmtx(j+6,i)
       end do
      end do
c        *** elastic stiffness****************
      do i=1,6
       do j=1,6
        gmtx(i,j)=gmtx(i,j)+gpmtx(i,j)
        do ii=1,3
         gmtx(i,j)=gmtx(i,j)+gpmtx(i,ii+6)*gmtx(ii+6,j)
        end do
       end do
      end do

      return
      end subroutine gmat








      subroutine mtxinv(n,a,ainv)

c last modified: 2010/12/30
c a general purpose matrix inverter by augmenting-pivoting technique:

c  a b c | 1 0 0     1 0 0 | j k l
c  d e f | 0 1 0 =>  0 1 0 | m n o
c  g h i | 0 0 1     0 0 1 | p q r

c based on a lecture by prof. mcfarland
c explanation of passed parameters:
c n: dimension of square matrix
c a: square matrix of dimension nxn to be inverted
c ainv: the inverted matrix

      implicit real*8 (a-h,o-z)
      common /io/ it,in
       dimension a(n,n),ainv(n,n),b(n,2*n)
        data tol /1.0e-30/
c make augmented matrix
      do i=1,n
      do j=1,n
      b(i,j)=0.0d0
      b(i,j+n)=0.0d0

      b(i,j)=a(i,j)
      if(i.eq.j) then
      b(i,j+n)=1.0d0
      end if
      end do
      end do

      do i=1,n

c choose the leftmost non-zero element as pivot
      do j=1,n
      if(dabs(b(i,j)).gt.tol)then
           pivot=b(i,j)
           exit
       else
           if(j.eq.n)then
               write(25,10)i
           end if
      end if
      end do
c step 1: change the chosen pivot into "1" by dividing
c the pivot's row by the pivot number
      do j=1,2*n
      b(i,j)=b(i,j)/pivot
      end do
      pivot=b(i,i) !update pivot value
c step 2: change the remainder of the pivot's column into 0's
c by adding to each row a suitable multiple of the pivot row
      do k=1,n !row
       if(k.ne.i) then
       xnum=b(k,i)/pivot !same column with the current pivot
       do j=1,2*n !col
      b(k,j)=b(k,j)-xnum*b(i,j)
      end do
      end if
      end do

      end do
c prepare the final inverted matrix
      do i=1,n
      do j=1,n
      ainv(i,j)=b(i,j+n)
      end do
      end do
c
10    format(/,8x,'cannot find a non-zero  pivot for row #:',i6,/)
      return
      end subroutine mtxinv !(n,a,ainv)




      subroutine suctsp
     1(almda,tao,taorm,sigma,ctmtx,ictflg)
c*******************************************************************************
c     given:  taorm, tao(1~6), sigma(7~9)
c     update: taorm, tao(7~9), sigma(1~6)
c*******************************************************************************
      implicit double precision (a-h,o-z)
      parameter (itmxmum=50)
      dimension
     1 tao(9), !input
     2 taorm(9), !input, output
     2 sigma(9),ctmtx(9,9), !output
     3 ptaorm(9),pr(3),gmtx(9,9),gpmtx(9,9),ginv(9,9),sigmabr(9),
     4 umtx(9,9),amtx(9,9),sigmabk(9),sigmaht(9),hmtx(9,9),anvct(9),
     5 amvct(9),vmtx(9,9),wmtx(9,9),rtao(9),bmtx(9,9),ajmtx(9,9),
     6 ajinv(9,9),deta(9,9),dtaorm(9),un(3),temp(9,9),qmtx(9,9),
     7 gbrmtx(9,9),anbrvct(9),albrvct(9),
     8 ct11mtx(6,6),ct12mtx(6,3),ct21mtx(3,6),ct22mtx(3,3),ct22inv(3,3)!local
!      common /matconst/ yong,poiss,dielc,dpz33,dpz31,dpz15,p0,e0,he0,ame
      common /errtol/ tolrtao,tolrphi,toldlmda,toldtaorm,tol01
      common /io/ it,in

c     initialization
      do i=1,9
      do j=1,9
          ctmtx(i,j)=0.0d0
      if(i.eq.j)then
              deta(i,j)=1.0d0
      else
              deta(i,j)=0.0d0
      endif
      end do
      end do

c     save taorm of the previous converged step
      do 10 i=1,9
      ptaorm(i)=taorm(i)
10    continue

c     ****** predictor step: compute the tiral state ******
      do i=1,3
      pr(i)=taorm(i+6)
      enddo
c     solve for sigma(1~6),tao(7~9)
      call gmat(pr,
     1          ginv,gpmtx,gmtx)
c
      do i=1,6
        sigma(i)=0.0d0
      do i1=1,6
        sigma(i)=sigma(i)+gpmtx(i,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        sigma(i)=sigma(i)+gpmtx(i,i1+6)*sigma(i1+6)
      end do
      end do
c
      do i=1,3
        tao(i+6)=taorm(i+6)
      do i1=1,6
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1+6)*sigma(i1+6)
      end do
      end do
c     sigmaht, phi
      call dgmat(pr,sigma,
     1           sigmabr,umtx,amtx)
      call flhdelc(pr,sigma,sigmabr,
     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
c     check consistency condition for tiral step
      if(phi.lt.0.0d0) then
c
          do i=1,9
              do j=1,9
                  ctmtx(i,j)=gmtx(i,j)
              end do
          end do
c
          return
      endif

c     ****** corrector step: solve the residual equations ******
c     set guess vale of pr when entering the corrector step for the first time
      prmnt=0.0d0
      do i=1,3
          prmnt=prmnt+pr(i)*pr(i)
      end do
      prmnt=dsqrt(prmnt)
      if(prmnt.lt.tol01)then
          beta=1.0d0/dielc
          epcnst=beta*he0/(beta+he0)
          stsmnt=0.0d0
          do i=1,3
              stsmnt=stsmnt+sigma(i+6)*sigma(i+6)
          end do
          stsmnt=dsqrt(stsmnt)
          dstsmnt=stsmnt-e0
          do i=1,3
              un(i)=sigma(i+6)/stsmnt
          end do
          edisp=dstsmnt/epcnst
          edispe=dstsmnt/beta
          prmnt=edisp-edispe
          do i=1,3
              pr(i)=prmnt*un(i)
          end do
          do i=1,3
              taorm(i+6)=pr(i)
          end do
      end if
      call dgmat(pr,sigma,
     1           sigmabr,umtx,amtx)
      call flhdelc(pr,sigma,sigmabr,
     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
c
c___________________________________________________
c___________________________________________________
      i100=1
c___________________________________________________
c___________________________________________________
      iter=0
500   iter=iter+1
c     initialization
      do i=1,9
      rtao(i)=0.0d0
      dtaorm(i)=0.0d0
      do j=1,9
      ajinv(i,j)=0.0d0
      end do
      end do
c     compute residuals
      do 30 i=1,9
      rtao(i)=rtao(i)-taorm(i)+ptaorm(i)+almda*anvct(i)
30    continue
      errphi=dabs(phi)/e0/e0
c     check for convergence
      anmstnrm=0.0d0
      anmprm=0.0d0
      do 40 i=1,6
      anmstnrm=anmstnrm+ptaorm(i)*ptaorm(i)
40    continue
      do 50 i=1,3
      anmprm=anmprm+ptaorm(i+6)*ptaorm(i+6)
50    continue
      anorm1=0.0d0
      anorm2=0.0d0
      do 60 i=1,6
      anorm1=anorm1+rtao(i)*rtao(i)
60    continue
      if (anmstnrm.gt.tol01)then
      anorm1=anorm1/anmstnrm
      endif
      do 70 i=1,3
      anorm2=anorm2+rtao(i+6)*rtao(i+6)
70    continue
      if(anmprm.gt.tol01)then
      anorm2=anorm2/anmprm
      endif
      errtao=anorm1+anorm2
      errtao=dsqrt(errtao)
c___________________________________________________
c___________________________________________________
      if(iter.eq.i100)then
c          write(25,2010)errtao,errphi,errdtaorm,errdlmda
          i100=i100+5
      end if
c___________________________________________________
c___________________________________________________
c
      if(iter.gt.1)then
      if((errtao.lt.tolrtao .and. errphi.lt.tolrphi) .and.
     1   (errdlmda.lt.toldlmda .and. errdtaorm.lt.toldtaorm))then
c__________________________________________________________________________________________
c       *** compute the consistent tangent stiffness matix ***
c  ********************************************************************
c   ictflg=1: consistent tangent stiffness using double inverse
c   ictflg=2: consistent tangent stiffness using single inverse
c   ictflg>2: continuum tangent stiffness
c  ********************************************************************
c--------------------------------------------------------------------------
          if(ictflg.eq.1)then
c--------------------------------------------------------------------------
          do i=1,9
              do j=1,9
                  vmtx(i,j)=almda*vmtx(i,j)
                  wmtx(i,j)=almda*wmtx(i,j)
              end do
          end do
              do i=1,9
                  do j=1,9
                      temp(i,j)=deta(i,j)-wmtx(i,j)
                      do i1=1,9
                          temp(i,j)=temp(i,j)+vmtx(i,i1)*
     1                              (hmtx(i1,j)-umtx(i1,j))
                      end do
                  end do
              end do
          call mtxinv(9,temp,qmtx) !!!******************
          do 71 i=1,9
          do 71 j=1,9
          temp(i,j)=ginv(i,j)
          do 71 i1=1,9
          do 71 i2=1,9
          do 71 i3=1,9
          temp(i,j)=temp(i,j)+(deta(i,i1)+amtx(i,i1))*qmtx(i1,i2)*
     1               vmtx(i2,i3)*(deta(i3,j)+amtx(j,i3))
71        continue
          call mtxinv(9,temp,gbrmtx) !!!******************
          do 72 i=1,9
          anbrvct(i)=0.0d0
          do 72 i1=1,9
          do 72 i2=1,9
          anbrvct(i)=anbrvct(i)+(deta(i,i1)+amtx(i,i1))*
     1                qmtx(i1,i2)*anvct(i2)
72        continue
          do 73 i=1,9
          albrvct(i)=0.0d0
          do 73 i1=1,9
          albrvct(i)=albrvct(i)+anvct(i1)*(deta(i1,i)+amtx(i,i1))
          do 73 i2=1,9
          do 73 i3=1,9
          albrvct(i)=albrvct(i)+amvct(i1)*qmtx(i1,i2)*vmtx(i2,i3)*
     1                (deta(i3,i)+amtx(i,i3))
          do 73 i4=1,9
          albrvct(i)=albrvct(i)-anvct(i1)*(hmtx(i1,i2)-umtx(i1,i2))*
     1                qmtx(i2,i3)*vmtx(i3,i4)*(deta(i4,i)+amtx(i,i4))
73        continue
          dbr=0.0d0
          do 74 i1=1,9
          do 74 i2=1,9
          dbr=dbr-amvct(i1)*qmtx(i1,i2)*anvct(i2)
          do 74 i3=1,9
          dbr=dbr+anvct(i1)*(hmtx(i1,i2)-umtx(i1,i2))*
     1             qmtx(i2,i3)*anvct(i3)
74        continue
          ctdinom=dbr
          do 75 i1=1,9
          do 75 i2=1,9
          ctdinom=ctdinom+albrvct(i1)*gbrmtx(i1,i2)*albrvct(i2)
75        continue
          do 76 i=1,9
          do 76 j=1,9
          ctmtx(i,j)=gbrmtx(i,j)
          do 76 i1=1,9
          do 76 i2=1,9
          ctmtx(i,j)=ctmtx(i,j)-gbrmtx(i,i1)*anbrvct(i1)*
     1                          albrvct(i2)*gbrmtx(i2,j)/ctdinom
76        continue
c--------------------------------------------------------------------------
          else if(ictflg.eq.2)then
c--------------------------------------------------------------------------
          call bmat(4,amtx,umtx,gmtx,hmtx,
     1              bmtx)
          do i=1,9
              do j=1,9
                  ajinv(i,j)=0.0d0
              end do
          end do
          do 81 i=1,9
          do 81 j=1,9
          ajinv(i,j)=ajinv(i,j)+deta(i,j)-almda*wmtx(i,j)
          do 81 ii=1,9
          ajinv(i,j)=ajinv(i,j)+almda*vmtx(i,ii)*bmtx(ii,j)
81        continue
          call mtxinv(9,ajinv,ajmtx)  !!!******************

          do i=1,9
              do j=1,9
                  vmtx(i,j)=almda*vmtx(i,j)
                  wmtx(i,j)=almda*wmtx(i,j)
              end do
          end do
          dbr=0.0d0
          do 82 i1=1,9
          do 82 i2=1,9
          dbr=dbr-amvct(i1)*ajmtx(i1,i2)*anvct(i2)
          do 82 i3=1,9
          dbr=dbr+anvct(i1)*bmtx(i1,i2)*ajmtx(i2,i3)*anvct(i3)
82        continue
          do 83 i=1,9
          anbrvct(i)=0.0d0
          albrvct(i)=0.0d0
          do 83 i1=1,9
          do 83 i2=1,9
          albrvct(i)=albrvct(i)+anvct(i1)*(deta(i1,i2)+amtx(i2,i1))*
     1                           gmtx(i2,i)
          do 83 i3=1,9
          anbrvct(i)=anbrvct(i)+gmtx(i,i1)*(deta(i1,i2)+amtx(i1,i2))*
     1                           ajmtx(i2,i3)*anvct(i3)
          do 83 i4=1,9
          albrvct(i)=albrvct(i)+amvct(i1)*ajmtx(i1,i2)*vmtx(i2,i3)*
     1                          (deta(i3,i4)+amtx(i4,i3))*gmtx(i4,i)
          do 83 i5=1,9
          albrvct(i)=albrvct(i)-anvct(i1)*bmtx(i1,i2)*ajmtx(i2,i3)*
     1               vmtx(i3,i4)*(deta(i4,i5)+amtx(i5,i4))*gmtx(i5,i)
83        continue
          do 84 i=1,9
          do 84 j=1,9
          ctmtx(i,j)=gmtx(i,j)-anbrvct(i)*albrvct(j)/dbr
          do 84 i1=1,9
          do 84 i2=1,9
          do 84 i3=1,9
          do 84 i4=1,9
          do 84 i5=1,9
          ctmtx(i,j)=ctmtx(i,j)-gmtx(i,i1)*(deta(i1,i2)+amtx(i1,i2))*
     1                          ajmtx(i2,i3)*vmtx(i3,i4)*
     2                          (deta(i4,i5)+amtx(i5,i4))*gmtx(i5,j)
84        continue
c--------------------------------------------------------------------------
          else
c--------------------------------------------------------------------------
          call bmat(4,amtx,umtx,gmtx,hmtx,
     1              bmtx)
          dbr=0.0d0
          do 91 i1=1,9
          dbr=dbr-amvct(i1)*anvct(i1)
          do 91 i2=1,9
          dbr=dbr+anvct(i1)*bmtx(i1,i2)*anvct(i2)
91        continue
          do 92 i=1,9
          anbrvct(i)=0.0d0
          do 92 i1=1,9
          do 92 i2=1,9
          anbrvct(i)=anbrvct(i)+gmtx(i,i1)*(deta(i1,i2)+amtx(i1,i2))*
     1                           anvct(i2)
92        continue
          do 93 i=1,9
          do 93 j=1,9
          ctmtx(i,j)=gmtx(i,j)-anbrvct(i)*anbrvct(j)/dbr
93        continue
c
          end if
c--------------------------------------------------------------------------
c     ctmtx for scalar potential
c--------------------------------------------------------------------------
          do i=1,6
          do j=1,6
            ct11mtx(i,j)=ctmtx(i,j)
          end do
          do j=1,3
            ct12mtx(i,j)=ctmtx(i,j+6)
          end do
          end do
c
          do i=1,3
          do j=1,6
            ct21mtx(i,j)=ctmtx(i+6,j)
          end do
          do j=1,3
            ct22mtx(i,j)=ctmtx(i+6,j+6)
          end do
          end do
c
          call mtxinv(3,ct22mtx,ct22inv)
          ctmtx=0.0d0
c
          do i=1,6
           do j=1,6
            ctmtx(i,j)=ct11mtx(i,j)
            do i1=1,3
             do i2=1,3
              ctmtx(i,j)=ctmtx(i,j)-ct12mtx(i,i1)*ct22inv(i1,i2)*
     1                                            ct21mtx(i2,j)
             end do
            end do
           end do
           do j=1,3
            do i1=1,3
             ctmtx(i,j+6)=ctmtx(i,j+6)+ct12mtx(i,i1)*ct22inv(i1,j)
            end do
           end do
          end do
c
          do i=1,3
           do j=1,6
            do i1=1,3
             ctmtx(i+6,j)=ctmtx(i+6,j)-ct22inv(i,i1)*ct21mtx(i1,j)
            end do
           end do
           do j=1,3
            ctmtx(i+6,j+6)=ct22inv(i,j)
           end do
          end do
c ___________________________________________________________________________________________
c
              write(25,1030)iter-1,errtao,errphi
              write(25,1040)iter-1,errdtaorm,errdlmda
              return
      endif
      end if
      if(iter.gt.itmxmum)then
              write(25,1010)itmxmum
              write(25,1020)errtao,errphi
              write(25,1050)errdtaorm,errdlmda
              stop
              return
      endif
c     compute inverted reduced jacobian
      call bmat(2,amtx,umtx,gmtx,hmtx,
     1          bmtx)
      do 80 i=1,9
      do 80 j=1,9
      ajinv(i,j)=ajinv(i,j)+deta(i,j)-almda*wmtx(i,j)
      do 80 ii=1,9
      ajinv(i,j)=ajinv(i,j)+almda*vmtx(i,ii)*bmtx(ii,j)
80    continue
      call mtxinv(9,ajinv,ajmtx)
c     calculate the increment to consistent parameter
      anumer=0.0d0
      adenom=0.0d0
      anumer=anumer+phi
      do 90 i=1,9
      do 90 j=1,9
      anumer=anumer+amvct(i)*ajmtx(i,j)*rtao(j)
      adenom=adenom-amvct(i)*ajmtx(i,j)*anvct(j)
      do 90 ii=1,9
      anumer=anumer-anvct(i)*bmtx(i,j)*ajmtx(j,ii)*rtao(ii)
      adenom=adenom+anvct(i)*bmtx(i,j)*ajmtx(j,ii)*anvct(ii)
90    continue
      dlmda=anumer/adenom
c     update variables
      almda=almda+dlmda
      do i=1,9
      do j=1,9
          dtaorm(i)=dtaorm(i)+ajmtx(i,j)*(rtao(j)+anvct(j)*dlmda)
      end do
          taorm(i)=taorm(i)+dtaorm(i)
      end do
      do i=1,3
      pr(i)=taorm(i+6)
      end do
c     update sigma(1~6),tao(7~9)
      call gmat(pr,
     1          ginv,gpmtx,gmtx)
c
      do i=1,6
        sigma(i)=0.0d0
      do i1=1,6
        sigma(i)=sigma(i)+gpmtx(i,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        sigma(i)=sigma(i)+gpmtx(i,i1+6)*sigma(i1+6)
      end do
      end do
c
      do i=1,3
        tao(i+6)=taorm(i+6)
      do i1=1,6
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1)*(tao(i1)-taorm(i1))
      end do
      do i1=1,3
        tao(i+6)=tao(i+6)+gpmtx(i+6,i1+6)*sigma(i1+6)
      end do
      end do
c     update sigmaht,phi
      call dgmat(pr,sigma,
     1           sigmabr,umtx,amtx)
      call flhdelc(pr,sigma,sigmabr,
     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
c     relative change of almda and taorm, for convergence-check use
      errdlmda=dabs(dlmda/almda)
      do i=1,9
       errdtaorm1=dtaorm(i)**2
       errdtaorm2=taorm(i)**2
      end do
      errdtaorm=errdtaorm1/errdtaorm2
      errdtaorm=dsqrt(errdtaorm)
c
      goto 500

1010  format(/,8x,'***** convergence criterion is not satisfied
     1 after muximum #:',i9,'   of iteration allowed is reached*****')
1020  format(/,8x,'***** error after itmxmum iteration is:  errtao:',
     1 e12.4,'  errphi:',e12.4)
1030  format(/,8x,'***** converged!!!   error after',i6,
     1'   iterations is,  errtao:',e12.4,'  errphi:',e12.4)
1040  format(/,8x,'***** converged!!!   error after',i6,
     1'   iterations is,  errdtaorm:',e12.4,'  errdlmda:',e12.4)
1050  format(/,8x,'***** error after itmxmum iteration is:  errdtaorm:',
     1 e12.4,'  errdlmda:',e12.4)
2010  format(3x,4e12.4)

      end subroutine suctsp







!
!      subroutine suctsp
!     1(almda,tao,taorm,sigma,ctmtx,ictflg)
!c*******************************************************************************
!c     given:  taorm, tao(1~6), sigma(7~9)
!c     update: taorm, tao(7~9), sigma(1~6)
!c*******************************************************************************
!      implicit double precision (a-h,o-z)
!      parameter (itmxmum=50)
!      dimension
!     1 tao(9), !input
!     2 taorm(9), !input, output
!     2 sigma(9),ctmtx(9,9), !output
!     3 ptaorm(9),pr(3),gmtx(9,9),gpmtx(9,9),ginv(9,9),sigmabr(9),
!     4 umtx(9,9),amtx(9,9),sigmabk(9),sigmaht(9),hmtx(9,9),anvct(9),
!     5 amvct(9),vmtx(9,9),wmtx(9,9),rtao(9),bmtx(9,9),ajmtx(9,9),
!     6 ajinv(9,9),deta(9,9),dtaorm(9),un(3),temp(9,9),qmtx(9,9),
!     7 gbrmtx(9,9),anbrvct(9),albrvct(9),
!     8 ct11mtx(6,6),ct12mtx(6,3),ct21mtx(3,6),ct22mtx(3,3),ct22inv(3,3)!local
!      common /matconst/ yong,poiss,dielc,dpz33,dpz31,dpz15,p0,e0,he0,ame
!      common /errtol/ tolrtao,tolrphi,toldlmda,toldtaorm,tol01
!      common /io/ it,in
!
!c     initialization
!      do i=1,9
!      do j=1,9
!          ctmtx(i,j)=0.0d0
!      if(i.eq.j)then
!              deta(i,j)=1.0d0
!      else
!              deta(i,j)=0.0d0
!      endif
!      end do
!      end do
!
!c     save taorm of the previous converged step
!      do 10 i=1,9
!      ptaorm(i)=taorm(i)
!10    continue
!
!c     ****** predictor step: compute the tiral state ******
!      do i=1,3
!      pr(i)=taorm(i+6)
!      enddo
!c     solve for sigma(1~6),tao(7~9)
!      call gmat(pr,
!     1          ginv,gpmtx,gmtx)
!c
!      do i=1,6
!        sigma(i)=0.0d0
!      do i1=1,6
!        sigma(i)=sigma(i)+gpmtx(i,i1)*(tao(i1)-taorm(i1))
!      end do
!      do i1=1,3
!        sigma(i)=sigma(i)+gpmtx(i,i1+6)*sigma(i1+6)
!      end do
!      end do
!c
!      do i=1,3
!        tao(i+6)=taorm(i+6)
!      do i1=1,6
!        tao(i+6)=tao(i+6)+gpmtx(i+6,i1)*(tao(i1)-taorm(i1))
!      end do
!      do i1=1,3
!        tao(i+6)=tao(i+6)+gpmtx(i+6,i1+6)*sigma(i1+6)
!      end do
!      end do
!c     sigmaht, phi
!      call dgmat(pr,sigma,
!     1           sigmabr,umtx,amtx)
!      call flhdelc(pr,sigma,sigmabr,
!     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
!c     check consistency condition for tiral step
!      if(phi.lt.0.0d0) then
!c
!          do i=1,9
!              do j=1,9
!                  ctmtx(i,j)=gmtx(i,j)
!              end do
!          end do
!c
!          return
!      endif
!
!c     ****** corrector step: solve the residual equations ******
!c     set guess vale of pr when entering the corrector step for the first time
!      prmnt=0.0d0
!      do i=1,3
!          prmnt=prmnt+pr(i)*pr(i)
!      end do
!      prmnt=dsqrt(prmnt)
!      if(prmnt.lt.tol01)then
!          beta=1.0d0/dielc
!          epcnst=beta*he0/(beta+he0)
!          stsmnt=0.0d0
!          do i=1,3
!              stsmnt=stsmnt+sigma(i+6)*sigma(i+6)
!          end do
!          stsmnt=dsqrt(stsmnt)
!          dstsmnt=stsmnt-e0
!          do i=1,3
!              un(i)=sigma(i+6)/stsmnt
!          end do
!          edisp=dstsmnt/epcnst
!          edispe=dstsmnt/beta
!          prmnt=edisp-edispe
!          do i=1,3
!              pr(i)=prmnt*un(i)
!          end do
!          do i=1,3
!              taorm(i+6)=pr(i)
!          end do
!      end if
!      call dgmat(pr,sigma,
!     1           sigmabr,umtx,amtx)
!      call flhdelc(pr,sigma,sigmabr,
!     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
!c
!c___________________________________________________
!c___________________________________________________
!      i100=1
!c___________________________________________________
!c___________________________________________________
!      iter=0
!500   iter=iter+1
!c     initialization
!      do i=1,9
!      rtao(i)=0.0d0
!      dtaorm(i)=0.0d0
!      do j=1,9
!      ajinv(i,j)=0.0d0
!      end do
!      end do
!c     compute residuals
!      do 30 i=1,9
!      rtao(i)=rtao(i)-taorm(i)+ptaorm(i)+almda*anvct(i)
!30    continue
!      errphi=dabs(phi)/e0/e0
!c     check for convergence
!      anmstnrm=0.0d0
!      anmprm=0.0d0
!      do 40 i=1,6
!      anmstnrm=anmstnrm+ptaorm(i)*ptaorm(i)
!40    continue
!      do 50 i=1,3
!      anmprm=anmprm+ptaorm(i+6)*ptaorm(i+6)
!50    continue
!      anorm1=0.0d0
!      anorm2=0.0d0
!      do 60 i=1,6
!      anorm1=anorm1+rtao(i)*rtao(i)
!60    continue
!      if (anmstnrm.gt.tol01)then
!      anorm1=anorm1/anmstnrm
!      endif
!      do 70 i=1,3
!      anorm2=anorm2+rtao(i+6)*rtao(i+6)
!70    continue
!      if(anmprm.gt.tol01)then
!      anorm2=anorm2/anmprm
!      endif
!      errtao=anorm1+anorm2
!      errtao=dsqrt(errtao)
!c___________________________________________________
!c___________________________________________________
!      if(iter.eq.i100)then
!c          write(25,2010)errtao,errphi,errdtaorm,errdlmda
!          i100=i100+5
!      end if
!c___________________________________________________
!c___________________________________________________
!c
!      if(iter.gt.1)then
!      if((errtao.lt.tolrtao .and. errphi.lt.tolrphi) .and.
!     1   (errdlmda.lt.toldlmda .and. errdtaorm.lt.toldtaorm))then
!c__________________________________________________________________________________________
!c       *** compute the consistent tangent stiffness matix ***
!c  ********************************************************************
!c   ictflg=1: consistent tangent stiffness using double inverse
!c   ictflg=2: consistent tangent stiffness using single inverse
!c   ictflg>2: continuum tangent stiffness
!c  ********************************************************************
!c--------------------------------------------------------------------------
!          if(ictflg.eq.1)then
!c--------------------------------------------------------------------------
!          do i=1,9
!              do j=1,9
!                  vmtx(i,j)=almda*vmtx(i,j)
!                  wmtx(i,j)=almda*wmtx(i,j)
!              end do
!          end do
!              do i=1,9
!                  do j=1,9
!                      temp(i,j)=deta(i,j)-wmtx(i,j)
!                      do i1=1,9
!                          temp(i,j)=temp(i,j)+vmtx(i,i1)*
!     1                              (hmtx(i1,j)-umtx(i1,j))
!                      end do
!                  end do
!              end do
!          call mtxinv(9,temp,qmtx) !!!******************
!          do 71 i=1,9
!          do 71 j=1,9
!          temp(i,j)=ginv(i,j)
!          do 71 i1=1,9
!          do 71 i2=1,9
!          do 71 i3=1,9
!          temp(i,j)=temp(i,j)+(deta(i,i1)+amtx(i,i1))*qmtx(i1,i2)*
!     1               vmtx(i2,i3)*(deta(i3,j)+amtx(j,i3))
!71        continue
!          call mtxinv(9,temp,gbrmtx) !!!******************
!          do 72 i=1,9
!          anbrvct(i)=0.0d0
!          do 72 i1=1,9
!          do 72 i2=1,9
!          anbrvct(i)=anbrvct(i)+(deta(i,i1)+amtx(i,i1))*
!     1                qmtx(i1,i2)*anvct(i2)
!72        continue
!          do 73 i=1,9
!          albrvct(i)=0.0d0
!          do 73 i1=1,9
!          albrvct(i)=albrvct(i)+anvct(i1)*(deta(i1,i)+amtx(i,i1))
!          do 73 i2=1,9
!          do 73 i3=1,9
!          albrvct(i)=albrvct(i)+amvct(i1)*qmtx(i1,i2)*vmtx(i2,i3)*
!     1                (deta(i3,i)+amtx(i,i3))
!          do 73 i4=1,9
!          albrvct(i)=albrvct(i)-anvct(i1)*(hmtx(i1,i2)-umtx(i1,i2))*
!     1                qmtx(i2,i3)*vmtx(i3,i4)*(deta(i4,i)+amtx(i,i4))
!73        continue
!          dbr=0.0d0
!          do 74 i1=1,9
!          do 74 i2=1,9
!          dbr=dbr-amvct(i1)*qmtx(i1,i2)*anvct(i2)
!          do 74 i3=1,9
!          dbr=dbr+anvct(i1)*(hmtx(i1,i2)-umtx(i1,i2))*
!     1             qmtx(i2,i3)*anvct(i3)
!74        continue
!          ctdinom=dbr
!          do 75 i1=1,9
!          do 75 i2=1,9
!          ctdinom=ctdinom+albrvct(i1)*gbrmtx(i1,i2)*albrvct(i2)
!75        continue
!          do 76 i=1,9
!          do 76 j=1,9
!          ctmtx(i,j)=gbrmtx(i,j)
!          do 76 i1=1,9
!          do 76 i2=1,9
!          ctmtx(i,j)=ctmtx(i,j)-gbrmtx(i,i1)*anbrvct(i1)*
!     1                          albrvct(i2)*gbrmtx(i2,j)/ctdinom
!76        continue
!c--------------------------------------------------------------------------
!          else if(ictflg.eq.2)then
!c--------------------------------------------------------------------------
!          call bmat(4,amtx,umtx,gmtx,hmtx,
!     1              bmtx)
!          do i=1,9
!              do j=1,9
!                  ajinv(i,j)=0.0d0
!              end do
!          end do
!          do 81 i=1,9
!          do 81 j=1,9
!          ajinv(i,j)=ajinv(i,j)+deta(i,j)-almda*wmtx(i,j)
!          do 81 ii=1,9
!          ajinv(i,j)=ajinv(i,j)+almda*vmtx(i,ii)*bmtx(ii,j)
!81        continue
!          call mtxinv(9,ajinv,ajmtx)  !!!******************
!
!          do i=1,9
!              do j=1,9
!                  vmtx(i,j)=almda*vmtx(i,j)
!                  wmtx(i,j)=almda*wmtx(i,j)
!              end do
!          end do
!          dbr=0.0d0
!          do 82 i1=1,9
!          do 82 i2=1,9
!          dbr=dbr-amvct(i1)*ajmtx(i1,i2)*anvct(i2)
!          do 82 i3=1,9
!          dbr=dbr+anvct(i1)*bmtx(i1,i2)*ajmtx(i2,i3)*anvct(i3)
!82        continue
!          do 83 i=1,9
!          anbrvct(i)=0.0d0
!          albrvct(i)=0.0d0
!          do 83 i1=1,9
!          do 83 i2=1,9
!          albrvct(i)=albrvct(i)+anvct(i1)*(deta(i1,i2)+amtx(i2,i1))*
!     1                           gmtx(i2,i)
!          do 83 i3=1,9
!          anbrvct(i)=anbrvct(i)+gmtx(i,i1)*(deta(i1,i2)+amtx(i1,i2))*
!     1                           ajmtx(i2,i3)*anvct(i3)
!          do 83 i4=1,9
!          albrvct(i)=albrvct(i)+amvct(i1)*ajmtx(i1,i2)*vmtx(i2,i3)*
!     1                          (deta(i3,i4)+amtx(i4,i3))*gmtx(i4,i)
!          do 83 i5=1,9
!          albrvct(i)=albrvct(i)-anvct(i1)*bmtx(i1,i2)*ajmtx(i2,i3)*
!     1               vmtx(i3,i4)*(deta(i4,i5)+amtx(i5,i4))*gmtx(i5,i)
!83        continue
!          do 84 i=1,9
!          do 84 j=1,9
!          ctmtx(i,j)=gmtx(i,j)-anbrvct(i)*albrvct(j)/dbr
!          do 84 i1=1,9
!          do 84 i2=1,9
!          do 84 i3=1,9
!          do 84 i4=1,9
!          do 84 i5=1,9
!          ctmtx(i,j)=ctmtx(i,j)-gmtx(i,i1)*(deta(i1,i2)+amtx(i1,i2))*
!     1                          ajmtx(i2,i3)*vmtx(i3,i4)*
!     2                          (deta(i4,i5)+amtx(i5,i4))*gmtx(i5,j)
!84        continue
!c--------------------------------------------------------------------------
!          else
!c--------------------------------------------------------------------------
!          call bmat(4,amtx,umtx,gmtx,hmtx,
!     1              bmtx)
!          dbr=0.0d0
!          do 91 i1=1,9
!          dbr=dbr-amvct(i1)*anvct(i1)
!          do 91 i2=1,9
!          dbr=dbr+anvct(i1)*bmtx(i1,i2)*anvct(i2)
!91        continue
!          do 92 i=1,9
!          anbrvct(i)=0.0d0
!          do 92 i1=1,9
!          do 92 i2=1,9
!          anbrvct(i)=anbrvct(i)+gmtx(i,i1)*(deta(i1,i2)+amtx(i1,i2))*
!     1                           anvct(i2)
!92        continue
!          do 93 i=1,9
!          do 93 j=1,9
!          ctmtx(i,j)=gmtx(i,j)-anbrvct(i)*anbrvct(j)/dbr
!93        continue
!c
!          end if
!c--------------------------------------------------------------------------
!c     ctmtx for scalar potential
!c--------------------------------------------------------------------------
!          do i=1,6
!          do j=1,6
!            ct11mtx(i,j)=ctmtx(i,j)
!          end do
!          do j=1,3
!            ct12mtx(i,j)=ctmtx(i,j+6)
!          end do
!          end do
!c
!          do i=1,3
!          do j=1,6
!            ct21mtx(i,j)=ctmtx(i+6,j)
!          end do
!          do j=1,3
!            ct22mtx(i,j)=ctmtx(i+6,j+6)
!          end do
!          end do
!c
!          call mtxinv(3,ct22mtx,ct22inv)
!          ctmtx=0.0d0
!c
!          do i=1,6
!           do j=1,6
!            ctmtx(i,j)=ct11mtx(i,j)
!            do i1=1,3
!             do i2=1,3
!              ctmtx(i,j)=ctmtx(i,j)-ct12mtx(i,i1)*ct22inv(i1,i2)*
!     1                                            ct21mtx(i2,j)
!             end do
!            end do
!           end do
!           do j=1,3
!            do i1=1,3
!             ctmtx(i,j+6)=ctmtx(i,j+6)+ct12mtx(i,i1)*ct22inv(i1,j)
!            end do
!           end do
!          end do
!c
!          do i=1,3
!           do j=1,6
!            do i1=1,3
!             ctmtx(i+6,j)=ctmtx(i+6,j)-ct22inv(i,i1)*ct21mtx(i1,j)
!            end do
!           end do
!           do j=1,3
!            ctmtx(i+6,j+6)=ct22inv(i,j)
!           end do
!          end do
!c ___________________________________________________________________________________________
!c
!              write(25,1030)iter-1,errtao,errphi
!              write(25,1040)iter-1,errdtaorm,errdlmda
!              return
!      endif
!      end if
!      if(iter.gt.itmxmum)then
!              write(25,1010)itmxmum
!              write(25,1020)errtao,errphi
!              write(25,1050)errdtaorm,errdlmda
!              stop
!              return
!      endif
!c     compute inverted reduced jacobian
!      call bmat(2,amtx,umtx,gmtx,hmtx,
!     1          bmtx)
!      do 80 i=1,9
!      do 80 j=1,9
!      ajinv(i,j)=ajinv(i,j)+deta(i,j)-almda*wmtx(i,j)
!      do 80 ii=1,9
!      ajinv(i,j)=ajinv(i,j)+almda*vmtx(i,ii)*bmtx(ii,j)
!80    continue
!      call mtxinv(9,ajinv,ajmtx)
!c     calculate the increment to consistent parameter
!      anumer=0.0d0
!      adenom=0.0d0
!      anumer=anumer+phi
!      do 90 i=1,9
!      do 90 j=1,9
!      anumer=anumer+amvct(i)*ajmtx(i,j)*rtao(j)
!      adenom=adenom-amvct(i)*ajmtx(i,j)*anvct(j)
!      do 90 ii=1,9
!      anumer=anumer-anvct(i)*bmtx(i,j)*ajmtx(j,ii)*rtao(ii)
!      adenom=adenom+anvct(i)*bmtx(i,j)*ajmtx(j,ii)*anvct(ii)
!90    continue
!      dlmda=anumer/adenom
!c     update variables
!      almda=almda+dlmda
!      do i=1,9
!      do j=1,9
!          dtaorm(i)=dtaorm(i)+ajmtx(i,j)*(rtao(j)+anvct(j)*dlmda)
!      end do
!          taorm(i)=taorm(i)+dtaorm(i)
!      end do
!      do i=1,3
!      pr(i)=taorm(i+6)
!      end do
!c     update sigma(1~6),tao(7~9)
!      call gmat(pr,
!     1          ginv,gpmtx,gmtx)
!c
!      do i=1,6
!        sigma(i)=0.0d0
!      do i1=1,6
!        sigma(i)=sigma(i)+gpmtx(i,i1)*(tao(i1)-taorm(i1))
!      end do
!      do i1=1,3
!        sigma(i)=sigma(i)+gpmtx(i,i1+6)*sigma(i1+6)
!      end do
!      end do
!c
!      do i=1,3
!        tao(i+6)=taorm(i+6)
!      do i1=1,6
!        tao(i+6)=tao(i+6)+gpmtx(i+6,i1)*(tao(i1)-taorm(i1))
!      end do
!      do i1=1,3
!        tao(i+6)=tao(i+6)+gpmtx(i+6,i1+6)*sigma(i1+6)
!      end do
!      end do
!c     update sigmaht,phi
!      call dgmat(pr,sigma,
!     1           sigmabr,umtx,amtx)
!      call flhdelc(pr,sigma,sigmabr,
!     1             sigmabk,sigmaht,hmtx,anvct,amvct,vmtx,wmtx,phi)
!c     relative change of almda and taorm, for convergence-check use
!      errdlmda=dabs(dlmda/almda)
!      do i=1,9
!       errdtaorm1=dtaorm(i)**2
!       errdtaorm2=taorm(i)**2
!      end do
!      errdtaorm=errdtaorm1/errdtaorm2
!      errdtaorm=dsqrt(errdtaorm)
!c
!      goto 500
!
!1010  format(/,8x,'***** convergence criterion is not satisfied
!     1 after muximum #:',i9,'   of iteration allowed is reached*****')
!1020  format(/,8x,'***** error after itmxmum iteration is:  errtao:',
!     1 e12.4,'  errphi:',e12.4)
!1030  format(/,8x,'***** converged!!!   error after',i6,
!     1'   iterations is,  errtao:',e12.4,'  errphi:',e12.4)
!1040  format(/,8x,'***** converged!!!   error after',i6,
!     1'   iterations is,  errdtaorm:',e12.4,'  errdlmda:',e12.4)
!1050  format(/,8x,'***** error after itmxmum iteration is:  errdtaorm:',
!     1 e12.4,'  errdlmda:',e12.4)
!2010  format(3x,4e12.4)
!
!      return
!      end






















      end module polarization_switching_





