program md  !simple MD program
  implicit none
  double precision:: deltat, t, tmax, vx(500), vy(500), xcord(500), ycord(500),deltax(500), deltay(500),f(500)
  double precision::  deltar(500),kenergy(500), totkenergy
  double precision,parameter:: mass=39.948, sigma=3.405, epsillon=0.01035, npart=500, lx=60 ,ly=60   
  real,external :: ran3
  
  write(*,*) '################################ Welcome to MD Simulation ###############################'
  write(*,*) '############################### 500 Argon particle in cell #########################'

  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) '###############################    simulation start    #########################################'
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''

  
  call init(xcord, ycord, vx, vy, totkenergy)  !initialization
  
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''

  
  t=0
  !tmax=0.00000000001
   tmax=0.00000000001
  write(*,*) tmax
  deltat=tmax/9000
  write(*,*) deltat

  !t=0
  !deltat=0
  !tmax=6
  !deltat=2
  
  
  do while (t<tmax) !MD loop
     call force(xcord, ycord, vx, vy,deltax, deltay, deltar,totkenergy,f,t)  !determine the forces
!write(*,*)'the force complete'
     !   call integrate(f,en) !integrate equations of motion
     t=t+deltat
     !write(*,*) t
     !  call sample  !sample averages
  enddo
  !  stop
  
  
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) ''
  write(*,*) '#####################################  simulation finish  ######################################'
end program md


!#############################################################################################################
!#############################################################################################################


subroutine init(xcord, ycord, vx, vy,totkenergy)
  implicit none
  double precision:: deltat, t, tmax, vx(500), vy(500), xcord(500), ycord(500),deltax(500), deltay(500), f(500)
  double precision::  deltar(500),kenergy(500), totkenergy
  double precision,parameter:: mass=39.948, sigma=3.405, epsillon=0.01035, npart=500,lx=60,ly=60   
  real,external :: ran3
  INTEGER:: idum,i
  
 

  write(*,*) '#################################### The Initialize Start ######################################'

  xcord=0
  ycord=0
  vx=0
  vy=0
  kenergy=0
  totkenergy=0
  
  OPEN(5,file='total.txt')
  
  idum= 40
  do  i=1,npart
     
     xcord(i)=ran3(idum)*60 !this scale to 0 to 60
     
     ycord(i)=ran3(idum)*60  !this scale to 0 to 60
     
     vx(i)=ran3(idum)
     
     vy(i)=ran3(idum)
     
     kenergy(i)=0.5*mass*(vx(i)**2+vy(i)**2)
     
     totkenergy=totkenergy+kenergy(i)

     10 format(I5,5F10.5)
     write(5,10) i, xcord(i), ycord(i), vx(i), vy(i), kenergy(i)
  enddo
  
  
  write(*,*) ''
  write(*,*) '' 
  write(*,*) 'total kinetic energy:',totkenergy
  write(*,*) ''
  
  write(*,*) ''
  write(*,*) ''
  
  write(*,*) '############################################# The Initialize Finish #############################'
  
end subroutine init

!#############################################################################################################
!#############################################################################################################



subroutine force(xcord, ycord, vx, vy,deltax, deltay, deltar,totkenergy,f,t)
  implicit none
  double precision:: deltat, t, tmax, vx(500), vy(500), xcord(500), ycord(500),deltax(500), deltay(500), f(500)
  double precision::  deltar(500),kenergy(500), totkenergy
  double precision,parameter:: mass=39.948, sigma=3.405, epsillon=0.01035, npart=500, lx=60 ,ly=60 
  INTEGER:: i, j , k
  double precision:: rm2,rm6,rm12,phi,dphi
  
  do k=1,npart
     f(k)=0            !set force to zero
  enddo
  
  
  do i=1,npart-1
     do j=i+1, npart
        !write(*,*)'############### the time',t,'npart',i,j,'equal   ##################'
        !write (*,*) 'ghabl'
        deltax(j)=xcord(i)-xcord(j)
        deltay(j)=ycord(i)-ycord(j)        

        !write (*,*) 'deltax',deltax(j),''
        
        !write (*,*)'deltay',deltay(j),''

        !periodic condition in x 
        
        if (deltax(j) < -30) then
           deltax(j)= deltax(j) +60
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) 'manfi'
           !write(*,*) ''
           !write(*,*) ''
           !write (*,*) deltax(j)
        elseif(deltax(j) > 30) then
           deltax(j)= deltax(j)-60
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) 'mosbat'
           !write(*,*) ''
           !write(*,*) ''
           !write (*,*) deltax(j)
        else
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) 'bedone taghir'
           !write(*,*) ''
           !write(*,*) ''
           !write (*,*) deltax(j)   
        endif

        !periodic condition in y
        
        !write (*,*)'deltay'
        
        if (deltay(j) < -30) then
           deltay(j)= deltay(j) +60
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) 'manfi'
           !write(*,*) ''
           !write(*,*) ''
           !write (*,*) deltay(j)
        elseif(deltay(j) > 30) then
           deltay(j)= deltay(j)-60
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) 'mosbat'
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) deltay(j)
        else
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) 'bedone taghir'
           !write(*,*) ''
           !write(*,*) ''
           !write(*,*) deltay(j)   
        endif
       
        deltar(j)= sqrt(deltax(j)**2+deltay(j)**2)
        
        !write(*,*) ''
        !write(*,*) ''
        !write(*,*) '2.5 sigma=',2.5*sigma,''
        !write(*,*) ''
        !write(*,*) '' 
        


        !write(*,*) ''
        !write(*,*) ''
        !write(*,*) 'deltar',deltar(j),''
        !write(*,*) ''
        !write(*,*) ''
        
        
        if(deltar(j) <2.5*sigma) then   ! compute Lennard-Jones potntl
           !write(*,*) "mohasebe shavad"
           rm2=1./deltar(j)  !1/r^2
           rm6=rm2**3        !1/r^6
           rm12=rm2**6       !1/r^12
           phi=4*epsillon*(rm12-rm6)             !4[1/r^12 - 1/r^6] - phi(Rc)
                                                 !  The following is dphi = -(1/r)(dV/dr)
           dphi=24.0*rm2*( 2.d0*rm12 - rm6 )     !  24[2/r^14 - 1/r^8]

           
        !   write(*,*) 'the force is',dphi,'' 
           
           
       else
           ! write(*,*) "mohasebe nashavad"
        end if
     enddo
  enddo
  ! if
  
end subroutine force


!#############################################################################################################
!#############################################################################################################     
subroutine integrate
end subroutine integrate
!#############################################################################################################
!#############################################################################################################
FUNCTION ran3(idum)
  INTEGER idum
  INTEGER MBIG,MSEED,MZ
  ! C     REAL MBIG,MSEED,MZ
  REAL ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
  ! C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
  INTEGER i,iff,ii,inext,inextp,k
  INTEGER mj,mk,ma(55)
  ! C     REAL mj,mk,ma(55)
  SAVE iff,inext,inextp,ma
  DATA iff /0/
  if(idum.lt.0.or.iff.eq.0)then
     iff=1
     mj=abs(MSEED-abs(idum))
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do 11 i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(ii)
11      continue
        do 13 k=1,4
           do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12            continue
13            continue
              inext=0
              inextp=31
              idum=1
           endif
           inext=inext+1
           if(inext.eq.56)inext=1
           inextp=inextp+1
           if(inextp.eq.56)inextp=1
           mj=ma(inext)-ma(inextp)
           if(mj.lt.MZ)mj=mj+MBIG
           ma(inext)=mj
           ran3=mj*FAC
           return
        
end FUNCTION ran3       
!#############################################################################################################
!#############################################################################################################
