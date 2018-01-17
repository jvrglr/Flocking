	subroutine dran_ini(iseed0)

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(np=14)

	parameter(nbit=31)

	parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)

	integer ix(ip)

	dimension g(0:m)



	data c0,c1,c2/2.515517,0.802853,0.010328/

	data d1,d2,d3/1.432788,0.189269,0.001308/

c

	common /ixx/ ix

	common /icc/ ic     

	common /gg/ g

c

	dseed=iseed0

    	do 200 i=1,ip

        ix(i)=0

        do 200 j=0,nbit-1

        if(rand_xx(dseed).lt.0.5) ix(i)=ibset(ix(i),j)

200	continue

	ic=0

c

	pi=4.0d0*datan(1.0d0)

	do 1 i=m/2,m

	p=1.0-real(i+1)/(m+2)

	t=sqrt(-2.0*log(p))

	x=t-(c0+t*(c1+c2*t))/(1.0+t*(d1+t*(d2+t*d3)))

	g(i)=x

	g(m-i)=-x

1	continue



	u2th=1.0-real(m+2)/m*sqrt(2.0/pi)*g(m)*exp(-g(m)*g(m)/2)

	u2th=nn1*sqrt(u2th)

	do 856 i=0,m

856	g(i)=g(i)/u2th



	return

	end



        subroutine dran_read(iunit)

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(np=14)

	parameter(m=2**np)

        integer ix(ip)

	dimension g(0:m)

        common /ixx/ ix

        common /icc/ ic

	common /gg/ g

        read (iunit,*) ic

        read (iunit,*) (ix(i),i=1,ip)

	read (iunit,*) (g(i),i=0,m)

        return

        end



        subroutine dran_write(iunit)

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(np=14)

	parameter(m=2**np)

        integer ix(ip)

	dimension g(0:m)

	common /ixx/ ix

	common /icc/ ic      

	common /gg/ g

        write (iunit,*) ic

        write (iunit,*) (ix(i),i=1,ip)

	write (iunit,*) (g(i),i=0,m)

        return

        end



    	function i_dran(n)

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(iq=418)

	parameter(is=ip-iq)

	integer ix(ip)

	common /ixx/ ix

	common /icc/ ic

	ic=ic+1

    	if(ic.gt.ip) ic=1

	if(ic.gt.iq) then

		ix(ic)=ieor(ix(ic),ix(ic-iq))

	else

        	ix(ic)=ieor(ix(ic),ix(ic+is))

    	endif

    	i_ran=ix(ic)

        if (n.gt.0) i_dran=mod(i_ran,n)+1

    	return

	end



	function dran_u()

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(iq=418)

	parameter(is=ip-iq)

    	parameter (rmax=2147483647.0)

	integer ix(ip)

	common /ixx/ ix

	common /icc/ ic

	ic=ic+1

    	if(ic.gt.ip) ic=1

	if(ic.gt.iq) then

		ix(ic)=ieor(ix(ic),ix(ic-iq))

	else

        	ix(ic)=ieor(ix(ic),ix(ic+is))

    	endif

    	dran_u=real(ix(ic))/rmax

	return

	end





	function dran_g()

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(iq=418)

	parameter(np=14)

	parameter(nbit=31)

	parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)

	parameter(is=ip-iq)



	integer ix(ip)

	dimension g(0:m)



	common /ixx/ ix

	common /icc/ ic

	common /gg/ g



	ic=ic+1

    	if(ic.gt.ip) ic=1

	if(ic.gt.iq) then

		ix(ic)=ieor(ix(ic),ix(ic-iq))

	else

        	ix(ic)=ieor(ix(ic),ix(ic+is))

    	endif

	i=ishft(ix(ic),-np1)

	i2=iand(ix(ic),nn)

	dran_g=i2*g(i+1)+(nn1-i2)*g(i)

	return

	end





	subroutine dran_bm(x1,x2)

	implicit double precision(a-h,o-z)

    	parameter(ip=1279)

	parameter(iq=418)

	parameter(is=ip-iq)

    	parameter (rmax=2147483647.0)

	integer ix(ip)

	common /ixx/ ix

	common /icc/ ic

	data pi2 /6.2831853072/

	ic=ic+1

    	if (ic.gt.ip) ic=1

	if(ic.gt.iq) then

		ix(ic)=ieor(ix(ic),ix(ic-iq))

	else

        	ix(ic)=ieor(ix(ic),ix(ic+is))

    	endif

	u=pi2*real(ix(ic))/rmax

	ic=ic+1

    	if(ic.gt.ip) ic=1

	if(ic.gt.iq) then

		ix(ic)=ieor(ix(ic),ix(ic-iq))

	else

        	ix(ic)=ieor(ix(ic),ix(ic+is))

    	endif

	v=(real(ix(ic))+0.5)/rmax

	v=sqrt(-2.0*log(v))

	x1=cos(u)*v

	x2=sin(u)*v

	return

	end





	subroutine dran_gv(u,n)

	implicit double precision(a-h,o-z)

	parameter(ip=1279)

	parameter(iq=418)

	parameter(np=14)

	parameter(nbit=31)

	parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)

	parameter(is=ip-iq)

	dimension g(0:m)

	dimension u(n)

	dimension ix(ip)

	common /gg/ g

	common /ixx/ ix

	common /icc/ic



        n1=0

10      if (ic.lt.iq) then

          kmax=min(n-n1,iq-ic)

	  do 99 k=1,kmax

  	  ic=ic+1

	  ix(ic)=ieor(ix(ic),ix(ic+is))

	  i=ishft(ix(ic),-np1)

	  i2=iand(ix(ic),nn)

	  u(n1+k)=i2*g(i+1)+(nn1-i2)*g(i)

99        continue

        else

          kmax=min(n-n1,ip-ic)

	  do 98 k=1,kmax

	  ic=ic+1

	  ix(ic)=ieor(ix(ic),ix(ic-iq))

	  i=ishft(ix(ic),-np1)

	  i2=iand(ix(ic),nn)

	  u(n1+k)=i2*g(i+1)+(nn1-i2)*g(i)

98        enddo

        endif

        if(ic.ge.ip) ic=0

        n1=n1+kmax

        if (n1.lt.n) goto 10

        

	return

	end





	subroutine dran_gvv(iv,u,n)

	implicit double precision(a-h,o-z)

	parameter(ip=1279)

	parameter(iq=418)

	parameter(np=14)

	parameter(nbit=31)

	parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)

	parameter(is=ip-iq)

        dimension iv(n+ip)

	dimension g(0:m)

	dimension u(n)

	dimension ix(ip)

	common /gg/ g

	common /ixx/ ix

	common /icc/ic



	do 10 i=ic+1,ip

10	iv(i-ic)=ix(i)

	do 20 i=1,ic

20	iv(ip-ic+i)=ix(i)	



cCDEC$ INIT_DEP_FWD

	do 99 k=1,n

	ir=ieor(iv(k+is),iv(k))

	i=ishft(ir,-np1)

	i2=iand(ir,nn)

	u(k)=i2*g(i+1)+(nn1-i2)*g(i)

	iv(k+ip)=ir

99	continue



	do 11 i=ic+1,ip

11	ix(i)=iv(n+i-ic)

	do 21 i=1,ic

21	ix(i)=iv(n+ip-ic+i)



	return

	end



      function rand_xx(dseed)

      double precision a,c,xm,rm,dseed,rand_xx

      parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)

      dseed=mod(dseed*a+c,xm)

      rand_xx=dseed*rm

      return

      end



