module array
    implicit none
    double precision, allocatable::mi(:),th(:),ph(:),rwork(:),evl_s(:)
    integer,dimension(:),allocatable::x(:),y(:),tag(:)
    complex*16,allocatable:: h(:,:),bb(:,:),work(:)
    double precision::nsum_t(1600),nn1(1600),ee1(1600),e_2(1600),dist1(1600),dist_mi_sum(1600)
    doubleprecision::dos_sum_a(1600),dos_sum_b(1600),dos_sum_c(1600)
end module array

module global
    implicit none
    double precision::t,ua,ub,uc,mu,m_sum,e_sum,u,esum_t,m1,ph1,th1,nsum,e,s_f,sf_A,sf_B_C
    double precision::lambda,n_a,n_b,n_c,n_a_up,n_b_up,n_c_up,n_a_dn,n_b_dn,n_c_dn
    integer::i,j,d,nos,ie,flag_isoc,unit_cells,sw,nms!nos=no of lattice points
endmodule global

program liebmain
    use array
    use global
    use mtmod
    implicit none
    integer::seed,site,it,accept,it_count
    double precision::temp,get_mu_s,de,ep,ec,random1,pr
    double precision::m_temp,th_temp,ph_temp,inpt,filling,fill_inpt


    open(8,file='input.dat',status='unknown')
    do i=1,8
    read(8,*)inpt
    if (i.eq.1) seed=int(inpt)
    if (i.eq.2) d=int(inpt)
    if (i.eq.3) t=dble(inpt)
    if (i.eq.4) ua=dble(inpt)
    if (i.eq.5) ub=dble(inpt)
    if (i.eq.6) uc=dble(inpt)
    if (i.eq.7) nms=int(inpt)
    if (i.eq.8) fill_inpt=dble(inpt)
    
    end do
    close(8)   
    call sgrnd(seed)

    unit_cells=(d/2)**2
    nos=3*unit_cells ! no of sites
    filling=fill_inpt*2*nos         
    allocate(x(nos),y(nos),h((2*nos),(2*nos)),mi(nos),th(nos),ph(nos),tag(nos))
    allocate(evl_s(2*nos))


    print*,'unit cells_____________=',unit_cells
    print*,'sites__________________=',nos
    print*,'matrix dim_____________=',2*nos,'X',2*nos
    print*,'hopping parameter______=',t
    print*,'Ua_____________________=',ua
    print*,'ub_____________________=',ub
    print*,'uc_____________________=',uc
    print*,'Monte-Carlo steps______=',nms
    print*,'filling________________=',fill_inpt,filling

    lambda=0.00d0
    Temp=0.650d0


    call lattice_labeling
    call auxf(0,site) 
    call matgen
    call diagonalization(h,1)
    mu=get_mu_s(filling,temp)
    call energy_cal(temp)
    Ep=esum_t
    it_count=0

    do iT=1,30
        if(iT.le.11) Temp=Temp-0.050d0
        if ((iT.gt.11).and.(iT.lt.21)) Temp=Temp-0.010d0
        if ((iT.gt.21)) Temp=Temp-0.0010d0

        accept=0
        do sw=1,nms
            do site=1,nos!no of site1
                m_temp=mi(site)         ! temporary  auxiliary fields
                th_temp=th(site)          
                ph_temp=ph(site)
                    
                call auxf(1,site)
                call matgen
                call diagonalization(h,1)
                mu=get_mu_s(filling,temp)
                call energy_cal(temp)
                Ec=esum_t
                dE=Ec-Ep
                if(de.le.0)then
                    Ep=Ec
                    accept=accept+1        
                else       
                    call rannum(random1)
                    pr = exp(-(dE)/Temp)
                    if(random1.lt.pr) then
                        Ep=Ec
                        accept=accept+1                    
                    else  
                        mi(site)=m_temp
                        th(site)=th_temp
                        ph(site)=ph_temp   
                    endif 
                endif       
            enddo 

            if((sw.gt.(nms/2)).and.(mod(sw,5).eq.0))then
                do i=1,nos
                    write(1000+it,*) Temp,mi(i),th(i),ph(i) 
                    flush(1000+it)
                enddo       
            endif

        enddo

        write(82,*)it,accept
        flush(82) 
    enddo     

end program liebmain

subroutine lattice_labeling
    use global
    use array
    implicit none
    integer::ix,iy,i1
    i1=-2
    do iy=1,d
        do ix=1,d
            if((mod(iy,2).ne.0).and.(mod(ix,2).ne.0))then   

                i1=i1+3
                x(i1)=ix                     !A
                y(i1)=iy
                tag(i1)=1

                x(i1+1)=ix+1                 !B
                y(i1+1)=iy
                tag(i1+1)=2

                x(i1+2)=ix                   !C
                y(i1+2)=iy+1
                tag(i1+2)=3

            endif
        enddo
    enddo
    i1=0
    do i=1,nos     
        write(12,*) i,x(i),y(i),tag(i)
        flush(12)    
    enddo

endsubroutine lattice_labeling


subroutine matgen
    use array
    use global
    implicit none
    integer::l,k,xi,yi,xd,yd,a,b
    h=complex(0.0d0,0.0d0)
    do l=1,nos
        if(tag(l).eq.1) then
            xi=1
            yi=1
            xd=1
            yd=1
            if(x(l).eq.d) xi=-d+1
            if(y(l).eq.d) yi=-d+1

            if(x(l).eq.1) xd=-d+1
            if(y(l).eq.1) yd=-d+1

            do k=1,nos
                if((tag(k).eq.2).or.((tag(k).eq.3)))then        
                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.y(l)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)+yi)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
        
                    if((x(k).eq.x(l)).and.(y(k).eq.(y(l)-yd)))then  
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=t;h(b,a)=t
                        h(a+1,b+1)=t;h(b+1,a+1)=t
                    endif
    
                endif
            enddo
        endif
    enddo
    call matgen_mcmf_u
    do l=1,2*nos
        do k=1,2*nos
            ! if(h(l,k).ne.0) write(34,*) l,k,h(l,k)
            if(h(l,k).ne.conjg(h(k,l))) print*,l,k,'not hermitian'
        enddo
    enddo 

endsubroutine matgen

subroutine matgen_mcmf_u
    use array
    use global
    use mtmod
    implicit none
    integer::l,a,b

    m_sum=0.0d0

    do l=1,nos
        if(tag(l).eq.1) u=ua
        if(tag(l).eq.2) u=ub
        if(tag(l).eq.3) u=uc

        a=2*l-1
        b=2*l-1
        h(a,b)=(u/2.0d0)*(-mi(l)*cos(th(l)))
        h(a,b+1)=(-u/2.0d0)*mi(l)*sin(th(l))*complex(cos(ph(l)),-sin(ph(l)))
        h(a+1,b)=conjg(h(a,b+1))
        h(a+1,b+1)=(u/2.0d0)*(-mi(l)*(-cos(th(l))))

        m_sum=m_sum+((mi(l))**2.0d0)*(U/4.0d0)
    enddo

endsubroutine matgen_mcmf_u



subroutine matgen_so
    use global
    use array
    implicit none
    integer::l,a,b,k,xi,yi,xd,yd

    do l=1,nos
        if(tag(l).eq.2)then
            xi=1
            yi=1
            xd=1
            yd=1
            if(x(l).eq.d) xi=-d+1
            if(y(l).eq.d) yi=-d+1

            if(x(l).eq.1) xd=-d+1
            if(y(l).eq.1) yd=-d+1
            do k=1,nos
                if((tag(k).eq.3))then

                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.(y(l)+yi)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=-lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))
                    endif

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.(y(l)+yi)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=-lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))

                    endif 

                    if((x(k).eq.(x(l)-xd)).and.(y(k).eq.(y(l)-yd)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=-lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))
                    endif 

                    if((x(k).eq.(x(l)+xi)).and.(y(k).eq.(y(l)-yd)))then 
                        a=2*l-1
                        b=2*k-1
                        h(a,b)=lambda*cmplx(0,1);h(b,a)=conjg(h(a,b))
                        h(a+1,b+1)=-lambda*cmplx(0,1);h(b+1,a+1)=conjg(h(a+1,b+1))

                    endif   

                endif
            enddo
        endif
    enddo
endsubroutine matgen_so

subroutine diagonalization(h_temp,flag_diag)
    use array
    use global
    implicit none
    integer::lda,lwmax,info,lwork,flag_diag
    complex*16::h_temp(2*nos,2*nos)

    allocate(rwork(3*(2*nos)-2),work(2*(2*nos)-1))
    lda=(2*nos)
    lwmax=(2*nos)
    lwork=(2*(2*nos)-1)
    evl_s=0.0d0

    if(flag_diag==1)then
       call zheev('n','u',2*nos,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(flag_diag==2)then
       call zheev('v','u',2*nos,h_temp,lda,evl_s,work,lwork,rwork,info)
    endif

    if(info.ne.0)then
        print*,'algorithm failed'  
    endif 
   
    !  do i=1,2*nos
    !      write(*,*) i, evl_s(i)
    !  enddo
    
    deallocate(rwork,work)  
endsubroutine diagonalization

subroutine auxf(flag,site1)
        use mtmod
        use array
        use global
        integer::flag,site1
        double precision::u1,v1,pi
        pi=4.0*atan(1.0)
        if(flag.eq.0)then
            do i=1,nos    
                call rannum(u1)
                m1=u1
                call rannum(u1)
                ph1=2*pi*u1
                call rannum(v1)
                th1=acos(2*v1-1)
                th(i)=th1
                ph(i)=ph1
                mi(i)=m1
            enddo
        endif
        if(flag.eq.1)then
            call rannum(u1)
            m1=u1
            call rannum(u1)
            ph1=2*pi*u1
            call rannum(v1)
            th1=acos(2*v1-1)
            th(site1)=th1
            ph(site1)=ph1
            mi(site1)=m1
        endif
endsubroutine auxf

subroutine mi_dist
    use array
    use global
    implicit none
    integer::np
    doubleprecision::mi_sum
    double precision::pi,eta,mi_min,mi_max,dist,d_mi,mi_
    pi=4.0*atan(1.0)
    mi_min=0
    mi_max=1
    np=1600
    d_mi=abs(mi_max-mi_min)/np
    d_mi=mi_min
    eta=0.10d0
    do ie=1,np
        mi_=mi_+d_mi
        mi_sum=0.0d0
        do j=1,nos
          mi_sum=mi_sum+((eta/pi)/((dist-mi(j))**2+((eta)**2)))
        enddo
        dist_mi_sum(ie)=dist_mi_sum(ie)+(mi_sum/(nos))
        if(sw.eq.nms) dist1(ie)=mi_
    enddo
endsubroutine mi_dist



subroutine rannum(r)
    use mtmod
    implicit none
    double precision::r
    r=grnd()
    return
endsubroutine rannum

double precision function get_mu_s(fill,temp2)
    use array
    use global
    implicit none
    double precision::f,fL2,fR,mR,mL,m_d,temp2
    double precision::fill

    mR = maxval(evl_s)       !right-side chemical potential
    fr=0.0d0
    do i=1,(2*nos)
        fr=fr+(1.0d0/(exp((evl_s(i)-mR)/Temp2)+1.0d0))
    end do

    mL = minval(evl_s)       !left-side chemical potential
    fL2=0.0d0

    do i=1,(2*nos)
        fL2=fL2+(1.0d0/(exp((evl_s(i)-mL)/Temp2)+1.0d0))
    end do

    m_d = 0.5d0*(mL+mR)    !middle chemical potential
    f=0.0d0
    
    do i=1,(2*nos)
        f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
    end do     
    !print*,f,fill
    do while((abs(f-fill)/(2*nos)).ge.1e-8)
        m_d = 0.5d0*(mL+mR)
         
        f=0.0d0     
        do i=1,(2*nos)
            f=f+(1.0d0/(exp((evl_s(i)-m_d)/Temp2)+1.0d0))
        end do
        if(f.gt.fill)then
            !if middle filling is above target, make it the new right bound.
            mR = m_d
            fR = f
        elseif(f.lt.fill)then
            !if middle filling is below target, make it the new left bound.
            mL = m_d
            fR = f
        endif
    enddo

    !Return the middle value
    get_mu_s = m_d
    return
end function get_mu_s

subroutine n_vs_mu(t11)
    use array
    use global
    implicit none
    integer::ip
    double precision::n_f,t11,m_u,dm_u


    m_u=-3.060d0
    dm_u=0.030d0
    do ip=1,201
       m_u=m_u+dm_u
       n_f=0.0d0
       do i=1,(2*nos)
          n_f=n_f+(1/(1+exp((evl_s(i)-m_u)/t11)))
       enddo
      !write(569,*) m_u,n_f
    enddo
end subroutine n_vs_mu


subroutine energy_cal(temp2)
    use array
    use global
    implicit none
    doubleprecision::temp2,esum
    esum=0.0d0
    do i=1,(2*nos)
        esum=esum+(evl_s(i)/(1.0d0+(exp((evl_s(i)-mu)/temp2))))
    enddo
    !print*,esum
    esum_t=esum+m_sum
  
endsubroutine energy_cal

subroutine dos
    use array
    use global
    implicit none
    double precision::pi,eta_,wmin,wmax,dw,w
    integer::nwp

    pi=4.0*atan(1.0)
    wmin=-10
    wmax=10
    nwp=1000
    dw=abs(wmax-wmin)/nwp
    w=wmin
    eta_=0.10d0

    do ie=1,nwp!no of interval bw bandwitdh of e	
        w=w+dw
        nsum=0.0d0

        do j=1,2*nos
           nsum=nsum+((eta_/pi)/((w-(evl_s(j)))**2+((eta_)**2)))
        enddo
   
        write(200,*) w,nsum/(2*nos)
        flush(200)
    enddo  
           
endsubroutine dos

subroutine struc_factor
    use array
    use global
    implicit none
    double precision::pi,m_x,m_y,m_z,mimj,qx,qy,rx,ry
    pi=4.0*atan(1.0)
    qx=0.0d0
    qy=0.0d0
    !s_f=0.0d0
    sf_A=0.0d0
    
    do i=1,nos ! nos is no of sites  
        if(tag(i).eq.1)then
            do j=1,nos
                if(tag(j).eq.1)then
                    
                    m_x=(mi(i)*sin(th(i))*cos(ph(i)))*(mi(j)*sin(th(j))*cos(ph(j)))
                    m_y=(mi(i)*sin(th(i))*sin(ph(i)))*(mi(j)*sin(th(j))*sin(ph(j)))
                    m_z=(mi(i)*cos(ph(i)))*(mi(j)*cos(ph(j)))
                    mimj=m_x+m_y+m_z
                    rx=qx*(x(i)-x(j))/2
                    ry=qy*(y(i)-y(j))/2
                    !sf_A=sf_A+(complex(cos(rx+ry),sin(rx+ry)))*(mimj)        
                    sf_A=sf_A+exp(cmplx(0.0d0,rx+ry))*(mimj)        

                endif
            enddo
        endif
    enddo
    sf_B_C=0.0d0
    do i=1,nos 
        if((tag(i).eq.2).or.(tag(i).eq.3))then
            do j=1,nos
                if((tag(j).eq.2).or.(tag(j).eq.3))then                  
                    m_x=(mi(i)*sin(th(i))*cos(ph(i)))*(mi(j)*sin(th(j))*cos(ph(j)))
                    m_y=(mi(i)*sin(th(i))*sin(ph(i)))*(mi(j)*sin(th(j))*sin(ph(j)))
                    m_z=(mi(i)*cos(ph(i)))*(mi(j)*cos(ph(j)))
                    mimj=m_x+m_y+m_z
                    rx=qx*(x(i)-x(j))/2
                    ry=qy*(y(i)-y(j))/2
                    !sf_B_C=sf_B_C+(complex(cos(rx+ry),sin(rx+ry)))*(mimj)
                    sf_B_C=sf_B_C+exp(cmplx(0.0d0,rx+ry))*(mimj)        
                    print*,i,j
                endif
            enddo
        endif
    enddo
endsubroutine struc_factor

subroutine spec_func
    use array
    use global
    implicit none
    double precision::mat_e,c1,pi,eta,spec_f,kx,ky
    complex*16::c2
    complex*16::mat_sum(2*nos,2*nos)
    integer::b1,l,k_x,k_y
    pi=4.0*atan(1.0)
    eta=0.025000d0   
    e=2.0d0  
    mat_sum=0.0d0
    do l=1,nos     
        do j=1,nos
            do b1=1,2*nos
                mat_e=((h(2*j,b1)*conjg(h(2*l,b1)))+(h(2*j-1,b1)*conjg(h(2*l-1,b1))))
                mat_sum(l,j)=mat_sum(l,j)+mat_e*((eta/(2*pi))/((e-evl_s(b1))**2+((eta/2)**2)))                     
            enddo
        enddo
    enddo
          
    do k_y=1,(d/2)   
        ky =(2*pi/(d/2))*k_y
        do k_x=1,(d/2)
            kx=(2*pi/(d/2))*k_x
            spec_f=0.0d0 
            do l=1,nos     
                do j=1,nos                                
                    c1=(kx*(x(j)-x(l))+ky*(y(j)-y(l)))
                    c2=complex(cos(c1),sin(c1))
                    spec_f=spec_f+c2*mat_sum(l,j)                      
                enddo   
            enddo    
            write(120,*) kx,ky,spec_f/nos
            flush(120)
        enddo       
    enddo
end subroutine spec_func

subroutine dos_m_c
    use array
    use global
    implicit none
    double precision::pi,eta,wmin,wmax,dw,w
    integer::nwp
    
    pi=4.0*atan(1.0)
    wmin=-14
    wmax=14
    nwp=1600
    dw=abs(wmax-wmin)/nwp
    w=wmin
    eta=0.050d0
    do ie=1,nwp!no of pointss bw bandwitdh of e
        w=w+dw
        nsum=0.0d0
        do j=1,2*nos
          nsum=nsum+((eta/pi)/((w-(evl_s(j)-mu))**2+((eta)**2)))
        enddo
        nsum_t(ie)=nsum_t(ie)+(nsum/(2*nos))
        if(sw.eq.nms) ee1(ie)=w
    enddo
       
endsubroutine dos_m_c

subroutine proj_dos_m_c(temp2)
    use global
    use array
    implicit none
    integer::nwp
    double precision::eta_,wmin,wmax,dw,w,p_dos_a,fermi_fn,p_dos_b,p_dos_c,mat_a,mat_b,mat_c,pi,temp2
    
    pi=4.0*atan(1.0) 
    wmin=-25
    wmax=25
    nwp=2000
    dw=abs(wmax-wmin)/nwp
    w=wmin
    eta_=0.10d0

    do ie=1,nwp
        
        w=w+dw
        p_dos_a=0.0d0
        p_dos_b=0.0d0
        p_dos_c=0.0d0
        
        
        
        do i=1,nos

            if(tag(i).eq.1)then
                do j=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))
                mat_a=((h(2*i,j)*conjg(h(2*i,j)))+(h(2*i-1,j)*conjg(h(2*i-1,j))))
                p_dos_a=p_dos_a+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))*mat_a
                 enddo
            endif
        enddo
        
        dos_sum_a(ie)=dos_sum_a(ie)+(p_dos_a/(2*nos))

        do i=1,nos
        
            if (tag(i).eq.2)then
                do j=1,2*nos
                    fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))
                    mat_b=((h(2*i,j)*conjg(h(2*i,j)))+(h(2*i-1,j)*conjg(h(2*i-1,j))))
                    p_dos_b=p_dos_b+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))*mat_b
                enddo
            endif
        enddo

        dos_sum_b(ie)=dos_sum_b(ie)+(p_dos_b/(2*nos))

        do i=1,nos
            if (tag(i).eq.3)then
                do j=1,2*nos
                    fermi_fn=(1/(1.0d0+(exp((evl_s(j)-mu)/temp2))))
                    mat_c=((h(2*i,j)*conjg(h(2*i,j)))+(h(2*i-1,j)*conjg(h(2*i-1,j))))
                    p_dos_c=p_dos_c+((eta_/pi)/((w-(evl_s(j)-mu))**2+((eta_)**2)))*mat_c
                enddo
            endif
        enddo
        dos_sum_c(ie)=dos_sum_c(ie)+(p_dos_c/(2*nos))

        !write(300,*)w-mu,p_dos_a/(2*nos),p_dos_b/(2*nos),p_dos_c/(2*nos),(p_dos_a+p_dos_b+p_dos_c)/(2*nos)
        
                                    
    enddo 
endsubroutine proj_dos_m_c



subroutine occupation(temp2)
    use array
    use global
    implicit none
    integer::l,k,a,b
    doubleprecision::fermi_fn,temp2

    n_a_up=0.0d0
    n_b_up=0.0d0
    n_c_up=0.0d0
    n_a_dn=0.0d0
    n_b_dn=0.0d0
    n_c_dn=0.0d0

    do l=1,nos
        if(tag(l).eq.1)then
            a=2*l-1
            do k=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(k)-mu)/temp2))))         
                b=k
                n_a_up=n_a_up+(h(a,b)*conjg(h(a,b)))*fermi_fn
                n_a_dn=n_a_dn+(h(a+1,b)*conjg(h(a+1,b)))*fermi_fn
            enddo
        endif
          
        if((tag(l).eq.2))then
            a=2*l-1
            do k=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(k)-mu)/temp2))))
                b=k
                n_b_up=n_b_up+(h(a,b)*conjg(h(a,b)))*fermi_fn
                n_b_dn=n_b_dn+(h(a+1,b)*conjg(h(a+1,b)))*fermi_fn
            enddo
        endif

        if(tag(l).eq.3)then
            a=2*l-1
            do k=1,2*nos
                fermi_fn=(1/(1.0d0+(exp((evl_s(k)-mu)/temp2))))
                b=k
                n_c_up=n_c_up+(h(a,b)*conjg(h(a,b)))*fermi_fn
                n_c_dn=n_c_dn+(h(a+1,b)*conjg(h(a+1,b)))*fermi_fn
            enddo
        endif
    enddo
endsubroutine occupation
