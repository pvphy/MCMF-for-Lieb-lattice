module array
    implicit none
    double precision, allocatable::mi(:),th(:),ph(:),rwork(:),evl_s(:)
    integer,dimension(:),allocatable::x(:),y(:),tag(:),x_cl,y_cl,tag_cl
    complex*16,allocatable:: h(:,:),bb(:,:),work(:)
    double precision::nsum_t(1600),nn1(1600),ee1(1600),e_2(1600)
    double precision, allocatable::mi_cl(:),th_cl(:),ph_cl(:)
    doubleprecision::dos_sum_a(1600),dos_sum_b(1600),dos_sum_c(1600),dist1(1600),dist_mi_sum(1600)
end module array

module global
    implicit none
    double precision::t,ua,ub,uc,mu,m_sum,e_sum,u,esum_t,m1,ph1,th1,nsum,e,s_f,sf_A,sf_B_C
    double precision::lambda,n_a,n_b,n_c,n_a_up,n_b_up,n_c_up,n_a_dn,n_b_dn,n_c_dn,filling_cl
    integer::d,nos,ie,flag_isoc,unit_cells,sw,nms
    integer::unit_cells_cl,nos_cl,dim_cl,d_cl
endmodule global

program liebmain
    use array
    use global
    implicit none
    integer::seed,it,i,j
    double precision::temp,temp1
    double precision::m_temp1,th_temp1,ph_temp1,avg_sf_A,avg_sf_B_C,inpt,filling,fill_inpt
 
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
    
    unit_cells=(d/2)**2
    nos=3*unit_cells ! no of sites
    filling=fill_inpt*2*nos
      
    print*,'unit cells_____________=',unit_cells
    print*,'sites system___________=',nos
    print*,'hopping parameter______=',t
    print*,'Ua_____________________=',ua
    print*,'ub_____________________=',ub
    print*,'uc_____________________=',uc
    print*,'Monte-Carlo steps______=',nms
    print*,'filling system_________=',filling

    allocate(x(nos),y(nos),h((2*nos),(2*nos)),mi(nos),th(nos),ph(nos),tag(nos))

    call lattice_labeling  
    Temp=0.650d0

    
    do iT=1,30
        open(1000+it)   
        if(iT.le.11) Temp=Temp-0.050d0
        if ((iT.gt.11).and.(iT.lt.21)) Temp=Temp-0.010d0
        if ((iT.gt.21)) Temp=Temp-0.0010d0

        avg_sf_A=0.0d0
        avg_sf_B_C=0.0d0
        dist_mi_sum=0.0d0
        do sw=1,nms/10
            j=0
            do i=(sw*nos-nos+1),sw*nos
                j=j+1
                read(1000+it,*) Temp1,m_temp1,th_temp1,ph_temp1 
                mi(j)=m_temp1
                th(j)=th_temp1
                ph(j)=ph_temp1
            enddo


            call mi_dist
            call struc_factor
            avg_sf_A=avg_sf_A+(sf_A/nos**2)
            avg_sf_B_C=avg_sf_B_C+(sf_B_C/nos**2)    
            
        enddo 
        
        write(93,*) temp,avg_sf_A/(nms/2),avg_sf_B_C/(nms/2)
        flush(93)
        do ie=1,1600  
            write(700+it,*) dist1(ie),dist_mi_sum(ie)/(nms/10)
            flush(700+it)
        enddo
    enddo   

end program liebmain

subroutine lattice_labeling
    use global
    use array
    implicit none
    integer::ix,iy,i1,i
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


subroutine mi_dist
    use array
    use global
    implicit none
    integer::np,j
    doubleprecision::mi_sum
    double precision::pi,eta,mi_min,mi_max,d_mi,mi_
    pi=4.0*atan(1.0)
    mi_min=0.0d0
    mi_max=1.0d0
    np=1600
    d_mi=abs(mi_max-mi_min)/np

    mi_=mi_min
    eta=0.10d0
    do ie=1,np
        mi_=mi_+d_mi
        mi_sum=0.0d0
        do j=1,nos
          mi_sum=mi_sum+((eta/pi)/((mi_-mi(j))**2+((eta)**2)))
        enddo
        dist_mi_sum(ie)=dist_mi_sum(ie)+(mi_sum/(nos))
        if(sw.eq.(nms/10)) dist1(ie)=mi_
    enddo
endsubroutine mi_dist

subroutine struc_factor
    use array
    use global
    implicit none
    double precision::pi,m_x,m_y,m_z,mimj,qx,qy,rx,ry
    integer::i,j
    pi=4.0*atan(1.0)
    qx=0.0d0
    qy=0.0d0
    !s_f=0.0d0
    sf_A=0.0d0
    sf_B_C=0.0d0
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

                endif
            enddo
        endif
    enddo
endsubroutine struc_factor

