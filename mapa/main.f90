program main
    use rutine


    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654

    integer :: nv, ns, ned, analiza, u
    real(dp) :: fc, ec, fs, es
    real(dp), allocatable :: xy(:,:), ed(:,:), xys(:,:), ds(:)

    real(dp) :: a0,a2



    !Branje podatkov
    open(newunit = u, file = "in.txt",status = "old")
    read(u,*) analiza, fc, ec, fs, es, nv, ns, ned
    close(u)

    if (ned == 0) then
        allocate(xy(2,nv),xys(2,ns),ds(ns))
        open(newunit = u, file = "in.txt",status = "old")
        read(u,*) analiza, fc, ec, fs, es, nv, ns, ned, xy, xys, ds
        close(u)
    else
        allocate(xy(2,nv),xys(2,ns),ed(3,ned),ds(ns))
        open(newunit = u, file = "in.txt",status = "old")
        read(u,*) analiza, fc, ec, fs, es, nv, ns, ned,  xy, xys, ds, ed
        close(u)
    end if



    !RAČUN KARAKTERISTIK IN PREMIK TEZISCA V IZODISCE
    block
        real(dp) :: sx=0.0_dp,sy = 0.0_dp,xyt(2,nv)
        a0 = 0.0_dp
        a2 = 0.0_dp


        xyt(1,:) = xy(2,:)
        xyt(2,:) = xy(1,:)

        call area_n(xy,nv,a0,0)
        call area_n(xy,nv,sx,1)
        call area_n(xyt,nv,sy,1)

        xy(2,:) = xy(2,:)-sx/a0
        xy(1,:) = xy(1,:)+sy/a0
        xys(2,:) = xys(2,:)-sx/a0
        xys(1,:) = xys(1,:)+sy/a0


        call area_n(xy,nv,a2,2)

    end block

    goto (10,20,30), analiza




    !RAČUN ELASTICNEGA
 10 block
        real(dp) :: rd(3,20402),phi, r, def_pl(2),y_extr(2),y_y, xyt(2,nv),i0,i1,i2,h,xyts(2,nv),ys_extr(2)
        real(dp) :: eps_y
        eps_ys = fs/es
        eps_yc = fc/ec


        !print *, eps_y

        do i=0,199
            phi = 2.0_dp*pi*(i/100.0_dp)

            xyt(1,:) = xy(1,:)*cos(phi) + xy(2,:) *sin(phi)
            xyt(2,:) = xy(2,:)*cos(phi) - xy(1,:) *sin(phi)

            xyts(1,:) = xys(1,:)*cos(phi) + xys(2,:) *sin(phi)
            xyts(2,:) = xys(2,:)*cos(phi) - xys(1,:) *sin(phi)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            h = y_extr(2)-y_extr(1)
            ys_extr = (/minval(xyts(2,:)),maxval(xyts(2,:))/)

            !Karakteristike rotiranega prereza
            i0 = 0.0_dp
            i1 = 0.0_dp
            i2 = 0.0_dp

            call area_n(xyt,nv,i0,0)
            call area_n(xyt,nv,i1,1)
            call area_n(xyt,nv,i2,2)

            do j = 0,202

                def_pl(2) = (-eps_yc-eps_ys)/(y_extr(2)-ys_extr(2))*((202-j)/202.0_dp)
                def_pl(1) = -eps_cy-def_pl(2)*y_extr(2)

                rd(:,(i)*200+(j+1)) = -ec*(/def_pl(1)*i0, sin(phi)*(def_pl(2)*i2) , cos(phi)*(def_pl(2)*i2)/)

                 do k=1,ns
                    rd(:,(i)*200+(j+1)) = rd(:,(i)*200+(j+1)) + ds(k)**2*pi/4.0_dp*(es)*(/(def_pl(1)+def_pl(2)*xyts(2,k)), cos(phi)*((def_pl(1)+def_pl(2)*xyts(2,k))*xyts(2,k)) , sin(phi)*(def_pl(1)+def_pl(2)*xyts(2,k))*xyts(2,k)/)
                end do

            end do
        end do


        call write_js(xy,nv,rd,40600,ed,ned,analiza,xys,ds,ns)
    end block
    goto 110


    !RAČUN PLASTICNEGA
20  block
        real(dp) :: y_extr(2), ys_extr(2), y_crac(3), phi, i0, iy1, iy2, iy3, ixy0, ixy1, ixy2
        real(dp) :: xy_eff(2,2*nv), xyt_eff(2,2*nv), xyt(2,nv), xyts(2,ns), rd(3,40600), rot(2,2)
        real(dp) :: def_pl(2), eps_s, eps_y
        eps_y = fs/es

        xy_eff(:,:) = 0.0_dp


        do i=0,199
            phi = 2.0_dp*pi*(i/200.0_dp)
            rot(1,:) = (/cos(phi),sin(phi)/)
            rot(2,:) = (/-sin(phi),cos(phi)/)

            xyt = matmul(rot,xy)
            xyts = matmul(rot,xys)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            ys_extr = (/minval(xyts(2,:)),maxval(xyts(2,:))/)

            !BETON V TLAKU
            do j = 0,20

                def_pl(2) = (-0.0035_dp)/(y_extr(2)-y_extr(1))*(j/20.0_dp)
                def_pl(1) = (-0.002_dp)+def_pl(2)*(-y_extr(2) + (1.0_dp-2.0_dp/3.5_dp)*(y_extr(2)-y_extr(1)))

                !Pomembe y vrednosti
                y_crac(1) = y_extr(1)-1.0_dp; y_crac(3) = y_extr(2)+1.0_dp
                if (def_pl(2) == 0.0_dp) then
                    y_crac(2)= y_crac(3)-0.5_dp
                else
                    y_crac(2)= (-0.002_dp-def_pl(1))/def_pl(2)
                end if



                i0 = 0.0_dp; iy1 = 0.0_dp; iy2 = 0.0_dp; iy3 = 0.0_dp
                ixy0 = 0.0_dp; ixy1 = 0.0_dp; ixy2 = 0.0_dp

                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_crac(2),y_crac(3),xy_eff)

                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,iy1,1)
                call area_ny1x(xy_eff,nv*2,ixy0,0)

                rd(:,(i)*200+(j+1)) = -fc*(/i0, iy1, ixy0/)



                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_crac(1),y_crac(2),xy_eff)

                i0 = 0.0_dp; iy1 = 0.0_dp; ixy0 = 0.0_dp

                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,iy1,1)
                call area_n(xy_eff,nv*2,iy2,2)
                call area_n(xy_eff,nv*2,iy3,3)
                !Integral y^2*x
                call area_ny1x(xy_eff,nv*2,ixy0,0)
                call area_ny1x(xy_eff,nv*2,ixy1,1)
                call area_ny1x(xy_eff,nv*2,ixy2,2)


                rd(:,(i)*200+(j+1)) = rd(:,(i)*200+(j+1)) + fc*(/(def_pl(2)/0.002_dp)**2*iy2 + 2*(1.0_dp+def_pl(1)/0.002_dp)*def_pl(2)/0.002_dp*iy1 - (1-(1+def_pl(1)/0.002_dp)**2)*i0,(def_pl(2)/0.002_dp)**2*iy3 + 2*(1.0_dp+def_pl(1)/0.002_dp)*def_pl(2)/0.002_dp*iy2 - (1-(1+def_pl(1)/0.002_dp)**2)*iy1 ,(def_pl(2)/0.002_dp)**2*ixy2 + 2*(1.0_dp+def_pl(1)/0.002_dp)*def_pl(2)/0.002_dp*ixy1 - (1-(1+def_pl(1)/0.002_dp)**2)*ixy0/)


                do k=1,ns
                    rd(:,(i)*200+(j+1)) = rd(:,(i)*200+(j+1)) + (fs)*(/1.0_dp,xyts(2,k),xyts(1,k)/)*sig_eps_s(def_pl(1) + xyts(2,k)*def_pl(2),eps_y)*ds(k)**2*pi/4.0_dp
                end do

            end do

            !UPOGIB
            do j = 0,180
                !Deformacija najnižje točke betona
                eps_s = 0.055_dp*j/180.0_dp

                def_pl(2) = (-0.0035_dp-eps_s)/(y_extr(2)-y_extr(1))
                def_pl(1) = (-0.0035_dp)- def_pl(2) *y_extr(2)


                !Pomembe y vrednosti
                y_crac(1) = (-def_pl(1))/def_pl(2)
                y_crac(2) = (-0.002_dp-def_pl(1))/def_pl(2)
                y_crac(3)= (-0.0035_dp-def_pl(1))/def_pl(2)

                i0 = 0.0_dp
                iy1 = 0.0_dp
                iy2 = 0.0_dp
                iy3 = 0.0_dp
                ixy0 = 0.0_dp
                ixy1 = 0.0_dp
                ixy2 = 0.0_dp

                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_crac(2),y_crac(3),xy_eff)

                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,iy1,1)
                call area_ny1x(xy_eff,nv*2,ixy0,0)

                rd(:,(i)*200+(21 + j+1)) = -fc*(/i0, iy1 , ixy0/)

                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_crac(1),y_crac(2),xy_eff)

                i0 = 0.0_dp
                iy1 = 0.0_dp
                ixy0 = 0.0_dp

                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,iy1,1)
                call area_n(xy_eff,nv*2,iy2,2)
                call area_n(xy_eff,nv*2,iy3,3)
                !Integral y^2*x
                call area_ny1x(xy_eff,nv*2,ixy0,0)
                call area_ny1x(xy_eff,nv*2,ixy1,1)
                call area_ny1x(xy_eff,nv*2,ixy2,2)



                rd(:,(i)*200+(21+j+1)) = rd(:,(i)*200+(21+j+1)) + fc*(/(def_pl(2)/0.002_dp)**2*iy2 + 2*(1.0_dp+def_pl(1)/0.002_dp)*def_pl(2)/0.002_dp*iy1 - (1-(1+def_pl(1)/0.002_dp)**2)*i0,(def_pl(2)/0.002_dp)**2*iy3 + 2*(1.0_dp+def_pl(1)/0.002_dp)*def_pl(2)/0.002_dp*iy2 - (1-(1+def_pl(1)/0.002_dp)**2)*iy1 , (def_pl(2)/0.002_dp)**2*ixy2 + 2*(1.0_dp+def_pl(1)/0.002_dp)*def_pl(2)/0.002_dp*ixy1 - (1-(1+def_pl(1)/0.002_dp)**2)*ixy0/)

                do k=1,ns
                    rd(:,(i)*200+(21+j+1)) = rd(:,(i)*200+(21+j+1)) +(fs)*(/1.0_dp,xyts(2,k),xyts(1,k)/)*sig_eps_s(def_pl(1) + xyts(2,k)*def_pl(2),eps_y)*ds(k)**2*pi/4.0_dp
                end do

            end do
            rd(2:3, i*200 + 1 : i*200 +203) =matmul(rot,rd(2:3, i*200 + 1 : i*200 +203))
        end do

        call write_js(xy,nv,rd,40600,ed,ned,analiza,xys,ds,ns)
    end block
goto 110




30  block
        real(dp) :: rd(6,20401),phi, r, def_pl(2), y_extr(2), y_y, xyt(2,nv), i0, i1, i2, xy_eff(2,2*nv)
        real(dp) :: i0el, i1el, i2el, h, eps_y

        eps_y = fy/es
        xy_eff(:,:) = 0.0_dp


        do i=0,100
            phi = 2.0_dp*pi*(i/100.0_dp)
            !print *, phi
            xyt(1,:) = xy(1,:)*cos(phi) + xy(2,:) *sin(phi)
            xyt(2,:) = xy(2,:)*cos(phi) - xy(1,:) *sin(phi)

            y_extr = (/minval(xyt(2,:)),maxval(xyt(2,:))/)
            h = y_extr(2)-y_extr(1)

            i0el = 0.0_dp
            i1el = 0.0_dp
            i2el = 0.0_dp

            call area_n(xyt,nv,i0el,0)
            call area_n(xyt,nv,i1el,1)
            call area_n(xyt,nv,i2el,2)


            do j = 0,100


            end do


            do j = 0,100

                !ELASTICNA
                def_pl(2) = 2.0_dp*eps_y/h*((100-j)/100.0_dp)
                def_pl(1) = eps_y-def_pl(2)*y_extr(2)
                rd(1:3,(i)*100+(j+1)) = es*(/def_pl(1)*i0el+def_pl(2)*i1el, sin(phi)*(def_pl(1)*i1el+def_pl(2)*i2el) , cos(phi)*( def_pl(1)*i1el+def_pl(2)*i2el)/)
                rd(1:3,10201+(i)*100+(j+1)) = -rd(1:3,(i)*100+(j+1))



                !PLASTICNA
                r = y_extr(1)*(j/100.0_dp) + y_extr(2)*((100-j)/100.0_dp)

                i0 = 0.0_dp
                i1 = 0.0_dp
                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,r,y_extr(2),xy_eff)
                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)
                rd(4:6,(i)*100+(j+1)) = fy*(/i0, cos(phi)*(i1) , sin(phi)*( i1)/)


                xy_eff(:,:) = 0.0_dp
                call eff_section(xyt,nv,y_extr(1),r,xy_eff)
                i0 = 0.0_dp
                i1 = 0.0_dp
                call area_n(xy_eff,nv*2,i0,0)
                call area_n(xy_eff,nv*2,i1,1)

                rd(4:6,(i)*100+(j+1)) = rd(4:6,(i)*100+(j+1)) - fy*(/i0, cos(phi)*(i1) , sin(phi)*( i1)/)
                !d(4:6,10201+(i)*100+(j+1)) = -rd(4:6,(i)*100+(j+1))
            end do
        end do

        call write_js(xy,nv,rd,20401,ed,ned,analiza,xys,ds,ns)
    end block
goto 110


110     print *," "
        print*,"    Račun končan."
        !print*," "








end program
