program zemi
    implicit none
    real   ,parameter::g=9.8,pi=3.141592,ng=0.025,roh=1020.0      ! 重力加速度，円周率，マニングの粗度係数，海水の密度
    real   ,parameter::dx=10,dt=0.005,amp=1.0,wp=50               ! 刻み値
    integer,parameter::nx=200,nt=80000,np=200                     ! ステップ数

    real m1(nx), m2(nx), e1(nx), e2(nx), zb(nx), D(nx), D1(nx), D2(nx),time
    real rh, ad, fr
    integer i, n

    m1=0.0;m2=0.0;e1=0.0;e2=0.0;zb=0.0;D=0.0

    do i=1, nx  !各地点での初期水位の設定（水面をz=0mとする）
        zb(i)=-10.0+0.05*i
    end do

    D=-zb

    open(10,file='e.dat')
    open(11,file='m.dat')

    do n=1, nt  !時間発展

        time=n*dt   !造波境界の設定
        m1(1)=amp*sin(2.*pi* time    /wp)
        m2(1)=amp*sin(2.*pi*(time+dt)/wp)
        
        do i=2,nx-1   !変動水位の計算
            rh=-(m1(i)-m1(i-1))/dx
            e2(i)=rh*dt+e1(i)
        end do

        do i=2,nx-2   !線流量の計算
            if (i.gt.3.and.i.le.nx-3) then  !移流項の計算
                ad = 0.0
                if (m1(i).gt.0) then
                    D1(i) = (D(i+1)+D(i))/2.
                    D2(i) = (D(i)+D(i-1))/2.
                    ad = -(m1(i)**2/D1(i)-m1(i-1)**2/D2(i))/dx
                else
                    D1(i) = (D(i+2)+D(i+1))/2.
                    D2(i) = (D(i+1)+D(i))/2.
                    ad = -(m1(i+1)**2/D1(i)-m1(i)**2/D2(i))/dx
                end if
            end if

            fr = -g*(ng**2)*(m1(i)**2.)/(D(i)**2.33)    !摩擦項の計算

            rh=(-g*D(i)*(e1(i+1)-e1(i))/dx)+ad+fr     !右辺の計算
            m2(i)=rh*dt+m1(i)
        end do

        m1=m2   !次の時間ステップの初期値
        e1=e2   !次の時間ステップの初期値

        if (mod(n,np).eq.1) then    !全時間ステップから200ずつ出力
            write(10,'(100e13.3)') e2(2:nx-1)
            write(11,'(100e13.3)') m2(1:nx-1)
        end if

    end do

    close(10)
    close(11)

    stop

end program zemi
