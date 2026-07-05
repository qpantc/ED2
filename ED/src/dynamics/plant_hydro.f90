!==========================================================================================!
!==========================================================================================!
! MODULE: PLANT_HYDRO
!! 驱动植物水动力学模型。利用气孔开度，解算水分如何从各层土壤通过根系、木质部管道源源不断地泵向叶片。
!! 根据当前的土壤水分状态和植物上一时间步的蒸腾操作，计算植物体内（叶片、木质部）的水势（Water Potential）以及
!! 土壤与植物之间的水分通量（Water Flux）
!! 这段代码实现的是经典的 土壤-植物-大气连续体（SPAC, Soil-Plant-Atmosphere Continuum） 理论。
!! 植物就像一根多孔的管道，水从土壤流向根部，通过木质部（wood），最后从叶片（leaf）蒸腾到大气中。
!! 水流动的驱动力是水势差（$\psi$）。
!> \brief Calculations of plant hydrodynamics at ED2 timestep, including various
!> utils for plant hydraulics
!> \details Util functions to perform conversions between psi, rwc and water_int
!> \author Xiangtao Xu, 26 Jan. 2018
!==========================================================================================!
!==========================================================================================!
module plant_hydro

   !------ Tolerance for minimum checks.  -------------------------------------------------!
   real(kind=4), parameter :: tol_buff   = 1.e-4
   real(kind=8), parameter :: tol_buff_d = dble(tol_buff)
   ! 解析： 定义单精度（kind=4）和双精度（kind=8，通过 dble 转换）的容差缓冲值（$10^{-4}$），用于后续数学计算的边界防御，防止数值下溢或除以 0。

   real(kind=4), parameter :: om_buff   = 1. - tol_buff
   real(kind=8), parameter :: om_buff_d = dble(om_buff)

   real(kind=4), parameter :: op_buff   = 1. + tol_buff
   real(kind=8), parameter :: op_buff_d = dble(op_buff)
   ! 解析： 定义略小于 1 (om_buff = 0.9999) 和略大于 1 (op_buff = 1.0001) 的乘数因子，通常用于限制变量的数值范围。

   real(kind=4), parameter :: mg_safe   = 0.05
   real(kind=4), parameter :: op_safe   = 1. + mg_safe
   real(kind=4), parameter :: om_safe   = 1. - mg_safe
   ! 解析： 定义数值积分的安全权重系数。mg_safe（5%）与 om_safe（95%）在后续的土壤水分截断中起到了加权缓冲作用。
   !---------------------------------------------------------------------------------------!

   contains
   ! 解析： 标志着模块的声明部分结束，下面开始定义具体的函数和子程序。
    
   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE PLANT_HYDRO_DRIVER
   !> \brief   Main driver to calculate plant hydrodynamics within a site.
   !> \details This subroutine works at DTLSM scale, similar to photosyn_driv.
   !> Keep in mind that this subroutine uses average transpiration from last timestep.
   !> Therefore, we should also use water potential at the start of the last
   !> timestep. \n
   !>   Alternatively, one can directly call calc_plant_water_flux
   !> subroutine within the rk4_derivs.f90, which will give an estimate of water
   !> flux within each integration timestep using the transpiration/water
   !> potential at the start of the integration timestep. However, this can
   !> incur extra computational cost.
   !>
   !> \author Xiangtao Xu, 30 Jan. 2018
   !---------------------------------------------------------------------------------------!
   subroutine plant_hydro_driver(csite,ipa,ntext_soil)
      ! 解析： 声明驱动子程序，接收三个核心参数：当前站点结构体（csite）、斑块索引（ipa）和土壤纹理类型数组（ntext_soil）。

      use ed_state_vars        , only : sitetype               & ! structure
                                      , patchtype              ! ! structure
      ! 解析： 从状态变量模块引用 sitetype（站点类型）和 patchtype（斑块类型）的结构体定义。
      use ed_misc_coms         , only : dtlsm                  & ! intent(in)
                                      , dtlsm_o_frqsum         & ! intent(in)
                                      , current_time           ! ! intent(in)
      ! 解析： 引用时间相关变量：dtlsm 是陆面过程的时间步长；
      ! dtlsm_o_frqsum 是当前步长占总输出频率时间的权重比例；
      ! current_time 包含当前模拟的绝对时间（年月日）。

      use soil_coms            , only : soil                   & ! intent(in)
                                      , matric_potential       & ! function
                                      , hydr_conduct           ! ! function
      ! 解析： 引用土壤物理模块：soil 存储土壤基础属性；
      ! matric_potential 函数用于计算土壤基质势；
      ! hydr_conduct 函数用于计算土壤导水率。
      use grid_coms            , only : nzg                    ! ! intent(in)
      use consts_coms          , only : pio4                   & ! intent(in)
                                      , wdns                   ! ! intent(in)
      use allometry            , only : dbh2sf                 ! ! function
      use physiology_coms      , only : plant_hydro_scheme     ! ! intent(in)
      ! 解析： * nzg：整个垂直剖面的土壤层数。pio4：几何常量 $\pi / 4$（用于计算圆面积）。
      ! wdns：水的密度（常数）。
      ! dbh2sf：异速生长方程函数，根据胸径（DBH）计算边材比例。
      ! plant_hydro_scheme：全局物理方案开关选项。
      use pft_coms             , only : C2B                    & ! intent(in)
                                      , leaf_water_cap         & ! intent(in)
                                      , leaf_psi_min           & ! intent(in)
                                      , small_psi_min          ! ! intent(in)
      ! 解析： 引用植物功能型（PFT）生理参数：C2B（碳到生物量的转换系数）；
      ! leaf_water_cap（叶片水容/热容）；
      ! leaf_psi_min 与 small_psi_min 分别表示成年树与幼树的极限凋萎水势。

      implicit none
      ! 解析： 强制显式声明所有变量类型，防止 Fortran 的隐式类型错误，是标准规范。
      ! 本地变量声明区 (Local Variables)
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      integer               , intent(in)  :: ipa
      integer,dimension(nzg), intent(in)  :: ntext_soil
      ! 解析： 定义输入参数的类型。csite 作为 target，意味着它允许被内部的指针指向。
      !----- Local Vars  ------------------------------------------------------------------!
      type(patchtype)       , pointer     :: cpatch      !< patch strcture
      real                                :: swater_min  !< Min. soil moisture for condct.
      real                                :: swater_max  !< Max. soil moisture for condct.
      real                                :: swater_use  !< soil moisture
      integer                             :: nsoil       !< soil type for soil
      integer                             :: k           !< iterator for soil lyr
      integer                             :: ico         !< iterator for cohort
      integer                             :: ipft        !< PFT index
      real ,dimension(nzg)                :: soil_psi    !< soil water potential   [      m]
      real ,dimension(nzg)                :: soil_cond   !< soil water conductance [kg/m2/s]
      real                                :: sap_frac    !< sapwood fraction       [    ---]
      real                                :: sap_area    !< sapwood area           [     m2]
      real                                :: bsap        !< sapwood biomass        [    kgC]
      real                                :: transp      !< transpiration rate     [   kg/s]
      real                                :: c_leaf      !< leaf capacitance       [   kg/m]
      logical                             :: track_hydraulics !< whether track hydraulics
      ! 解析： 声明内部局部变量，主要包括单层土壤的水分约束范围（swater_min/max/use）、土壤层循环因子 k、
      ! 植物群落循环因子 ico，以及临时存放计算所得的水势（soil_psi）、导水率（soil_cond）、边材面积（sap_area）和
      ! 蒸腾速率（transp）的变量。
      !----- Variables for debugging purposes ---------------------------------------------!
      integer, parameter                  :: dco        = 0 ! the cohort to debug
      logical, dimension(3)               :: error_flag
      logical, parameter                  :: debug_flag = .false.
      character(len=13)     , parameter   :: efmt       = '(a,1x,es12.5)'
      character(len=9)      , parameter   :: ifmt       = '(a,1x,i5)'
      character(len=9)      , parameter   :: lfmt       = '(a,1x,l1)'
      !----- External functions. ----------------------------------------------------------!
      logical               , external    :: isnan_real
      ! 解析： Debug 专用变量和格式化输出字符串。isnan_real 是一个用于判断数值是否为 NaN（非法无效值）的外部函数。
      !------------------------------------------------------------------------------------!


      !-- Point to the cohort structures --------------------------------------------------!
      cpatch => csite%patch(ipa)
      ! 解析： 将局部指针 cpatch 指向当前正在计算的土地斑块。这样在后续代码中可以直接使用 cpatch%，省去了冗长的多级结构体寻址，能大幅提高模型执行速度。
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Decide whether or not to solve dynamic plant hydraulics.                      !
      !------------------------------------------------------------------------------------!
      select case (plant_hydro_scheme)
      case (0)
         !------ Compatible with original ED-2.2, do not track plant hydraulics. ----------!
         !! （传统模式）： 关闭动态水动力学。将所有水分通量（wflux）设为 0，
         !! 叶片和木质部的相对含水量（rwc）设为 1.0（永远饱和）。
         do ico = 1, cpatch%ncohorts
             ipft = cpatch%pft(ico)

             cpatch%wflux_wl        (ico)    = 0.
             cpatch%wflux_gw        (ico)    = 0.
             cpatch%wflux_gw_layer(:,ico)    = 0.

             cpatch%leaf_rwc        (ico)    = 1.0
             cpatch%wood_rwc        (ico)    = 1.0
             cpatch%leaf_psi        (ico)    = 0.
             cpatch%wood_psi        (ico)    = 0.
         end do
        !----------------------------------------------------------------------------------!
      case default
         !! （动态模式）： 激活水动力学方案，执行后续的复杂计算。
         !---------------------------------------------------------------------------------!
         !    Dynamic plant hydraulics.                                                    !
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Calculate water potential and conductance in each soil layer in preparation
         ! for later calculations.
         !---------------------------------------------------------------------------------!
         do k = 1,nzg
            !! 遍历了所有土壤层（$nzg$），对于每一层土壤，首先根据土壤类型（ntext_soil(k)）获取该层的土壤参数（nsoil）。
            nsoil = ntext_soil(k)

            !------------------------------------------------------------------------------!
            !      Get bounded soil moisture.                                              !
            !  MLO.  The lower bound used to be air-dry soil moisture.  This causes issues !
            !  in the RK4 integrator if the soil moisture is just slightly above air-dry   !
            !  and dtlsm is long.  For the time being, I am assuming that soil             !
            !  conductivity is halted just below the permanent wilting point.  Similarly,  !
            !  I am assuming that matric potential cannot exceed a value slightly less     !
            !  than the bubbling point.                                                    !
            !------------------------------------------------------------------------------!
            swater_min = mg_safe * soil(nsoil)%soilcp  + om_safe * soil(nsoil)%soilwp
            swater_max = mg_safe * soil(nsoil)%sfldcap + om_safe * soil(nsoil)%slmsts
            ! 解析： 数值物理边界防御。 为防止数值积分器（如后文提到的 RK4）在土壤极端干旱或极度暴雨饱和时崩溃，
            ! 这里计算了安全的上下限：
            ! 最低有效水分 swater_min 介于土壤空气干燥点（soilcp）与永久凋萎点（soilwp）之间。
            ! 最高有效水分 swater_max 介于田间持水量（sfldcap）与饱和含水量（slmsts）之间。
            swater_use = max( swater_min                                                   &
                            , min(swater_max                                               &
                                 ,csite%soil_water(k,ipa) * csite%soil_fracliq(k,ipa) ) )
            ! 解析： 提取当前层真实的液态有效水分（总水分 soil_water $\times$ 液态水比例 soil_fracliq），
            ! 并使用 min 和 max 函数将其死死强制截断在上面算出的 [swater_min, swater_max] 安全区间内，赋值给 swater_use。
            !------------------------------------------------------------------------------!


            !----- Clapp & Hornberger curves. ---------------------------------------------!
            soil_psi(k)  = matric_potential(nsoil,swater_use)
            ! 解析： 调用经验公式函数（如 Clapp-Hornberger 曲线），
            ! 根据截断后的安全水分计算当前层土壤的基质势（soil_psi，单位通常为米水柱）。
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    In the model, soil can't get drier than residual soil moisture.  Ensure   !
            ! that hydraulic conductivity is effectively zero in case soil moisture        !
            ! reaches this level or drier.                                                 !
            !------------------------------------------------------------------------------!
            if (csite%soil_water(k,ipa) < swater_min) then
            ! 如果当前层真实的土壤水分已经低于 swater_min，意味着已经干旱到断流，直接令土壤导水率 soil_cond(k) = 0.；
               soil_cond(k) = 0.
            else
            ! 否则调用 hydr_conduct 计算当前水分下的土壤导水率，并乘以密度常量 wdns 转换单位。
               soil_cond(k) = wdns * hydr_conduct(k,nsoil,csite%soil_water(k,ipa)          &
                                                 ,csite%soil_fracliq(k,ipa))
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Loop over cohorts, calculate plant hydraulic fluxes.                       !
         !---------------------------------------------------------------------------------!
         ! B. 植被群落循环与水势状态解耦。对于每个植物群落（cohort），根据其叶片和木质部的可解析状态，决定是否需要追踪水动力学。
         cohortloop: do ico = 1, cpatch%ncohorts
            ipft = cpatch%pft(ico)
            ! 解析： 开始遍历当前斑块下的每一个植物群落（Cohort）。

            !------------------------------------------------------------------------------!
            !     Track the plant hydraulics when either leaf or wood are resolvable. Leaf !
            ! become un-resolvable in dry season when all leaves are shed. In this scena-  !
            ! rio, we still nedd to track plant hydraulics to update wood_psi so that the  !
            ! model knows when to reflush the leaves. Otherwise, soil water will never re- !
            ! fill wood water pool.                                                        !
            !     Special case: leaf is not resolvable when bleaf is on allometry while    !
            ! wood can be still resolvable. Transpiration is always zero, whereas soil     !
            ! water can still flow into the stem. This will ultimately leads to unreal-    !
            ! istically high wood_psi and even positive psi due to numerical erros. There- !
            ! fore, we do not track plant hydraulics in this case.                         !
            !------------------------------------------------------------------------------!

            track_hydraulics = cpatch%leaf_resolvable(ico) .or.                            &
                              (cpatch%wood_resolvable(ico) .and.                           & 
                               .not. (cpatch%elongf(ico) == 1.0 .and.                      &
                                      .not. cpatch%leaf_resolvable(ico) ))
            ! 解析： 判定当前群落是否需要追踪水动力学。 * 满足以下条件之一即可追踪：
            ! 1. 叶片可解析（即树上有叶子）；
            ! 2. 木质部可解析（即树干存在），且剔除掉“处于异速生长极端状态且无叶”的特殊数值异常。
            ! 这样处理的科学逻辑在于：
            ! 干旱季节大树掉光了叶子（叶片不可解析），但仍要追踪树干的水势，以便雨季来临时树干能吸水，触发模型后续的“春回长叶（reflush）”信号。


            if (track_hydraulics) then
               !----- Prepare input for plant water flux calculations. --------------------!
               sap_frac    = dbh2sf(cpatch%dbh(ico),ipft)                    ! m2
               sap_area    = sap_frac * pio4 * (cpatch%dbh(ico) / 100.) ** 2 ! m2
               bsap        = ( cpatch%bdeada   (ico) + cpatch%bdeadb   (ico)               &
                             + cpatch%bsapwooda(ico) + cpatch%bsapwoodb(ico) ) * sap_frac
               ! 解析： 为当前群落计算结构和形态参数。利用胸径 dbh 算出边材比例 sap_frac，
               ! 进而算出大树截面的边材面积 sap_area。
               !接着把各类枯木和活木生物量相加，乘上边材比例，得到具备传水能力的边材总生物量 bsap。
               transp      = ( cpatch%fs_open(ico) * cpatch%psi_open(ico)                  &
                             + (1. - cpatch%fs_open(ico)) * cpatch%psi_closed(ico) )       &
                           * cpatch%lai(ico) / cpatch%nplant(ico)            ! kg / s
               ! 计算群落的总蒸腾速率（transp）。基于气孔开（fs_open）、闭状态下的蒸腾通量加权平均，
               ! 再乘以叶面积指数 lai，并除以单株植物密度 nplant，转换为单位时间内单株大树的蒸腾失水速率（$kg/s$）。
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Note  that the current leaf_water_int has already deducted losses     !
               ! through transpiration (but not sapflow).  Consequently, leaf psi can be   !
               ! very low.  To get meaningful leaf_psi, leaf_psi and leaf_water_int        !
               ! become temporarily decoupled: transpiration is added back to              !
               ! leaf_water_int.  Therefore, leaf_psi represents the water potential at    !
               ! the START of the timestep.                                                !
               !---------------------------------------------------------------------------!
               call rwc2psi(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico),ipft                 &
                           ,cpatch%leaf_psi(ico),cpatch%wood_psi(ico))
               ! 解析： 调用函数 rwc2psi，将当前的叶片和木质部相对含水量（rwc）转换为初始的瞬时水势（psi）。
               c_leaf = leaf_water_cap(ipft) * C2B * cpatch%bleaf(ico)
               if (c_leaf > 0.) then
                  cpatch%leaf_psi(ico) = cpatch%leaf_psi(ico)  & ! m
                                       + transp * dtlsm        & ! kgH2O
                                       / c_leaf                ! ! kgH2O/m
                  !------------------------------------------------------------------------!
               else
                  !----- No leaves, set leaf_psi the same as wood_psi - hite. -------------!
                  cpatch%leaf_psi(ico) = cpatch%wood_psi(ico) - cpatch%hite(ico)
                  !------------------------------------------------------------------------!
               end if
               ! 解析： 核心数学处理（时间解耦）。 
               ! * 物理上，模型在之前的流水线中已经把当前步的蒸腾失水给扣掉了，直接导致计算出的叶片水势偏低。
               ! * 为了让后续的水分通量方程拿到当前时间步开始（START）时的基准水势，必须在这里进行解耦反推。
               ! * 算出叶片总水容 c_leaf（单位叶片水容 $\times$ 转换常数 $\times$ 叶片生物量 bleaf）。
               ! * 如果叶片水容大于 0，则通过 $\Delta \psi = \frac{\text{transp} \times \text{dtlsm}}{c\_leaf}$，
               ! * 将当前步消耗的水量临时加回到叶片水势上。如果无叶（c_leaf <= 0），
               ! * 则根据水力静压平衡，直接令叶片水势等于木质部水势减去树高 hite。
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Run sanity check.  The code will crash if any of these happen.       !
               !                                                                           !
               ! 1.  If leaf_psi is invalid (run the debugger, the problem may be else-    !
               !     where)                                                                !
               ! 2.  If leaf_psi is positive (non-sensical)                                !
               ! 3.  If leaf_psi is too negative (also non-sensical)                       !
               !---------------------------------------------------------------------------!
               ! C. 生态物理合理性检查 (Sanity Check Block)。对于计算得到的叶片水势（leaf_psi），进行一系列合理性检查：
               error_flag(1) = isnan_real(cpatch%leaf_psi(ico)) ! NaN values
               ! 检查是否产生了无效悬空数字 NaN。
               error_flag(2) = cpatch%leaf_psi(ico) > 0.        ! Positive potential
               ! 检查水势是否变成正数（植物体内由于张力，水势必然为负，正数说明计算彻底失真）。
               error_flag(3) = merge( cpatch%leaf_psi(ico) < small_psi_min(ipft)           &
                                    , cpatch%leaf_psi(ico) < leaf_psi_min (ipft)           &
                                    , cpatch%is_small(ico)                        )
               ! 检查水势是否过于负值，超过了幼树或成年树的极限凋萎水势（根据小树标志 is_small 选择对应的阈值）。
               if ((debug_flag .and. (dco == 0 .or. ico == dco)) .or. any(error_flag)) then
                  write (unit=*,fmt='(a)') ' '
                  write (unit=*,fmt='(92a)') ('=',k=1,92)
                  write (unit=*,fmt='(92a)') ('=',k=1,92)
                  write (unit=*,fmt='(a)'  )                                               &
                     ' Invalid leaf_psi detected.'
                  write (unit=*,fmt='(92a)') ('-',k=1,92)
                  write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : '    &
                                                  ,current_time%year,current_time%month    &
                                                  ,current_time%date,current_time%time
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa
                  write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
                  write (unit=*,fmt=ifmt   ) ' + PFT              =',ipft
                  write (unit=*,fmt=ifmt   ) ' + KRDEPTH          =',cpatch%krdepth(ico)
                  write (unit=*,fmt=efmt   ) ' + HEIGHT           =',cpatch%hite(ico)
                  write (unit=*,fmt=lfmt   ) ' + SMALL            =',cpatch%is_small(ico)

                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=lfmt   ) ' + FINITE           =',.not. error_flag(1)
                  write (unit=*,fmt=lfmt   ) ' + NEGATIVE         =',.not. error_flag(2)
                  write (unit=*,fmt=lfmt   ) ' + BOUNDED          =',.not. error_flag(3)

                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + LEAF_PSI_MIN     =',leaf_psi_min (ipft)
                  write (unit=*,fmt=efmt   ) ' + SMALL_PSI_MIN    =',small_psi_min(ipft)

                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + BLEAF            =',cpatch%bleaf(ico)
                  write (unit=*,fmt=efmt   ) ' + LAI              =',cpatch%lai(ico) 
                  write (unit=*,fmt=efmt   ) ' + NPLANT           =',cpatch%nplant(ico) 
                  write (unit=*,fmt=efmt   ) ' + BSAPWOOD (Hydro) =',bsap 
                  write (unit=*,fmt=efmt   ) ' + BSAPWOOD (Allom) ='                       &
                                                                  , cpatch%bsapwooda(ico)  &
                                                                  + cpatch%bsapwoodb(ico)

                  write (unit=*,fmt=efmt   ) ' + BROOT            =',cpatch%broot(ico)
                  write (unit=*,fmt=efmt   ) ' + SAPWOOD_AREA     =',sap_area
  
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + TRANSP           =',transp
                  write (unit=*,fmt=efmt   ) ' + C_LEAF           =',c_leaf
                  write (unit=*,fmt=efmt   ) ' + PSI_OPEN         =',cpatch%psi_open(ico)
                  write (unit=*,fmt=efmt   ) ' + PSI_CLOSED       =',cpatch%psi_closed(ico)
                  write (unit=*,fmt=efmt   ) ' + FS_OPEN          =',cpatch%fs_open(ico)
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + LEAF_PSI         =',cpatch%leaf_psi(ico)
                  write (unit=*,fmt=efmt   ) ' + LEAF_RWC         =',cpatch%leaf_rwc(ico)
                  write (unit=*,fmt=efmt   ) ' + LEAF_WATER_INT   ='                       &
                                                               ,cpatch%leaf_water_int(ico)
                  write (unit=*,fmt=efmt   ) ' + LEAF_WATER_IM2   ='                       &
                                                               ,cpatch%leaf_water_im2(ico)
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + WOOD_PSI         =',cpatch%wood_psi(ico)
                  write (unit=*,fmt=efmt   ) ' + WOOD_RWC         =',cpatch%wood_rwc(ico)
                  write (unit=*,fmt=efmt   ) ' + WOOD_WATER_INT   ='                       &
                                                               ,cpatch%wood_water_int(ico)
                  write (unit=*,fmt=efmt   ) ' + WOOD_WATER_IM2   ='                       &
                                                               ,cpatch%wood_water_im2(ico)
                  write (unit=*,fmt='(a)'  ) ' '
                  write (unit=*,fmt=efmt   ) ' + WFLUX_GW (LAST)  =',cpatch%wflux_gw(ico) 
                  write (unit=*,fmt=efmt   ) ' + WFLUX_WL (LAST)  =',cpatch%wflux_wl(ico)


                  write (unit=*,fmt='(a)'        ) ' '
                  write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
                  write (unit=*,fmt='(a,2(1x,a))') '    K','    SOIL_PSI','WFLUX_GW_LYR'
                  write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
                  do k = 1, nzg
                     write (unit=*,fmt='(i5,2(1x,es12.5))')                                &
                                                k,soil_psi(k),cpatch%wflux_gw_layer(k,ico)
                  end do
                  write (unit=*,fmt='(92a)'   ) ('-',k=1,92)
                  write (unit=*,fmt='(a)'     ) ' '
                  write (unit=*,fmt='(92a)'   ) ('=',k=1,92)
                  write (unit=*,fmt='(92a)'   ) ('=',k=1,92)
                  write (unit=*,fmt='(a)'     ) ' '

                  if (any(error_flag)) then 
                     call fatal_error('Plant Hydrodynamics is off-track.'                  &
                                     ,'plant_hydro_driver','plant_hydro.f90')
                  end if
                  !------------------------------------------------------------------------!
               end if
               ! 解析： 如果触发了上述任意一项物理红线（any(error_flag) = .true.），或者开启了调试模式，
               ! 程序会用极其详细的 write 语句将当前虚拟世界发生错误的时间（年月日）、这棵树的生物量、LAI、
               ! 土壤各层水势全部格式化打印出来。随后，直接调用 fatal_error 强行自毁终止整个地球系统模拟，
               ! 以防错误数值像滚雪球一样污染后续的生态数据。
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Find water fluxes.  Note that transp is from last timestep's psi_open  !
               ! and psi_closed.                                                           !
               !---------------------------------------------------------------------------!
               ! D. 水分通量方程组正向求解
               call calc_plant_water_flux(                            &
                        dtlsm                                         &!input
                       ,sap_area,cpatch%nplant(ico),ipft              &!input
                       ,cpatch%is_small(ico),cpatch%krdepth(ico)      &!input
                       ,cpatch%bleaf(ico),bsap,cpatch%broot(ico)      &!input
                       ,cpatch%hite(ico),transp                       &!input
                       ,cpatch%leaf_psi(ico),cpatch%wood_psi(ico)     &!input
                       ,soil_psi,soil_cond,ipa,ico                    &!input
                       ,cpatch%wflux_wl(ico),cpatch%wflux_gw(ico)     &!output
                       ,cpatch%wflux_gw_layer(:,ico))                 !!output
               !! 解析： 真正执行水动力学核心计算的算子。
               ! 检查完全通过后，将准备好的所有参数输入到 calc_plant_water_flux 子程序中。
               ! 该函数在底层解算植物内部的水力结构方程组，并最终计算并更新（输出）三个物理量：
               ! wflux_wl（木质部到叶片的流量）、
               ! wflux_gw（根系总吸水量）、
               ! wflux_gw_layer（根系从每一层土壤中分别吸走了多少水）。
               !---------------------------------------------------------------------------!
            else
               !----- Neither leaves nor wood are resolvable.  Assume zero flow. ----------!
               ! 解析： 如果最开始的 track_hydraulics 判定为假（例如死树或无法解析的群落），则直接假设水分断流，
               ! 所有通量赋为 0。至此，整个群落大循环 cohortloop 以及模式选择 select case 全部正常结束。
               cpatch%wflux_wl(ico) = 0.
               cpatch%wflux_gw(ico) = 0.
               cpatch%wflux_gw_layer(:,ico)  = 0.
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do cohortloop
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the most frequent timescale averages.                                   !
      !------------------------------------------------------------------------------------!
      do ico = 1, cpatch%ncohorts
         ! 解析： 再次启动一个群落循环，用来更新“高频时间尺度平均值（Fast timescale averages）”。
         ! 将本次瞬时计算出的叶片水势、木质部水势、叶片/木质部内部水量（_water_int），
         ! 分别乘以当前时间步长的权重 dtlsm_o_frqsum 并累加。这些数据最终会被用来输出日平均或月平均的生态学报表。
         cpatch%fmean_leaf_psi      (ico) = cpatch%fmean_leaf_psi      (ico)               &
                                          + cpatch%leaf_psi            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_psi      (ico) = cpatch%fmean_wood_psi      (ico)               &
                                          + cpatch%wood_psi            (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_leaf_water_int(ico) = cpatch%fmean_leaf_water_int(ico)               &
                                          + cpatch%leaf_water_int      (ico)               &
                                          * dtlsm_o_frqsum
         cpatch%fmean_wood_water_int(ico) = cpatch%fmean_wood_water_int(ico)               &
                                          + cpatch%wood_water_int      (ico)               &
                                          * dtlsm_o_frqsum
         if (cpatch%dmax_leaf_psi(ico) == 0.) then
             cpatch%dmax_leaf_psi(ico) =  cpatch%leaf_psi(ico)
         else
             cpatch%dmax_leaf_psi(ico) =  max( cpatch%dmax_leaf_psi(ico)                   &
                                             , cpatch%leaf_psi     (ico) )
         end if
         ! 解析： 利用 max 函数，滚动更新每个群落今天经历过的最大（最高）叶片水势 dmax_leaf_psi（通常发生在清晨水分完全恢复时）。
         ! 后续的代码（未贴出部分）对 dmin_leaf_psi（最小水势，常发生于正午暴晒暴蒸腾时）、
         ! 木质部最大/最小水势（dmax_wood_psi / dmin_wood_psi）执行了完全相同的滚动更新操作。
         if (cpatch%dmin_leaf_psi(ico) == 0.) then
             cpatch%dmin_leaf_psi(ico) =  cpatch%leaf_psi(ico)
         else
             cpatch%dmin_leaf_psi(ico) =  min( cpatch%dmin_leaf_psi(ico)                   &
                                             , cpatch%leaf_psi     (ico) )
         end if
         if (cpatch%dmax_wood_psi(ico) == 0.) then
             cpatch%dmax_wood_psi(ico) =  cpatch%wood_psi(ico)
         else
             cpatch%dmax_wood_psi(ico) =  max( cpatch%dmax_wood_psi(ico)                   &
                                             , cpatch%wood_psi     (ico) )
         end if
         if (cpatch%dmin_wood_psi(ico) == 0.) then
             cpatch%dmin_wood_psi(ico) =  cpatch%wood_psi(ico)
         else
             cpatch%dmin_wood_psi(ico) =  min( cpatch%dmin_wood_psi(ico)                   &
                                             , cpatch%wood_psi     (ico) )
         end if
         !---------------------------------------------------------------------------------!
       end do

      return

   end subroutine plant_hydro_driver
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! SUBROUTINE CALC_PLANT_WATER_FLUX
   !> \brief Calculate water flow within plants driven by hydraulic laws
   !> \details This subroutine calculates ground to wood/root, and wood/root to
   !> leaf water flow from initial plant water potential and transpiration rates
   !> based on hydraulic laws. To simplify calculation and reduce numerial error,
   !> It is assumed that within one timestep, leaf water pool << wood water pool
   !> << soil water pool. Special cases are handled when these
   !> assumptions are incorrect such as grasses and very small seedlings.\n
   !>   Note that this subroutine only update fluxes but not the water potential
   !> or water content, which will be updated within the RK4 Integrator.
   !>   Currently, this subroutine works at DTLSM level, but it is
   !> straightforwad to implement the plant hydrodynamics within the RK4
   !> integration scheme. Yet, it is not yet tested how much extra computational
   !> cost would it incur\n
   !> References:\n
   !> [X16] Xu X, Medvigy D, Powers JS, Becknell JM , Guan K. 2016. Diversity in plant 
   !>       hydraulic traits explains seasonal and inter-annual variations of vegetation
   !>       dynamics in seasonally dry tropical forests. New Phytol. 212: 80-95. 
   !>       doi:10.1111/nph.14009.
   !>
   !> [K03] Katul G, Leuning R , Oren R. 2003. Relationship between plant hydraulic and
   !>       biochemical properties derived from a steady-state coupled water and carbon 
   !>       transport model. Plant Cell Environ. 26: 339-350. 
   !>       doi:10.1046/j.1365-3040.2003.00965.x.
   !>
   !> \author Xiangtao Xu, 29 Jan. 2018
   !---------------------------------------------------------------------------------------!
   !! 这个子程序 calc_plant_water_flux 是植物水动力学模型的核心算子。
   !! 它的主要功能是：根据当前时间步长开始时的植物内部水势（叶片、木质部）和土壤水势，
   !! 正向求解“地下$\rightarrow$根系/木质部”以及“木质部$\rightarrow$叶片”的水分流动通量（Flux）。
   subroutine calc_plant_water_flux(dt                                  & !timestep
               ,sap_area,nplant,ipft,is_small,krdepth                   & !plant input
               ,bleaf,bsap,broot,hite                                   & !plant input
               ,transp,leaf_psi,wood_psi                                & !plant input
               ,soil_psi,soil_cond                                      & !soil  input
               ,ipa,ico                                                 & !debug input
               ,wflux_wl,wflux_gw,wflux_gw_layer)                       ! !flux  output
      !!  解析： 定义子程序名和形参列表。输入参数包括时间步长 dt、边材面积、植物密度、植物功能型（ipft）、
      !! 是否为幼树标志（is_small）、最大根深索引（krdepth）、叶/茎/根生物量、
      !! 树高 hite、瞬时蒸腾 transp、叶片水势 leaf_psi、木质部水势 wood_psi 以及各层土壤水势和导水率。
      !!输出为三个水分通量。
      use soil_coms       , only : slz8                 & ! intent(in)
                                 , dslz8                ! ! intent(in)
      use grid_coms       , only : nzg                  ! ! intent(in)
      use consts_coms     , only : pi18                 & ! intent(in)
                                 , lnexp_min8           ! ! intent(in)
      use rk4_coms        , only : tiny_offset          ! ! intent(in)
      ! 解析： 引用外部模块常量。slz8 是土壤层界面深度，dslz8 是土壤层厚度（均为双精度）；nzg 为总土壤层数；pi18 是 $4\pi$ 或相关几何常数；lnexp_min8 是防指数下溢的极小值（防止 exp(x) 的 x 太小导致系统崩溃）；tiny_offset 用于单双精度安全转换。
      use pft_coms        , only : leaf_water_cap       & ! intent(in) 
                                 , wood_water_cap       & ! intent(in)
                                 , leaf_psi_min         & ! intent(in)
                                 , wood_psi_min         & ! intent(in)
                                 , small_psi_min        & ! intent(in)
                                 , wood_psi50           & ! intent(in)
                                 , wood_Kmax            & ! intent(in)
                                 , wood_Kexp            & ! intent(in)
                                 , vessel_curl_factor   & ! intent(in)
                                 , root_beta            & ! intent(in)
                                 , SRA                  & ! intent(in)
                                 , C2B                  ! ! intent(in)
      use ed_misc_coms    , only : current_time         ! ! intent(in)
      implicit none
      ! 解析： 引用植物功能型的生理和水力参数（如饱和导水率 wood_Kmax、脆性曲线半失水水势 wood_psi50 等）。implicit none 强制显式声明变量。
      !----- Arguments --------------------------------------------------------------------!
      real   ,                 intent(in)  :: dt             !time step           [      s]
      real   ,                 intent(in)  :: sap_area       !sapwood_area        [     m2]
      real   ,                 intent(in)  :: nplant         !plant density       [  pl/m2]
      integer,                 intent(in)  :: ipft           !plant funct. type   [    ---]
      integer,                 intent(in)  :: krdepth        !Max. rooting depth  [    ---]
      logical,                 intent(in)  :: is_small       !Small cohort?       [    T|F]
      real   ,                 intent(in)  :: bleaf          !leaf biomass        [    kgC]
      real   ,                 intent(in)  :: bsap           !sapwood biomass     [ kgC/pl]
      real   ,                 intent(in)  :: broot          !fine root biomass   [ kgC/pl]
      real   ,                 intent(in)  :: hite           !plant height        [      m]
      real   ,                 intent(in)  :: transp         !transpiration       [   kg/s]
      real   ,                 intent(in)  :: leaf_psi       !leaf water pot.     [      m]
      real   ,                 intent(in)  :: wood_psi       !wood water pot.     [      m]
      real   , dimension(nzg), intent(in)  :: soil_psi       !soil water pot.     [      m]
      real   , dimension(nzg), intent(in)  :: soil_cond      !soil water cond.    [kg/m2/s]
      integer,                 intent(in)  :: ipa            !Patch index         [    ---]
      integer,                 intent(in)  :: ico            !Cohort index        [    ---]
      real   ,                 intent(out) :: wflux_wl       !wood-leaf flux      [   kg/s]
      real   ,                 intent(out) :: wflux_gw       !ground-wood flux    [   kg/s]
      real   , dimension(nzg), intent(out) :: wflux_gw_layer !wflux_gw for each soil layer
      !----- Temporary double precision variables (input/output). -------------------------!
      ! 解析： 声明形参和局部变量。为了保证积分和指数运算的数值精确度，防止舍入误差导致的数值振荡，该子程序内部所有的物理计算全部采用双精度变量（变量名后带 _d）。
      real(kind=8)                 :: dt_d
      real(kind=8)                 :: sap_area_d
      real(kind=8)                 :: bleaf_d
      real(kind=8)                 :: bsap_d
      real(kind=8)                 :: broot_d
      real(kind=8)                 :: nplant_d
      real(kind=8)                 :: hite_d
      real(kind=8)                 :: transp_d
      real(kind=8)                 :: leaf_psi_d
      real(kind=8)                 :: wood_psi_d
      real(kind=8), dimension(nzg) :: soil_psi_d
      real(kind=8), dimension(nzg) :: soil_cond_d
      real(kind=8)                 :: wflux_wl_d
      real(kind=8)                 :: wflux_gw_d
      real(kind=8), dimension(nzg) :: wflux_gw_layer_d
      !----- Temporary double precision variables (PFT parameters). -----------------------!
      real(kind=8)                 :: leaf_psi_min_d
      real(kind=8)                 :: wood_psi_min_d
      real(kind=8)                 :: leaf_psi_lwr_d
      real(kind=8)                 :: wood_psi_lwr_d
      real(kind=8)                 :: root_beta_d
      real(kind=8)                 :: SRA_d
      real(kind=8)                 :: wood_psi50_d
      real(kind=8)                 :: wood_Kexp_d
      real(kind=8)                 :: wood_Kmax_d
      real(kind=8)                 :: vessel_curl_factor_d
      !----- Auxiliary variables. ---------------------------------------------------------!
      real(kind=8)                          :: exp_term             !exponent term
      real(kind=8)                          :: ap                   ![s-1]
      real(kind=8)                          :: bp                   ![m s-1]
      real(kind=8)                          :: stem_cond            !stem conductance
      real(kind=8)                          :: plc                  !plant loss of conductance
      real(kind=8)                          :: c_leaf               !leaf water capacitance
      real(kind=8)                          :: c_stem               !stem water capacitance
      real(kind=8)                          :: RAI                  !root area index
      real(kind=8)                          :: root_frac            !fraction of roots
      real(kind=8)                          :: proj_leaf_psi        !projected leaf water pot.
      real(kind=8)                          :: proj_wood_psi        !projected wood water pot. 
      real(kind=8)                          :: gw_cond              !g->w water conductivity
      real(kind=8)                          :: org_wood_psi         !used for small tree
      real(kind=8)                          :: org_leaf_psi         !used for small tree
      real(kind=8)                          :: weighted_soil_psi
      real(kind=8)                          :: weighted_gw_cond
      real(kind=8)                          :: above_layer_depth
      real(kind=8)                          :: current_layer_depth
      real(kind=8)                          :: total_water_supply
      real(kind=8)      , dimension(nzg)    :: layer_water_supply
      !----- Counters. --------------------------------------------------------------------!
      integer                               :: k
      !----- Boolean flags. ---------------------------------------------------------------!
      logical                               :: zero_flow_wl
      logical                               :: zero_flow_gw
      logical           , dimension(5)      :: error_flag
      !----- Local constants. -------------------------------------------------------------!
      character(len=13) , parameter         :: efmt       = '(a,1x,es12.5)'
      character(len=9)  , parameter         :: ifmt       = '(a,1x,i5)'
      character(len=9)  , parameter         :: lfmt       = '(a,1x,l1)'
      integer           , parameter         :: dco        = 0
      logical           , parameter         :: debug_flag = .false.
      !----- External function ------------------------------------------------------------!
      real(kind=4)      , external          :: sngloff       ! Safe dble 2 single precision
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Convert all input state vars and some PFT-dependent parameters to double 
      ! precision.
      !------------------------------------------------------------------------------------!
      ! 解析： 将输入的单精度状态变量和参数转换为双精度变量（如 dt_d、sap_area_d）。同时根据植物功能型 ipft 获取根系分布参数 root_beta_d 和比根面积 SRA_d。
      dt_d                 = dble(dt                      )
      sap_area_d           = dble(sap_area                )
      bleaf_d              = dble(bleaf                   )
      bsap_d               = dble(bsap                    )
      broot_d              = dble(broot                   )
      nplant_d             = dble(nplant                   )
      hite_d               = dble(hite                    )
      transp_d             = dble(transp                  )
      leaf_psi_d           = dble(leaf_psi                )
      wood_psi_d           = dble(wood_psi                )
      soil_psi_d           = dble(soil_psi                )
      soil_cond_d          = dble(soil_cond               )
      root_beta_d          = dble(root_beta         (ipft))
      SRA_d                = dble(SRA               (ipft))
      !----- Minimum threshold depends on whether the plant is small or large. ------------!
      ! 解析： 生理阈值判定：如果是小树或草本（is_small），其叶片和木质部的极限生存水势采用 small_psi_min；
      ! 如果是大树，则采用标准的 leaf_psi_min 和 wood_psi_min。
      ! 然后利用缓冲因子计算出允许水势跌落的物理绝对下限 _lwr_d。
      if (is_small) then
         leaf_psi_min_d    = dble(small_psi_min     (ipft))
         wood_psi_min_d    = dble(small_psi_min     (ipft))
      else
         leaf_psi_min_d    = dble(leaf_psi_min      (ipft))
         wood_psi_min_d    = dble(wood_psi_min      (ipft))
      end if 
      leaf_psi_lwr_d       = om_buff_d * leaf_psi_min_d
      wood_psi_lwr_d       = om_buff_d * wood_psi_min_d
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update plant hydrodynamics.
      !
      !     The regular solver assumes stem water pool is way larger than the leaf water
      ! pool. In cases where leaf water pool is of similar magnitude to stem water pool 
      ! (seedlings and grasses), leaf water potential is forced to be the same as stem 
      ! water potential, to maintain numerical stability.  This, however, may bias the 
      ! water potential estimates for these plants.
      !
      ! Water flow is calculated from canopy to roots
      ! Positive flux means upward flow (g->w, w->l, l->air)
      !------------------------------------------------------------------------------------!


      !----- Initialise proj_psi as the starting psi. Also save the initial psi values. ---!
      !! 解析： 将当前的双精度水势赋值给“预测水势（proj_...）”以及“原始水势（org_...）”，作为迭代和积分的起点。
      proj_leaf_psi = leaf_psi_d 
      proj_wood_psi = wood_psi_d
      org_wood_psi  = wood_psi_d
      org_leaf_psi  = leaf_psi_d
      !------------------------------------------------------------------------------------!


      !! 三、 第一部分：计算“木质部 $\rightarrow$ 叶片”水分通量 (Lines 123 - 247)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! 1. Calculate wood/stem/root to leaf water flow
      !! 1. 计算水容与小树/草本特殊处理
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! First, check the relative magnitude of leaf and sapwood water pool
      ! If it is a small tree/grass, force psi_leaf to be the same as psi_stem
      !------------------------------------------------------------------------------------!
      c_leaf = dble(leaf_water_cap(ipft) * C2B) * bleaf_d            ! kg H2O / m
      c_stem = dble(wood_water_cap(ipft) * C2B) * (broot_d + bsap_d) ! kg H2O / m
      ! 解析： 计算群落中植物个体的总叶片水容（c_leaf）和总茎干根系水容（c_stem）。
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If cohort is considered small (see stable_cohorts.f90 for flag), then we do
      ! not distinguish between leaves and wood when calculating hydrodynamics.
      !------------------------------------------------------------------------------------!
      if (is_small) then
         ! 解析： 草本与幼苗的特异性物理简化。
         ! 对于矮小植物，其叶片水容与茎干水容在同一数量级。
         ! 为了防止两片小水力质点之间在短时间步长内发生数值剧烈振荡，模型强制将叶片水势与木质部水势融合成一个均匀的整体水势（质量加权平均）。
         ! 此时把叶片和木质部视作一个单一的水容（c_stem = c_stem + c_leaf），并令木质部向叶片的流量临时等于蒸腾速率。
         !---------------------------------------------------------------------------------!
         !   1.1.  Small tree, force leaf_psi to be the same as wood_psi.  Calculate the 
         !         new veg_psi of mixing leaf and wood
         !---------------------------------------------------------------------------------!
         wood_psi_d   = (c_leaf * leaf_psi_d + c_stem * wood_psi_d) / (c_leaf+c_stem)
         leaf_psi_d   = wood_psi_d
         !---------------------------------------------------------------------------------!


         !----- Use c_leaf+c_stem as the total capacitance to solve stem_psi later. -------!
         c_stem = c_stem + c_leaf
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    In this case, we temporarily assign transpiration as wflux_wl_d since leaves 
         ! and wood are treated as a single entity.  The value will be recalculated once 
         ! we obtain the projected water potential.
         !---------------------------------------------------------------------------------!
         wflux_wl_d   = transp_d
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Set zero_flow_wl to .false., so it appears correctly in the error message in !
         ! case the model crashes.                                                         !
         !---------------------------------------------------------------------------------!
         zero_flow_wl = .false.
         !---------------------------------------------------------------------------------!

      else
         !---------------------------------------------------------------------------------!
         ! 1.2.  Regular case, big trees.
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Special cases in which flow between leaves and wood should be zero.  Perhaps
         ! there are better ways to systematically avoid these traps, but currently, this
         ! is done in a case-by-case manner.
         !
         ! Case 0.  Negative wflux_wl overchange wood storage?  This should never occur
         !          if we initialise leaf water potential slightly lower than the stem
         !          water potential.
         !
         ! Case 1.  Cohort has no leaves.
         !
         ! Case 2.  Both wood and leaves are very dry, and either (a) wood cannot support 
         !          upward sapflow or (b) leaf cannot support downward flow.  This could 
         !          happen for dying trees experiencing extreme drought. 
         !
         ! Case 3.  The cohort just grows out of 'small tree status'.  Their leaves can be 
         !          over-charged with water because gravitational effect was not 
         !          considered for leaf water potential of small trees.  As a result, this 
         !          can lead to a down-ward sapflow, and potentially over-charging the 
         !          sapwood. We need to zero the flow in this case as well, until 
         !          leaf_psi_d drops below wood_psi_d - hite_d.
         !---------------------------------------------------------------------------------!
         zero_flow_wl = ( c_leaf == 0.d0                            ) .or.  & ! Case 1
                        ( leaf_psi_d >= (wood_psi_d - hite_d) .and.         &
                          leaf_psi_d <= leaf_psi_lwr_d              ) .or.  & ! Case 2a
                        ( leaf_psi_d <= (wood_psi_d - hite_d) .and.         &
                          wood_psi_d <= wood_psi_lwr_d              ) .or.  & ! Case 2b
                        ( leaf_psi_d >  (wood_psi_d - hite_d)       )       ! ! Case 3
         ! 解析： 对大树进行断流条件检查（zero_flow_wl）：
         ! Case 1： 没叶子（c_leaf == 0），无法传水。
         ! Case 2a/2b： 极度干旱。当植物已经干枯到极限水势（_lwr_d）以下时，茎干无法再向上输水，或叶片无法逆向输水。
         ! Case 3： 树木刚从小树长成大树，或者叶片水势高于水力静压平衡点（木质部水势减去重力水头 hite_d）。重力势能导致水无法向上流动，必须断流，直到叶片水势因蒸腾掉落到木质部以下。
         !---------------------------------------------------------------------------------




         !---------------------------------------------------------------------------------!
         !    Decide whether or not to calculate sapflow.
         !---------------------------------------------------------------------------------!
         if (zero_flow_wl) then
            !------------------------------------------------------------------------------!
            ! 1.2.1. No need to calculate sapflow
            !------------------------------------------------------------------------------!
            wflux_wl_d = 0.d0

            !------ Proj_leaf_psi is only dependent upon transpiration. -------------------!
            if (c_leaf > 0.) then
                proj_leaf_psi = leaf_psi_d - transp_d * dt_d / c_leaf
            else
                proj_leaf_psi = leaf_psi_d
            end if
            ! 解析： 如果判定断流，则木质部到叶片通量为 0。
            ! 预测时间步结束时的叶片水势 proj_leaf_psi 纯粹由叶片内部残存的水分被蒸腾消耗（$\Delta \psi = \text{transp} \times dt / c\_leaf$）来决定。
            !------------------------------------------------------------------------------!

         else
            !------------------------------------------------------------------------------!
            !     We do need to calculate sapflow.  First convert some PFT-dependent 
            ! parameters to double precision.
            !------------------------------------------------------------------------------!
            wood_psi50_d         = dble(wood_psi50        (ipft))
            wood_Kexp_d          = dble(wood_Kexp         (ipft))
            wood_Kmax_d          = dble(wood_Kmax         (ipft))
            vessel_curl_factor_d = dble(vessel_curl_factor(ipft))
            !------------------------------------------------------------------------------!

            !----- Calculate plant loss of conductivity [dimensionless]. ------------------!
            plc = 1.d0 / (1.d0 + (wood_psi_d / wood_psi50_d) ** wood_Kexp_d)
            ! 解析： 正常传水情况。利用 Sperry 导水率脆性曲线方程 计算植物由于木质部栓塞导致导水率丧失的比例（plc）。
            !------------------------------------------------------------------------------!



            !----- Calculate stem conductance [kg / s]. -----------------------------------!
            stem_cond = wood_Kmax_d * plc                 & ! kg/m/s
                      * sap_area_d                        & ! conducting area m2
                      / (hite_d * vessel_curl_factor_d)   ! ! conducting length m
            ! 解析： 根据达西定律计算整株树木木质部导水度 stem_cond（单位 $kg/s$）。
            ! $$\text{stem\_cond} = \frac{\text{最大导水率} \times \text{栓塞残余系数} \times \text{边材面积}}{\text{树高} \times \text{导管弯曲因子}}$$
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Find sapflow.
            !------------------------------------------------------------------------------!
            if (stem_cond == 0.) then
               !---- 1.2.2. Zero flux because stem conductivity is also zero. -------------!
               wflux_wl_d = 0.d0
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               ! 1.2.3. "Normal case", with positive c_leaf and positive stem_cond.  Check
               !        reference X16 for derivation of the equations.
               !---------------------------------------------------------------------------!
               ap = - stem_cond / c_leaf                                            ! [1/s]
               bp = ((wood_psi_d - hite_d) * stem_cond - transp_d) / c_leaf         ! [m/s]

               !----- Project the final leaf psi. -----------------------------------------!
               exp_term      = exp(max(ap * dt_d,lnexp_min8))
               proj_leaf_psi = max( leaf_psi_lwr_d                                         &
                                  , ((ap * leaf_psi_d + bp) * exp_term - bp) / ap )
               !---------------------------------------------------------------------------!


               !----- Calculate the average sapflow rate within the time step [kgH2O/s]. --!
               wflux_wl_d = (proj_leaf_psi - leaf_psi_d) * c_leaf / dt_d + transp_d
               !---------------------------------------------------------------------------!
               ! 解析： 一阶线性常微分方程的解析解求解（核心）。
               ! 叶片水分平衡方程为：$c\_leaf \frac{d\psi_l}{dt} = \text{stem\_cond}(\psi_w - \psi_l - h) - \text{transp}$。
               ! 代码将其化简为标准微分形式 $\frac{d\psi_l}{dt} = ap \cdot \psi_l + bp$。
               ! 利用解析解公式 $\psi_l(t) = \frac{(ap \cdot \psi_{l0} + bp)e^{ap \cdot dt} - bp}{ap}$ 直接精准计算出时间步结束时的叶片预测水势 proj_leaf_psi，并带有最干水势截断保护。最后反推求出此时间步内的平均流速 wflux_wl_d。
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! 2.  Calculate ground -> wood/stem/root water flow
      ! 四、 第二部分：计算“土壤 $\rightarrow$ 根系/木质部”水分通量 (Lines 248 - 378)
      !------------------------------------------------------------------------------------!
      weighted_soil_psi  = 0.d0
      weighted_gw_cond   = 0.d0
      layer_water_supply = 0.d0
      total_water_supply = 0.d0

      !----- Loop over all soil layers to get the aggregated water conductance. -----------!
      do k = krdepth,nzg
         ! 解析： 遍历从最大根深层 krdepth 到最底层的各土壤层。
         ! 使用 Gale 和 Grigal 模型（基于 root_beta 的指数分布）计算当前垂直土壤层 k 内所分配到的细根比例 root_frac。
         !---------------------------------------------------------------------------------!
         !   Define layer edges 
         !
         !---------------------------------------------------------------------------------!
         current_layer_depth = -slz8(k)
         above_layer_depth   = -slz8(k+1)
         !---------------------------------------------------------------------------------!



         !----- Calculate the root fraction of this layer. --------------------------------!
         root_frac = ( root_beta_d ** (above_layer_depth   / (-slz8(krdepth)))             &
                     - root_beta_d ** (current_layer_depth / (-slz8(krdepth))) )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !  Calculate RAI in each layer.                                                   !
         !---------------------------------------------------------------------------------!
         RAI = broot_d * SRA_d * root_frac * nplant_d  ! m2/m2
         ! 解析： 计算当前层内的根面积指数（RAI）
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !    Calculate soil-root water conductance kg H2O/m/s based on reference [K03].
         !---------------------------------------------------------------------------------!
         gw_cond = soil_cond_d(k) * sqrt(RAI) / (pi18 * dslz8(k))  & ! kg H2O / m3 / s
                 / nplant_d                                        ! ! conducting area  m2
         ! 解析： 依据 Katul 等人 (2003) 的经典公式，综合当前层土壤自身的导水率 soil_cond_d 和根系密度（$\sqrt{\text{RAI}}$），
         ! 解析出土壤到根系间的导水度 gw_cond。
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Disable hydraulic redistribution.  Assume roots will shut down if they are 
         ! going to lose water to soil.
         !---------------------------------------------------------------------------------!
         if (soil_psi_d(k) <= wood_psi_d) then
            gw_cond = 0.d0
         end if
         ! 解析： 关闭水力提升/水力重分配功能。 
         ! 如果当前层的土壤水势比植物木质部还低（土壤比树干干），植物根系会自动关闭该层的通道，防止植物体内的水分逆向倒流回土壤。
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Calculate weighted conductance, weighted psi, and water_supply_layer_frac.
         !---------------------------------------------------------------------------------!
         weighted_gw_cond      = weighted_gw_cond + gw_cond                  ! kgH2O/m/s
         weighted_soil_psi     = weighted_soil_psi + gw_cond * soil_psi_d(k) ! kgH2O/s
         layer_water_supply(k) = gw_cond * (soil_psi_d(k) - wood_psi_d)      ! kgH2O/s
         ! 解析： 在整个垂直土壤剖面上进行水力集成累加，得到群落的总土壤传导度 weighted_gw_cond、由传导度加权的土壤总水势、以及每一层理论上能为植物提供的水分贡献量 layer_water_supply(k)。
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Now we can calculate ground->wood water flow.
      ! First we handle special cases
      ! 2. 求解根系吸水总通量
      !------------------------------------------------------------------------------------!
      zero_flow_gw = (c_stem           == 0.d0) .or.  & ! No sapwood or fine  roots
                     (weighted_gw_cond == 0.d0)       ! ! soil is drier than wood

      if (zero_flow_gw) then
      ! 解析： 如果植物没有生物量或者所有土壤层都极干（weighted_gw_cond == 0），则判定地下无法吸水（zero_flow_gw = .true.）。此时地下总总通量 wflux_gw_d = 0。木质部预测水势 proj_wood_psi 纯粹由刚才算出的木质部向叶片的输水消耗来决定。
         !---------------------------------------------------------------------------------!
         !     No need to calculate water flow: wood psi is only dependent upon sapflow.
         !---------------------------------------------------------------------------------!
         wflux_gw_d    = 0.d0
         if (c_stem > 0.) then
            !----- Make sure that projected wood psi will be bounded. ---------------------!
            wflux_wl_d    = min(wflux_wl_d, (wood_psi_d - wood_psi_lwr_d) * c_stem / dt_d )
            proj_wood_psi = wood_psi_d - wflux_wl_d * dt_d / c_stem
            !------------------------------------------------------------------------------!
         else
            proj_wood_psi = wood_psi_d
         end if
         !---------------------------------------------------------------------------------!
      else
         ! 解析： 正常吸水情况。采用与前面叶片完全相同的解析解机制，求解木质部的水分平衡微分方程：$c\_stem \frac{d\psi_w}{dt} = \sum [\text{gw\_cond} \cdot (\psi_{sk} - \psi_w)] - \text{wflux\_wl}$。通过解析解得出预测的木质部终点水势 proj_wood_psi，并反推算出地下总吸水通量 wflux_gw_d。
         !---------------------------------------------------------------------------------!
         !     Calculate the average soil water uptake. Check reference X16 for derivation
         ! of the equations.
         !---------------------------------------------------------------------------------!
         ap = - weighted_gw_cond  / c_stem  ! ! 1/s
         bp = (weighted_soil_psi - wflux_wl_d) / c_stem ! m/s
         !---------------------------------------------------------------------------------!

         !----- Project the final wood psi, but ensure it will be bounded. ----------------!
         exp_term        = exp(max(ap * dt_d,lnexp_min8))
         proj_wood_psi   = max( wood_psi_lwr_d                                             &
                              , ((ap * wood_psi_d + bp) * exp_term - bp) / ap )
         !---------------------------------------------------------------------------------!


         !----- Calculate the average root extraction within the time step [kgH2O/s]. -----!
         wflux_gw_d     = (proj_wood_psi - wood_psi_d) * c_stem  / dt_d + wflux_wl_d
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Re-calculate the water fluxes in the cases of small cohorts.
      !------------------------------------------------------------------------------------!
      if (is_small) then
         ! 解析： 针对小树/草本的最终通量修正。 
         ! 由于在最前面将小树的叶片和木质部强制融合了，它们的最终预测水势必须保持一致。
         !此处令小树的叶片水势等于刚刚解出的木质部水势，并重新计算出小树体内真实的木质部$\rightarrow$叶片水分流动通量 wflux_wl_d。
         !---------------------------------------------------------------------------------!
         !     Ground->wood flux (wflux_gw_d) is correct, no need to update.  However, we
         ! do need to update wood->leaf (wflux_wl_d).
         !---------------------------------------------------------------------------------!
         proj_leaf_psi = proj_wood_psi
         wflux_wl_d    = (proj_leaf_psi - org_leaf_psi)  * c_leaf / dt_d + transp_d
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Now estimate the water uptake from each layer based on layer_water_supply.
      !------------------------------------------------------------------------------------!
      if (sum(layer_water_supply) == 0.d0) then
         wflux_gw_layer_d = 0.d0
      else
         wflux_gw_layer_d = layer_water_supply / sum(layer_water_supply) * wflux_gw_d
      end if
      ! 解析： 按比例分配吸水通量。 将总吸水量 wflux_gw_d 按照各土壤层刚才算出的供水能力比例（layer_water_supply / 总供水），
      ! 精准拆分到每一层土壤上（wflux_gw_layer_d）。这决定了下一步模型将从哪层土壤里扣除多少水。
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!



      !! 六、 物理安全性检查 (Lines 411 - 511)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! 3.  Sanity check.  Stop the simulation in case anything went wrong.   Sympotms for
      !     things going wrong:
      !     a.  NaN values --- Run the debugger
      !     b.  Projected leaf/wood potential is positive
      !     c.  Current leaf/wood potential is positive
      !     d.  Projected leaf/wood potential is less than minimum acceptable
      !     e.  Current leaf/wood potential is less than minimum acceptable
      !------------------------------------------------------------------------------------!
      error_flag(1) = isnan(wflux_wl_d)              .or. isnan(wflux_gw_d)
      error_flag(2) = proj_leaf_psi > 0.             .or. proj_wood_psi > 0.
      error_flag(3) = leaf_psi_d    > 0.             .or. wood_psi_d    > 0.
      error_flag(4) = proj_leaf_psi < leaf_psi_min_d .or. proj_wood_psi < wood_psi_min_d
      error_flag(5) = leaf_psi_d    < leaf_psi_min_d .or. wood_psi_d    < wood_psi_min_d
      ! 解析： 构建全方位的崩溃防御墙。 检查以下任何一项水动力学物理状态是否失真：产生 NaN、水势变成正数、或者水势低于植物毁灭性凋萎极限。

      if ( (debug_flag .and. (dco == 0 .or. ico == dco)) .or. any(error_flag)) then
         ! 解析： 只要触发了上述任意一个严重错误，系统会启动最大化输出，把当前所有的形态、环境、瞬时参数全部报错打印，然后调用 fatal_error 强行中止全球模拟，防止错误的物理量在时间步迭代中向后续系统扩散。
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(a)'  ) ' Plant hydrodynamics inconsistency detected!!'
         write (unit=*,fmt='(92a)') ('-',k=1,92)
         write (unit=*,fmt='(a,i4.4,2(1x,i2.2),1x,f6.0)') ' TIME           : '             &
                                                     ,current_time%year,current_time%month &
                                                     ,current_time%date,current_time%time
         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=ifmt   ) ' + IPA              =',ipa
         write (unit=*,fmt=ifmt   ) ' + ICO              =',ico
         write (unit=*,fmt=ifmt   ) ' + PFT              =',ipft
         write (unit=*,fmt=ifmt   ) ' + KRDEPTH          =',krdepth
         write (unit=*,fmt=efmt   ) ' + HEIGHT           =',hite

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=lfmt   ) ' + IS_SMALL         =',is_small
         write (unit=*,fmt=lfmt   ) ' + ZERO_FLOW_WL     =',zero_flow_wl
         write (unit=*,fmt=lfmt   ) ' + ZERO_FLOW_GW     =',zero_flow_gw

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + BLEAF            =',bleaf
         write (unit=*,fmt=efmt   ) ' + BSAPWOOD         =',bsap
         write (unit=*,fmt=efmt   ) ' + BROOT            =',broot
         write (unit=*,fmt=efmt   ) ' + SAPWOOD_AREA     =',sap_area

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=lfmt   ) ' + Finite fluxes     =',.not. error_flag(1)
         write (unit=*,fmt=lfmt   ) ' + Negative Proj Psi =',.not. error_flag(2)
         write (unit=*,fmt=lfmt   ) ' + Negative Curr Psi =',.not. error_flag(3)
         write (unit=*,fmt=lfmt   ) ' + Bounded Proj Psi  =',.not. error_flag(4)
         write (unit=*,fmt=lfmt   ) ' + Bounded Curr Psi  =',.not. error_flag(5)

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI_MIN      =',leaf_psi_min (ipft)
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI_MIN      =',wood_psi_min (ipft)
         write (unit=*,fmt=efmt   ) ' + SMALL_PSI_MIN     =',small_psi_min(ipft)

         write (unit=*,fmt='(a)'  ) ' '
         write (unit=*,fmt=efmt   ) ' + TRANSP           =',transp
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI (INPUT) =',leaf_psi
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI (INPUT) =',wood_psi
         write (unit=*,fmt=efmt   ) ' + LEAF_PSI (PROJ.) =',proj_leaf_psi
         write (unit=*,fmt=efmt   ) ' + WOOD_PSI (PROJ.) =',proj_wood_psi
         write (unit=*,fmt=efmt   ) ' + WFLUX_GW         =',wflux_gw_d
         write (unit=*,fmt=efmt   ) ' + WFLUX_WL         =',wflux_wl_d


         write (unit=*,fmt='(a)'        ) ' '
         write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
         write (unit=*,fmt='(a,2(1x,a))') '    K','    SOIL_PSI','WFLUX_GW_LYR'
         write (unit=*,fmt='(92a)'      ) ('-',k=1,92)
         do k = 1, nzg
            write (unit=*,fmt='(i5,2(1x,es12.5))') k,soil_psi(k),wflux_gw_layer_d(k)
         end do
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(92a)') ('=',k=1,92)
         write (unit=*,fmt='(a)'  ) ' '

         if (any(error_flag)) then 
            call fatal_error('Plant Hydrodynamics is off-track.'                           &
                            ,'calc_plant_water_flux','plant_hydro.f90')
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
 

      ! 七、 精度还原与变量输出
      !------------------------------------------------------------------------------------!
      !     Copy all the results to output variables.
      !------------------------------------------------------------------------------------!
      wflux_wl = sngloff(wflux_wl_d,tiny_offset)
      do k = 1, nzg
         wflux_gw_layer(k) = sngloff(wflux_gw_layer_d(k),tiny_offset)
      end do
      wflux_gw = sum(wflux_gw_layer)
      ! 解析： 将数据交还主程序。 使用安全单精度转换函数 sngloff（结合 tiny_offset 偏移量防止截断变为零值），将内部高精度的双精度通量结果还原为单精度，并赋值给外部的输出形参 wflux_wl 和 wflux_gw_layer。地下总吸水量 wflux_gw 则是各层单精度吸水量的总和。
      !------------------------------------------------------------------------------------!


      return
   end subroutine calc_plant_water_flux
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !   Util functions/subroutines for plant hydrodynamic calculations                      !
   !=======================================================================================!
   !=======================================================================================!
   ! 这段代码是陆地生态系统陆面过程模型中，专门为了处理植物水动力学（Plant Hydrodynamics）而编写的一组工具函数/子程序（Utility Functions）。






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2RWC           
   !> \breif Convert water potential of leaf and wood to relative water content
   !> \details Here we assume a constant hydraulic capacitance for both leaf and
   !> wood. From the definition of hydraulic capacitance we have \n
   !>       hydro_cap = delta_water_content / delta_psi \n
   !> Since psi = 0. when water_content is at saturation, we have \n
   !>       hydro_cap = (1. - rwc) * water_content_at_saturation / (0. - psi) \n
   !> Reorganize the equation above, we can get \n
   !>       rwc = 1. + psi * hydro_cap / water_content_at_saturation
   !=======================================================================================!
   subroutine psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
   !! 水势 $\rightarrow$ 相对含水量
   !! 输入叶片和木质部的水势（psi），计算出它们各自的相对含水量（rwc）。
   !! 物理原理： 假设植物组织的水容（Hydraulic Capacitance）是恒定的。
   !! 由于水势降得越低，植物组织脱水越严重，因此相对含水量 rwc 就会越低。
   !! 公式为：$$RWC = 1 + \frac{\psi \times \text{水容}}{\text{饱和含水量}}$$
      use pft_coms          ,   only : leaf_water_cap       & ! intent(in)
                                     , wood_water_cap       & ! intent(in)
                                     , leaf_water_sat       & ! intent(in)
                                     , wood_water_sat       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi    ! Water potential of leaves            [  m]
      real      , intent(in)    ::  wood_psi    ! Water potential of wood              [  m]
      integer   , intent(in)    ::  ipft        ! plant functional type                [  -]
      real      , intent(out)   ::  leaf_rwc    ! Relative water content of leaves     [0-1]
      real      , intent(out)   ::  wood_rwc    ! Relative water content of wood       [0-1]

      ! first caculate for leaf
      leaf_rwc  =   1.  + leaf_psi                  & ! [m]
                        * leaf_water_cap(ipft)      & ! [kg H2O/kg biomass/m]
                        / leaf_water_sat(ipft)      ! ! [kg H2O/kg biomass]

      ! same for wood
      wood_rwc  =   1.  + wood_psi                  & ! [m]
                        * wood_water_cap(ipft)      & ! [kg H2O/kg biomass/m]
                        / wood_water_sat(ipft)      ! ! [kg H2O/kg biomass]

      return
   end subroutine psi2rwc
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: RWC2PSI
   !> \brief Convert relative water content to water potential
   !> \details The inverse of psi2rwc
   !=======================================================================================!
   subroutine rwc2psi(leaf_rwc,wood_rwc,ipft,leaf_psi,wood_psi)
      !! rwc2psi —— 相对含水量 $\rightarrow$ 水势
      !! 核心功能： psi2rwc 的逆运算。输入叶片和木质部的相对含水量，反推它们的水势。
      !! 应用场景： 当积分器更新了植物体内的含水量后，模型需要用这个函数把水量转回“水势压力”，以此来驱动下一个时间步的水分流动。
      use pft_coms          ,   only : leaf_water_cap       & ! intent(in)
                                     , wood_water_cap       & ! intent(in)
                                     , leaf_water_sat       & ! intent(in)
                                     , wood_water_sat       ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_rwc    ! Relative water content of leaves     [0-1]
      real      , intent(in)    ::  wood_rwc    ! Relative water content of wood       [0-1]
      integer   , intent(in)    ::  ipft        ! plant functional type                [  -]
      real      , intent(out)   ::  leaf_psi    ! Water potential of leaves            [  m]
      real      , intent(out)   ::  wood_psi    ! Water potential of wood              [  m]

      ! first caculate for leaf
      leaf_psi  =   (leaf_rwc - 1.) * leaf_water_sat(ipft) / leaf_water_cap(ipft)
      ! same for wood
      wood_psi  =   (wood_rwc - 1.) * wood_water_sat(ipft) / wood_water_cap(ipft)

      return
   end subroutine rwc2psi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: RWC2TW            
   !> \brief Convert relative water content to total water for both leaf and wood
   !> \details total water content = \n
   !>  relative water content * water_content_at_saturation * biomass
   !> \warning In the hydro version, bsapwood is set as 0 and bdead is assumedto contain
   !> both sapwood and heart wood. Root is counted as wood.  When dynamic hydraulics is
   !> turned off, we account for water in sapwood using the sapwood biomass definition.
   !=======================================================================================!
   subroutine rwc2tw(leaf_rwc,wood_rwc,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh   &
                    ,ipft,leaf_water_int,wood_water_int)
   !! rwc2tw —— 相对含水量 $\rightarrow$ 绝对含水量
   !! 核心功能： 输入含水量比例（rwc），结合树木的生物量，算出单株树木体内到底存了多少公斤（$kg/plant$）的液体水。
   !! 重要细节（Case分类）：
   !! 传统模式（Case 0）： 粗暴地认为水只存在于活着的细根和边材（Sapwood）生物量里。
   !!动态水动力学模式（Default）： 更加科学。它先通过胸径（dbh）调用异速生长方程算出心材和边材的比例（sap_frac），然后认为所有的活组织（细根、边材、心材中的导管）共同组成了单株植物的木质部储水库。
      use pft_coms       , only : leaf_water_sat      & ! intent(in)
                                , wood_water_sat      & ! intent(in)
                                , C2B                 ! ! intent(in)
      use allometry      , only : dbh2sf              ! ! function
      use physiology_coms, only : plant_hydro_scheme  ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_rwc       ! Relative water content of leaves  [0-1]
      real      , intent(in)    ::  wood_rwc       ! Relative water content of wood    [0-1]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(out)   ::  wood_water_int ! Total internal water of wood      [ kg]
      !----- Local variables. -------------------------------------------------------------!
      real                      ::  sap_frac    ! Fraction of sapwood to basal area    [0-1]
      !------------------------------------------------------------------------------------!


      !----- Leaf.  This is the same, regardless of the plant hydraulic scheme. -----------!
      leaf_water_int    =   leaf_rwc * leaf_water_sat(ipft) * bleaf * C2B 
      !------------------------------------------------------------------------------------!


      !----- Wood.  Check the scheme to decide between sapwood fraction or biomass. -------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Use sapwood biomass to obtain sapwood internal water. ---------------------!
         wood_water_int = wood_rwc * wood_water_sat(ipft) * C2B                            &
                        * (broot + bsapwooda + bsapwoodb )
         !---------------------------------------------------------------------------------!
      case default
         !----- Find the sapwood fraction. ------------------------------------------------!
         sap_frac       = dbh2sf(dbh,ipft)
         !---------------------------------------------------------------------------------!

         !----- Total water only includes live biomass (fine roots and sapwood). ----------!
         wood_water_int = wood_rwc * wood_water_sat(ipft) * C2B                            &
                        * ( broot + (bdeada + bdeadb + bsapwooda + bsapwoodb) * sap_frac )
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      return
   end subroutine rwc2tw
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TW2RWC            
   !> \brief Convert total water to relative water content for both leaf and wood
   !> \details the inverse of rwc2tw \n
   !=======================================================================================!
   subroutine tw2rwc(leaf_water_int,wood_water_int,is_small,bleaf,bsapwooda,bsapwoodb      &
                    ,bdeada,bdeadb,broot,dbh,ipft,leaf_rwc,wood_rwc)
   !! tw2rwc —— 绝对含水量 $\rightarrow$ 相对含水量
   !! 核心功能： rwc2tw 的逆运算。输入单株植物的实际总水量（$kg$），
   !! 除以该株植物达到完全饱和时的最大持水量（tot_water_sat），重新得到 0.0 ~ 1.0 之间的相对含水量比例。
   !! 防御机制： 代码中考虑了极干或没有生物量的极端情况（分母为 0），
   !! 此时会强制给一个安全的最低含水量下限（small_rwc_min / leaf_rwc_min），防止程序因除以 0 而崩溃。
      use pft_coms       , only : leaf_water_sat     & ! intent(in)
                                , wood_water_sat     & ! intent(in)
                                , leaf_rwc_min       & ! intent(in)
                                , wood_rwc_min       & ! intent(in)
                                , small_rwc_min      & ! intent(in)
                                , C2B                ! ! intent(in)
      use allometry      , only : dbh2sf             ! ! function
      use physiology_coms, only : plant_hydro_scheme ! ! intent(in)
      use consts_coms    , only : tiny_num           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      logical   , intent(in)    ::  is_small       ! Small/large plant flag            [T|F]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_rwc       ! Relative water content of leaves  [0-1]
      real      , intent(out)   ::  wood_rwc       ! Relative water content of wood    [0-1]
      !----- Local variables --------------------------------------------------------------!
      real                      :: tot_water_sat
      real                      :: sap_frac        ! Fraction of sapwood to basal area [0-1]
      !------------------------------------------------------------------------------------!


      !----- Leaf.  This is the same, regardless of the plant hydraulic scheme. -----------!
      tot_water_sat = leaf_water_sat(ipft) * C2B * bleaf
      if (tot_water_sat > tiny_num) then
         leaf_rwc  = leaf_water_int / tot_water_sat
      elseif (is_small) then
         leaf_rwc  = op_buff * small_rwc_min(ipft)
      else
         leaf_rwc  = op_buff * leaf_rwc_min(ipft)
      end if
      !------------------------------------------------------------------------------------!


      !----- Wood.  Check the scheme to decide between sapwood fraction or biomass. -------!
      select case (plant_hydro_scheme)
      case (0)
         !----- Use sapwood biomass to obtain sapwood internal water. ---------------------!
         tot_water_sat = wood_water_sat(ipft) * C2B                                        &
                       * (broot + bsapwooda + bsapwoodb)
         !---------------------------------------------------------------------------------!
      case default
         !----- Find the sapwood fraction. ------------------------------------------------!
         sap_frac = dbh2sf(dbh,ipft)
         !---------------------------------------------------------------------------------!


         !----- Find the sapwood fraction. ------------------------------------------------!
         tot_water_sat = wood_water_sat(ipft) * C2B                                        &
                       * (broot + (bdeada + bdeadb + bsapwooda + bsapwoodb) * sap_frac)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------ Make sure the denominator is not zero. --------------------------------------!
      if (tot_water_sat > tiny_num) then
         wood_rwc = wood_water_int / tot_water_sat
      elseif (is_small) then
         wood_rwc = op_buff * small_rwc_min(ipft)
      elseif (is_small) then
         wood_rwc = op_buff * wood_rwc_min(ipft)
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine tw2rwc
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2TW            
   !> \brief Convert water potential to total water for both leaf and wood
   !=======================================================================================!
   subroutine psi2tw(leaf_psi,wood_psi,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh   &
                    ,ipft,leaf_water_int,wood_water_int)
      !! psi2tw：直接把水势（压力）变成总水量（公斤）。内部逻辑是先调用 psi2rwc，紧接着把结果传给 rwc2tw。
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi       ! Water potential of leaves         [  m]
      real      , intent(in)    ::  wood_psi       ! Water potential of wood           [  m]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(out)   ::  wood_water_int ! Total internal water of wood      [ kg]
      !----- Local Variables --------------------------------------------------------------!
      real                      ::  leaf_rwc       ! Relative water content of leaf    [  -]
      real                      ::  wood_rwc       ! Relative water content of wood    [  -]
      !------------------------------------------------------------------------------------!

      ! first convert to rwc
      call psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
      ! second convert to tw
      call rwc2tw(leaf_rwc,wood_rwc,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh,ipft &
                 ,leaf_water_int,wood_water_int)

      return

   end subroutine psi2tw
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TW2PSI            
   !> \brief Convert total water to water potential for both leaf and wood
   !> \details the inverse of psi2tw \n
   !=======================================================================================!
   subroutine tw2psi(leaf_water_int,wood_water_int,is_small,bleaf,bsapwooda,bsapwoodb      &
                    ,bdeada,bdeadb,broot,dbh,ipft,leaf_psi,wood_psi)
      !! tw2psi：直接把总水量（公斤）变成水势（压力）。内部逻辑是先调用 tw2rwc，再调用 rwc2psi。
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf      [ kg]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood      [ kg]
      logical   , intent(in)    ::  is_small       ! Small/large plant flag            [T|F]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf                   [kgC]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass       [kgC]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass       [kgC]
      real      , intent(in)    ::  bdeada         ! Aboveground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  bdeadb         ! Belowground (heart)wood biomass   [kgC]
      real      , intent(in)    ::  broot          ! Biomass of fine root              [kgC]
      real      , intent(in)    ::  dbh            ! Diameter at breast height         [ cm]
      integer   , intent(in)    ::  ipft           ! Plant functional type             [  -]
      real      , intent(out)   ::  leaf_psi       ! Water potential of leaves         [  m]
      real      , intent(out)   ::  wood_psi       ! Water potential of wood           [  m]
      !----- Local Variables --------------------------------------------------------------!
      real                      ::  leaf_rwc       ! Relative water content of leaf    [  -]
      real                      ::  wood_rwc       ! Relative water content of wood    [  -]
      !------------------------------------------------------------------------------------!

      ! first convert to rwc
      call tw2rwc(leaf_water_int,wood_water_int,is_small,bleaf,bsapwooda,bsapwoodb,bdeada  &
                 ,bdeadb,broot,dbh,ipft,leaf_rwc,wood_rwc)
      ! second convert to psi
      call rwc2psi(leaf_rwc,wood_rwc,ipft,leaf_psi,wood_psi)

      return

   end subroutine tw2psi
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE twi2twe
   !> \brief  Intensive to extensive internal water converter.
   !> \details This subroutine converts intensive internal water (kg/plant) to extensive
   !>          water content (kg/m2).  To avoid energy leaks, we assume that all water
   !>          stored in sapwood is included in the heat capacity.  This is not the most
   !>          elegant solution, and in the future, we should make it only the fraction
   !>          associated with branches, but this requires additional changes in the
   !>          budget checks.
   !> \author Marcos Longo 08 Sep 2019
   !---------------------------------------------------------------------------------------!
   subroutine twi2twe(leaf_water_int,wood_water_int,nplant,leaf_water_im2,wood_water_im2)
   !!  尺度转换（单株 $\leftrightarrow$ 每平方米）
   !! Intensive（单株尺度）： 单位是 $kg / plant$（每棵树多少公斤水）。
   !! Extensive（每平方米尺度）： 单位是 $kg / m^2$（每平方米地面上有多少公斤水）。

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_water_int ! Total internal water of leaf  [ kg/pl]
      real      , intent(in)    ::  wood_water_int ! Total internal water of wood  [ kg/pl]
      real      , intent(in)    ::  nplant         ! Stem density                  [ pl/m2]
      real      , intent(out)   ::  leaf_water_im2 ! Extensive leaf internal water [ kg/m2]
      real      , intent(out)   ::  wood_water_im2 ! Water potential of wood       [ kg/m2]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Convert water from intensive to extensive.                                    !
      !------------------------------------------------------------------------------------!
      leaf_water_im2 = nplant * leaf_water_int
      wood_water_im2 = nplant * wood_water_int
      !------------------------------------------------------------------------------------!

      return
   end subroutine twi2twe
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE twe2twi
   !> \brief  Extensive to extensive internal water converter.
   !> \details This subroutine converts extensive internal water (kg/m2) to extensive 
   !>          water content (kg/plant).  To avoid energy leaks, we assume that all water
   !>          stored in wood is included in the heat capacity.  In the future, we should
   !>          make it only the fraction associated with branches, but this requires 
   !>          additional changes in the budget checks.
   !> \author Marcos Longo 08 Sep 2019
   !---------------------------------------------------------------------------------------!
   subroutine twe2twi(leaf_water_im2,wood_water_im2,nplant,leaf_water_int,wood_water_int)
   !!  尺度转换（单株 $\leftrightarrow$ 每平方米）
   !! Intensive（单株尺度）： 单位是 $kg / plant$（每棵树多少公斤水）。
   !! Extensive（每平方米尺度）： 单位是 $kg / m^2$（每平方米地面上有多少公斤水）。
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)  ::  leaf_water_im2 ! Extensive leaf internal water   [ kg/m2]
      real      , intent(in)  ::  wood_water_im2 ! Water potential of wood         [ kg/m2]
      real      , intent(in)  ::  nplant         ! Stem density                    [ pl/m2]
      real      , intent(out) ::  leaf_water_int ! Total internal water of leaf    [ kg/pl]
      real      , intent(out) ::  wood_water_int ! Total internal water of wood    [ kg/pl]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Convert water from intensive to extensive.                                    !
      !------------------------------------------------------------------------------------!
      leaf_water_int = leaf_water_im2 / nplant
      wood_water_int = wood_water_im2 / nplant
      !------------------------------------------------------------------------------------!

      return
   end subroutine twe2twi
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: PSI2TWE
   !> \breif Convert water potential of leaf and wood to extensive water storage
   !> \details This sub-routine is useful when we need to go from water potential to 
   !>          storage, but we don't need the intermediate quantities (relative and
   !>          individual water contents).
   !=======================================================================================!
   subroutine psi2twe(leaf_psi,wood_psi,ipft,nplant,bleaf,bsapwooda,bsapwoodb,bdeada       &
                     ,bdeadb,broot,dbh,leaf_water_im2,wood_water_im2)
      !! 水势 $\rightarrow$ 斑块尺度总水量
      ! 核心功能： 这是一个三合一的高级封装函数。当外部主循环只需要知道“当前水势在每平方米的土地上对应多少储水量”，
      ! 而不需要中间的单株和比例变量时，直接调用它。
      ! 内部流水线： 它在内部一气呵成地执行了：psi2rwc（水势变比例） $\rightarrow$ rwc2tw（比例变单株公斤数） $\rightarrow$ twi2twe（单株变每平方米公斤数）。
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      real      , intent(in)    ::  leaf_psi       ! Water potential of leaves     [     m]
      real      , intent(in)    ::  wood_psi       ! Water potential of wood       [     m]
      integer   , intent(in)    ::  ipft           ! Plant functional type         [     -]
      real      , intent(in)    ::  nplant         ! Stem density                  [ pl/m2]
      real      , intent(in)    ::  bleaf          ! Biomass of leaf               [kgC/pl]
      real      , intent(in)    ::  bsapwooda      ! Aboveground sapwood biomass   [kgC/pl]
      real      , intent(in)    ::  bsapwoodb      ! Belowground sapwood biomass   [kgC/pl]
      real      , intent(in)    ::  bdeada         ! Aboveground heartwood biomass [kgC/pl]
      real      , intent(in)    ::  bdeadb         ! Belowground heartwood biomass [kgC/pl]
      real      , intent(in)    ::  broot          ! Biomass of fine root          [kgC/pl]
      real      , intent(in)    ::  dbh            ! Diameter at breast height     [    cm]
      real      , intent(out)   ::  leaf_water_im2 ! Extensive leaf internal water [ kg/m2]
      real      , intent(out)   ::  wood_water_im2 ! Extensive wood internal water [ kg/m2]
      !----- Local variables. -------------------------------------------------------------!
      real                      ::  leaf_rwc       ! Relative leaf water content   [    --]
      real                      ::  wood_rwc       ! Relative wood water content   [    --]
      real                      ::  leaf_water_int ! Intensive leaf internal water [ kg/pl]
      real                      ::  wood_water_int ! Intensive wood internal water [ kg/pl]
      !------------------------------------------------------------------------------------!

      !----- 1. Potential -> relative water content. --------------------------------------!
      call psi2rwc(leaf_psi,wood_psi,ipft,leaf_rwc,wood_rwc)
      !----- 2. Relative water content -> Intensive internal water. -----------------------!
      call rwc2tw(leaf_rwc,wood_rwc,bleaf,bsapwooda,bsapwoodb,bdeada,bdeadb,broot,dbh,ipft &
                 ,leaf_water_int,wood_water_int)
      !----- 3. Intensive internal water -> Extensive internal water. ---------------------!
      call twi2twe(leaf_water_int,wood_water_int,nplant,leaf_water_im2,wood_water_im2)
      !------------------------------------------------------------------------------------!


      return
   end subroutine psi2twe
   !=======================================================================================!
   !=======================================================================================!

   !=======================================================================================!
   !  SUBROUTINE: UPDATE_PLC
   !> \breif update percentage loss of xylem conductance using daily minimum leaf psi
   !> \details This subroutine is called at daily time scale in growth_balive.f90
   !> Daily minimum leaf psi is used because upper branch, which has similar
   !> water potential as leaf, should be the most vulnerable section along the hydraulic
   !> pathway
   !=======================================================================================!
   subroutine update_plc(cpatch,ico)
   ! 核心功能： 在每日尺度上，根据植物今天经历过的正午极限最低叶片水势（dmin_leaf_psi），计算并更新植物木质部由于栓塞（Embolism）导致的导水率丧失百分比（PLC, Percentage Loss of Conductance）。
   ! 科学逻辑： * 树木在白天由于烈日暴晒、蒸腾过猛，体内的水分张力极大，水势跌入谷底（最负值）。这时候植物的维管束最容易“拉断”产生气泡（即栓塞，类似于人类的血栓），导致输水能力永久受损。
   ! 该函数计算出今天的栓塞受损程度 plc_today 后，将其根据月内天数权重（ndaysi）累加到月度变量 plc_monthly 中，用来模拟干旱对森林长期生态健康、甚至枯死（Mortality）的累积胁迫效应。
      use ed_state_vars,  only : patchtype             ! ! structure
     use ed_misc_coms   , only : current_time          & ! intent(in)
                               , simtime               ! ! structure
      use physiology_coms,only : plant_hydro_scheme    ! ! intent(in)
      use pft_coms,       only : wood_psi50            & ! intent(in)
                               , wood_Kexp             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype), target       :: cpatch
      integer        , intent(in)   :: ico
      !----- Locals    --------------------------------------------------------------------!
      real                          :: plc_today
      integer                       :: ipft
      real                          :: ndaysi
      type(simtime)                 :: lastmonth
      !------------------------------------------------------------------------------------!
      
      ! No need to update PLC if we are not tracking hydro-dynamics
      if (plant_hydro_scheme == 0) return

      call lastmonthdate(current_time,lastmonth,ndaysi)
      ipft = cpatch%pft(ico)
      
      plc_today  =  max(0., 1. - 1. /                                         &
                        (1. + (cpatch%dmin_leaf_psi(ico)                      &
                        / wood_psi50(ipft)) ** wood_Kexp(ipft)))
      cpatch%plc_monthly   (13,ico) = cpatch%plc_monthly   (13,ico) + plc_today * ndaysi

    return

   end subroutine update_plc
   !=======================================================================================!
   !=======================================================================================!



end module plant_hydro

!==========================================================================================!
!==========================================================================================!
