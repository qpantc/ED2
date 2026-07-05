!==========================================================================================!
!==========================================================================================!
!     This module contains the wrappers for the Runge-Kutta integration scheme.            !
!! rk4_driver.F90（主驱动/多线程调度）,它的主要职责是管理和调度：
!! 这段代码是陆地生态系统模型（ED2）中极其重要的物理模块 —— rk4_driver（四阶龙格-库塔驱动模块）。
!! 它的核心任务是利用 OpenMP（多线程并行），在极短的时间步长内，为森林的每一个斑块（Patch）求解:
!!    能量平衡、水分流动和碳交换的常微分方程组。
!==========================================================================================!
!==========================================================================================!
module rk4_driver

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      Main driver of short-time scale dynamics of the Runge-Kutta integrator           !
   !      for the land surface model.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_timestep(cgrid)
      use rk4_coms               , only : integration_vars           & ! structure
                                        , rk4patchtype               & ! structure
                                        , zero_rk4_patch             & ! subroutine
                                        , zero_rk4_cohort            & ! subroutine
                                        , integration_buff           ! ! intent(out)
      use ed_para_coms           , only : nthreads                   ! ! intent(in)
      use ed_state_vars          , only : edtype                     & ! structure
                                        , polygontype                & ! structure
                                        , sitetype                   ! ! structure
      use met_driver_coms        , only : met_driv_state             ! ! structure
      use grid_coms              , only : nzg                        ! ! intent(in)
      use ed_misc_coms           , only : current_time               & ! intent(in)
                                        , dtlsm                      ! ! intent(in)
      use budget_utils           , only : update_cbudget_committed   & ! function
                                        , compute_budget             ! ! function
      use soil_respiration       , only : soil_respiration_driver    ! ! sub-routine
      use stem_resp_driv         , only : stem_respiration           ! ! function
      use photosyn_driv          , only : canopy_photosynthesis      ! ! sub-routine
      use rk4_misc               , only : sanity_check_veg_energy    ! ! sub-routine
      use rk4_copy_patch         , only : copy_rk4patch_init         ! ! sub-routine
      use rk4_integ_utils        , only : copy_met_2_rk4site         ! ! sub-routine
      use update_derived_utils   , only : update_patch_derived_props ! ! sub-routine
      use plant_hydro            , only : plant_hydro_driver         ! ! sub-routine
      use therm_lib              , only : tq2enthalpy                ! ! function
      !$ use omp_lib
      implicit none
      ! 解析：定义 rk4_driver 模块并进入子程序。通过 use 导入龙格-库塔算法缓冲区（integration_buff）、
      ! 并行线程数（nthreads）以及光合、水文、物候等核心科学计算函数。!$ use omp_lib 是 OpenMP 语句，在编译时引入多线程并行库。
      ! implicit none 强制所有变量显式声明。

      !----- Arguments --------------------------------------------------------------------!
      ! cgrid 是主输入网格，带有 target 属性，允许内部指针指向它。
      type(edtype)              , target      :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)         , pointer     :: cpoly
      type(sitetype)            , pointer     :: csite
      type(met_driv_state)      , pointer     :: cmet
      ! 声明 cpoly（多边形）、csite（站点）、cmet（气象场）三大核心空间树状结构的指针。
      type(rk4patchtype)       , pointer      :: initp
      type(rk4patchtype)       , pointer      :: yscal
      type(rk4patchtype)       , pointer      :: y
      type(rk4patchtype)       , pointer      :: dydx
      type(rk4patchtype)       , pointer      :: yerr
      type(rk4patchtype)       , pointer      :: ytemp
      type(rk4patchtype)       , pointer      :: ak2
      type(rk4patchtype)       , pointer      :: ak3
      type(rk4patchtype)       , pointer      :: ak4
      type(rk4patchtype)       , pointer      :: ak5
      type(rk4patchtype)       , pointer      :: ak6
      type(rk4patchtype)       , pointer      :: ak7
      ! 声明从 initp 到 ak7 的一系列 rk4patchtype 指针，这些是龙格-库塔法（RK4）在自适应步长积分时，
      ! 存放中间斜率和误差项（如经典的 $k_1, k_2, \dots, k_7$ 步长）的高速缓存区。
      integer                                 :: ipy
      integer                                 :: isi
      integer                                 :: ipa
      integer                                 :: nsteps
      integer                                 :: imon
      real                                    :: wcurr_loss2atm
      real                                    :: ecurr_netrad
      real                                    :: ecurr_loss2atm
      real                                    :: co2curr_loss2atm
      real                                    :: wcurr_loss2drainage
      real                                    :: ecurr_loss2drainage
      real                                    :: wcurr_loss2runoff
      real                                    :: ecurr_loss2runoff
      real                                    :: co2curr_denseffect
      real                                    :: ecurr_denseffect
      real                                    :: wcurr_denseffect
      real                                    :: ecurr_prsseffect
      real                                    :: old_can_prss
      real                                    :: old_can_enthalpy
      real                                    :: old_can_temp
      real                                    :: old_can_shv
      real                                    :: old_can_co2
      real                                    :: old_can_rhos
      real                                    :: old_can_dmol
      real                                    :: patch_vels
      real                                    :: rshort_tot
      ! 声明一系列用于结算水分（wcurr）、能量（ecurr）和二氧化碳（co2curr）在森林和大气、土壤径流间通量的临时浮点数。
      integer                                 :: ibuff
      integer                                 :: npa_thread
      integer                                 :: ita
      !----- Local constants. -------------------------------------------------------------!
      logical                   , parameter   :: test_energy_sanity = .false.
      !----- Functions --------------------------------------------------------------------!
      real                      , external    :: walltime
      !------------------------------------------------------------------------------------!


      !! 宏观地理网格嵌套循环与降水统计
      polygonloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)
            !! 解析：开启两层空间大循环。指针 cpoly、csite、cmet 分别指向当前正在处理的多边形、站点和对应的时间步长气象驱动数据。


            !----- Find the number of patches per thread. ---------------------------------!
            npa_thread = ceiling(real(csite%npatches) / real(nthreads))
            ! 解析：计算每个 CPU 线程平均需要分摊处理的斑块（Patch）数量。使用 ceiling 向上取整，确保所有斑块都能被分配到线程中。
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Update the monthly rainfall.                                             !
            !------------------------------------------------------------------------------!
            imon                             = current_time%month
            cpoly%avg_monthly_pcpg(imon,isi) = cpoly%avg_monthly_pcpg(imon,isi)            &
                                             + cmet%pcpg * dtlsm
            ! 解析：获取当前模拟月份 imon。将当前步长的降水速率（cmet%pcpg）乘以步长时长（dtlsm），累加到该站点的月度总降水量累积变量中。
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Copy the meteorological variables to the rk4site structure.               !
            !------------------------------------------------------------------------------!
            call copy_met_2_rk4site(nzg,cmet%atm_ustar,cmet%atm_theiv,cmet%atm_vpdef       &
                                   ,cmet%atm_theta,cmet%atm_tmp,cmet%atm_shv,cmet%atm_co2  &
                                   ,cmet%geoht,cmet%exner,cmet%pcpg,cmet%qpcpg,cmet%dpcpg  &
                                   ,cmet%prss,cmet%rshort,cmet%rlong,cmet%par_beam         &
                                   ,cmet%par_diffuse,cmet%nir_beam,cmet%nir_diffuse        &
                                   ,cmet%geoht,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)      &
                                   ,cpoly%green_leaf_factor(:,isi),cgrid%lon(ipy)          &
                                   ,cgrid%lat(ipy),cgrid%cosz(ipy))
            ! 解析：调用函数，将气象强迫场数据（如风速、摩擦速度、水汽压亏缺、气温、湿度、CO2浓度、气压以及短波直射/散射辐射等）打包复制到 RK4 求解器的专用站点结构体中。
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !  MLO - Changed the parallel do loop to account for cases in which the number !
            !        of threads is less than the number of patches.                        !
            !------------------------------------------------------------------------------!
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(                                     &
            !$OMP  ipa,ita,initp,yscal,y,dydx,yerr,ytemp,ak2,ak3,ak4,ak5,ak6,ak7           &
            !$OMP ,patch_vels,old_can_prss,old_can_enthalpy,old_can_temp,old_can_shv       &
            !$OMP ,old_can_co2,old_can_rhos,old_can_dmol,ecurr_netrad,wcurr_loss2atm       &
            !$OMP ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage,ecurr_loss2drainage &
            !$OMP ,wcurr_loss2runoff,ecurr_loss2runoff,co2curr_denseffect,ecurr_denseffect &
            !$OMP ,wcurr_denseffect,ecurr_prsseffect,rshort_tot,nsteps)
            !! 解析：极其核心的并行控制段。
            ! 启动 OpenMP 并行循环。DEFAULT(SHARED) 规定大部分空间网格数据对所有线程可见。
            ! PRIVATE(...) 非常关键：将 RK4 缓存指针（initp 到 ak7）、斑块索引（ipa）、积分步数（nsteps）以及各项通量流定义为线程私有。防止多个 CPU 核心在并行计算不同斑块时发生数据相互覆盖（数据竞争）。
            ! 开启 threadloop，让各可用线程领走自己的缓存区编号 ibuff。
            threadloop: do ibuff=1,nthreads
               !------ Update pointers. ---------------------------------------------------!
               initp => integration_buff(ibuff)%initp
               yscal => integration_buff(ibuff)%yscal
               y     => integration_buff(ibuff)%y
               dydx  => integration_buff(ibuff)%dydx
               yerr  => integration_buff(ibuff)%yerr
               ytemp => integration_buff(ibuff)%ytemp
               ak2   => integration_buff(ibuff)%ak2
               ak3   => integration_buff(ibuff)%ak3
               ak4   => integration_buff(ibuff)%ak4
               ak5   => integration_buff(ibuff)%ak5
               ak6   => integration_buff(ibuff)%ak6
               ak7   => integration_buff(ibuff)%ak7
               ! 解析：将当前私有线程的指针，精准指向专门为其开辟的独立内存缓冲区 integration_buff(ibuff) 中，彻底实现并行内存隔离。
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Loop through tasks.  We don't assign contiguous blocks of patches to  !
               ! each thread because patches are sorted by age and older patches have more !
               ! cohorts and are likely to be slower.                                      !
               !---------------------------------------------------------------------------!
               taskloop: do ita=1,npa_thread
                  ! 解析：线程内部开启任务循环，开始解算分配给当前线程的斑块。
                  !------------------------------------------------------------------------!
                  !     Find out which patch to solve.  In case the number of patches      !
                  ! is not a perfect multiple of number of threads, some patch numbers     !
                  ! will exceed csite%npatches in the last iteration, in which we can      !
                  ! terminate the loop.                                                    !
                  !------------------------------------------------------------------------!
                  ipa = ibuff + (ita - 1) * nthreads
                  if (ipa > csite%npatches) exit taskloop
                  ! 解析：精妙的负载均衡（Load Balancing）算法。
                  ! ED2 中的斑块是按年龄排序的，老斑块里植物个体（Cohorts）多，计算极慢；新斑块计算极快。如果采用常规的按块分配（例如线程 1 算 1-10 号，线程 2 算 11-20 号），会导致线程 1 被老斑块累死，线程 2 早早闲置。
                  ! 这里采用交错式分配：若有 4 个线程，线程 1 算 1, 5, 9 号；线程 2 算 2, 6, 10 号。这样完美地把计算量巨大的老斑块均匀分摊给每个 CPU。如果计算出的斑块号 ipa 超出了该站点实际的总斑块数，则安全跳出。
                  !------------------------------------------------------------------------!

                  !----- Reset all buffers to zero, as a safety measure. ------------------!
                  call zero_rk4_patch(initp)
                  call zero_rk4_patch(yscal)
                  call zero_rk4_patch(y)
                  call zero_rk4_patch(dydx)
                  call zero_rk4_patch(yerr)
                  call zero_rk4_patch(ytemp)
                  call zero_rk4_patch(ak2)
                  call zero_rk4_patch(ak3)
                  call zero_rk4_patch(ak4)
                  call zero_rk4_patch(ak5)
                  call zero_rk4_patch(ak6)
                  call zero_rk4_patch(ak7)
                  call zero_rk4_cohort(initp)
                  call zero_rk4_cohort(yscal)
                  call zero_rk4_cohort(y)
                  call zero_rk4_cohort(dydx)
                  call zero_rk4_cohort(yerr)
                  call zero_rk4_cohort(ytemp)
                  call zero_rk4_cohort(ak2)
                  call zero_rk4_cohort(ak3)
                  call zero_rk4_cohort(ak4)
                  call zero_rk4_cohort(ak5)
                  call zero_rk4_cohort(ak6)
                  call zero_rk4_cohort(ak7)
                  ! 解析：出于数值计算的高安全性要求，在正式开跑前，调用子程序把当前斑块和个体层面的所有积分历史缓存全部强制清零。
                  !------------------------------------------------------------------------!

                  !----- Get velocity for aerodynamic resistance. -------------------------!
                  if (csite%can_theta(ipa) < cmet%atm_theta) then
                     patch_vels = cmet%vels_stab
                  else
                     patch_vels = cmet%vels_unstab
                  end if
                  ! 解析：比较林冠空气的虚位温（can_theta）与大气背景虚位温（atm_theta）。
                  ! 如果林冠更冷，说明空气层结稳定，采用稳定风速通量传输系数（vels_stab）；
                  ! 反之采用不稳定系数（vels_unstab），用以精确解算空气动力学阻力。
                  !------------------------------------------------------------------------!


                  !----- Save the previous thermodynamic state. ---------------------------!
                  old_can_prss     = csite%can_prss(ipa)
                  old_can_enthalpy = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa)    &
                                                ,.true.)
                  old_can_temp     = csite%can_temp(ipa)
                  old_can_shv      = csite%can_shv (ipa)
                  old_can_co2      = csite%can_co2 (ipa)
                  old_can_rhos     = csite%can_rhos(ipa)
                  old_can_dmol     = csite%can_dmol(ipa)
                  ! 解析：备份当前步长初始时刻林冠的物理状态（气压、温度、比湿、CO2浓度、密度等）。
                  ! 其中利用 tq2enthalpy 函数将温度和比湿换算为热力学焓（Enthalpy），这是能量守恒的核心代数变量。
                  !------------------------------------------------------------------------!



                  !----- Find incoming radiation used by the radiation driver. ------------!
                  if (cpoly%nighttime(isi)) then
                     rshort_tot = 0.0
                  else
                     rshort_tot = cmet%par_beam * csite%fbeam(ipa) + cmet%par_diffuse      &
                                + cmet%nir_beam * csite%fbeam(ipa) + cmet%nir_diffuse
                  end if
                  ! 解析：计算当前斑块接收到的总短波辐射 rshort_tot。如果是黑夜则直接为 0；
                  ! 如果是白天，则将光合有效辐射（PAR）和近红外辐射（NIR）的直射分量（乘以地形消光系数 fbeam）与散射分量全部加总。
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !      Test whether temperature and energy are reasonable.               !
                  !------------------------------------------------------------------------!
                  if (test_energy_sanity) then
                     call sanity_check_veg_energy(csite,ipa)
                  end if
                  ! 解析：Debug 检查。若开启安全测试，调用函数检查当前植被能量、温度是否在地球常识范围内
                  ! （防止产生诸如绝对零度或几万度等数值发散导致的诡异数据）。
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Get plant water flow driven by plant hydraulics.  This must be     !
                  ! placed before canopy_photosynthesis, because plant_hydro_driver needs  !
                  ! fs_open from the previous timestep.                                    !
                  !------------------------------------------------------------------------!
                  call plant_hydro_driver(csite,ipa,cpoly%ntext_soil(:,isi))
                  !! 解析：解算步骤 1：
                  !! 驱动植物水动力学模型。利用上一步气孔开度，解算水分如何从各层土壤通过根系、木质部管道源源不断地泵向叶片。
                  !------------------------------------------------------------------------!


                  !----- Get photosynthesis, stomatal conductance, and transpiration. -----!
                  call canopy_photosynthesis(csite,cmet,nzg,ipa,ibuff                      &
                                            ,cpoly%ntext_soil(:,isi)                       &
                                            ,cpoly%leaf_aging_factor(:,isi)                &
                                            ,cpoly%green_leaf_factor(:,isi))
                  !! 解析：解算步骤 2：调用核心林冠光合作用子程序。计算当前斑块内所有植物的光合速率、气孔导度、以及由于蒸腾作用散发的水分。
                  !------------------------------------------------------------------------!

                  !----- Compute stem respiration. ----------------------------------------!
                  call stem_respiration(csite,ipa)
                  !------------------------------------------------------------------------!


                  !----- Compute root and heterotrophic respiration. ----------------------!
                  call soil_respiration_driver(csite,ipa,nzg,cpoly%ntext_soil(:,isi))
                  !! 解析：解算步骤 3 & 4：计算植物群落的树干呼吸（维持呼吸）以及由土壤微生物和根系贡献的土壤呼吸（异养+自养呼吸）。
                  !------------------------------------------------------------------------!


                  !----- Update the committed carbon change pool. -------------------------!
                  call update_cbudget_committed(csite,ipa)
                  !! 解析：将上述计算得到的碳汇变动更新至“已承诺碳库”中（周转记账）。
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Set up the integration patch.                                      !
                  !------------------------------------------------------------------------!
                  call copy_rk4patch_init(csite,ipa,ibuff,initp,patch_vels                 &
                                         ,old_can_enthalpy,old_can_rhos,old_can_dmol       &
                                         ,ecurr_prsseffect)
                  !!解析：把当前斑块的所有初始物理量及刚刚算出来的生理通量（光合、呼吸、水流）一次性打包打包进积分器入口结构体 initp 中。
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    This is the driver for the integration process...                   !
                  !------------------------------------------------------------------------!
                  call integrate_patch_rk4(csite,initp,ipa,isi,ibuff                       &
                                          ,cpoly%nighttime(isi),wcurr_loss2atm             &
                                          ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm    &
                                          ,wcurr_loss2drainage,ecurr_loss2drainage         &
                                          ,wcurr_loss2runoff,ecurr_loss2runoff             &
                                          ,co2curr_denseffect,ecurr_denseffect             &
                                          ,wcurr_denseffect,nsteps)
                  !! 解析：全程序最核心的微分方程求解：
                  !! 调用下方的子程序，驱使龙格-库塔数值求解器对当前斑块的微分方程组实施时间积分（详见下文）。
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Add the number of steps into the step counter. Workload            !
                  ! accumulation is order-independent, so this can stay shared.            !
                  !------------------------------------------------------------------------!
                  cgrid%workload(13,ipy) = cgrid%workload(13,ipy) + real(nsteps)
                  ! 解析：性能监控：由于 RK4 是自适应步长，如果环境剧烈变化，单步内它可能在内部迭代解算了成百上千次（nsteps）。
                  ! 这里把实际解算步数累加到 workload 中，以便后续统计各多边形的计算工作量（负载）。
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !   Update the minimum monthly temperature, based on canopy temperature. !
                  !------------------------------------------------------------------------!
                  if (cpoly%site(isi)%can_temp(ipa) < cpoly%min_monthly_temp(isi)) then
                     cpoly%min_monthly_temp(isi) = cpoly%site(isi)%can_temp(ipa)
                  end if
                  ! 解析：生态学历史记录：如果积分出来的林冠温度刷新了本月最低纪录，
                  !  则更新本月的最低温度（用于后续植物因极端低温冻死或触发落叶的生存环境判定）。
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Update roughness and canopy depth.  This should be done after the   !
                  ! integration.                                                           !
                  !------------------------------------------------------------------------!
                  call update_patch_derived_props(csite,ipa,.false.)
                  ! 解析：积分完成后，根据最新的植被状态更新斑块的派生属性（如最新的地表粗糙度、林冠深度等）。
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Compute the residuals.                                             !
                  !------------------------------------------------------------------------!
                  call compute_budget(csite,cpoly%lsl(isi),cmet%pcpg,cmet%qpcpg            &
                                     ,rshort_tot,cmet%rlong,ipa,wcurr_loss2atm             &
                                     ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm         &
                                     ,wcurr_loss2drainage,ecurr_loss2drainage              &
                                     ,wcurr_loss2runoff,ecurr_loss2runoff                  &
                                     ,co2curr_denseffect,ecurr_denseffect,wcurr_denseffect &
                                     ,ecurr_prsseffect,cpoly%area(isi)                     &
                                     ,cgrid%cbudget_nep(ipy),old_can_prss                  &
                                     ,old_can_enthalpy,old_can_temp,old_can_shv            &
                                     ,old_can_co2,old_can_rhos,old_can_dmol)
                  !! 解析：最后一道收支核验：结合老状态和新状态，调用 compute_budget 重新清算并校验当前斑块的水分、
                  !! 能量、碳是否满足严格的物理守恒定律，并结算净生态系统生产力（NEP）。
                  !------------------------------------------------------------------------!

               end do taskloop
               !---------------------------------------------------------------------------!
            end do threadloop
            !$OMP END PARALLEL DO
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polygonloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine rk4_timestep
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the integration process.                               !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_patch_rk4(csite,initp,ipa,isi,ibuff,nighttime,wcurr_loss2atm       &
                                 ,ecurr_netrad,ecurr_loss2atm,co2curr_loss2atm             &
                                 ,wcurr_loss2drainage,ecurr_loss2drainage                  &
                                 ,wcurr_loss2runoff,ecurr_loss2runoff,co2curr_denseffect   &
                                 ,ecurr_denseffect,wcurr_denseffect,nsteps)
      !! 这个内部子程序是真正的“微分方程推手”，它只负责一件纯粹的数学任务：
      !!    把斑块的物理方程丢进积分器，并在解算完成后将新状态安全的交还给模型主框架。
      use rk4_integ_utils , only : odeint               ! ! sub-routine
      ! 进入子程序。从外部工具库导入 odeint（常微分方程积分器）。
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use rk4_coms        , only : integration_vars     & ! structure
                                 , rk4patchtype         & ! structure
                                 , zero_rk4_patch       & ! subroutine
                                 , zero_rk4_cohort      & ! subroutine
                                 , tbeg                 & ! intent(inout)
                                 , tend                 & ! intent(inout)
                                 , dtrk4i               ! ! intent(inout)
      use rk4_copy_patch  , only : initp2modelp         ! ! sub-routine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      type(rk4patchtype)    , target      :: initp
      integer               , intent(in)  :: ipa
      integer               , intent(in)  :: isi
      integer               , intent(in)  :: ibuff
      logical               , intent(in)  :: nighttime
      real                  , intent(out) :: wcurr_loss2atm
      real                  , intent(out) :: ecurr_netrad
      real                  , intent(out) :: ecurr_loss2atm
      real                  , intent(out) :: co2curr_loss2atm
      real                  , intent(out) :: wcurr_loss2drainage
      real                  , intent(out) :: ecurr_loss2drainage
      real                  , intent(out) :: wcurr_loss2runoff
      real                  , intent(out) :: ecurr_loss2runoff
      real                  , intent(out) :: co2curr_denseffect
      real                  , intent(out) :: ecurr_denseffect
      real                  , intent(out) :: wcurr_denseffect
      integer               , intent(out) :: nsteps
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Zero the canopy-atmosphere flux values.  These values are updated every dtlsm, !
      ! so they must be zeroed at each call.                                               !
      !------------------------------------------------------------------------------------!
      initp%upwp = 0.d0
      initp%tpwp = 0.d0
      initp%qpwp = 0.d0
      initp%cpwp = 0.d0
      initp%wpwp = 0.d0
      ! 解析：将林冠与大气之间的各项湍流交换通量项（动量通量 upwp、热通量 tpwp、潜热/水分通量 qpwp、二氧化碳通量 cpwp、液态水通量 wpwp）
      ! 在积分前初始化为双精度 0.0 (0.d0)。

      !----- Go into the ODE integrator. --------------------------------------------------!
      call odeint(csite,ipa,isi,ibuff,nsteps)
      ! 解析：数学运算暴风眼：正式调用 odeint（常微分方程求解器）。
      ! 模型会在内部通过高阶数值矩阵运算，对林冠和土壤的温度、湿度、二氧化碳偏微分方程进行强行演进，
      ! 并将最终演进出的正确状态写回 initp 中。nsteps 传回实际迭代次数。
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Normalize canopy-atmosphere flux values.  These values are updated every      !
      ! dtlsm, so they must be normalized every time.                                      !
      !------------------------------------------------------------------------------------!
      initp%upwp = initp%can_rhos * initp%upwp * dtrk4i
      initp%tpwp = initp%can_rhos * initp%tpwp * dtrk4i
      initp%qpwp = initp%can_rhos * initp%qpwp * dtrk4i
      initp%cpwp = initp%can_dmol * initp%cpwp * dtrk4i
      initp%wpwp = initp%can_rhos * initp%wpwp * dtrk4i
      ! 解析：通量物理归一化。积分器传出来的通量属于时间和密度的累加项。
      ! 这里将它们乘以林冠空气密度（can_rhos 或摩尔密度 can_dmol），
      ! 再乘以当前实际总积分时间步长的倒数（dtrk4i = $1 / \Delta t$）。
      ! 经过这几行运算，这些变量从“无量纲的数学积分积分值”被还原成了具有明确物理意义的、真实的、单位时间标准物理通量。


      !------------------------------------------------------------------------------------!
      ! Move the state variables from the integrated patch to the model patch.             !
      !------------------------------------------------------------------------------------!
      call initp2modelp(tend-tbeg,initp,csite,ipa,nighttime,wcurr_loss2atm,ecurr_netrad    &
                       ,ecurr_loss2atm,co2curr_loss2atm,wcurr_loss2drainage                &
                       ,ecurr_loss2drainage,wcurr_loss2runoff,ecurr_loss2runoff            &
                       ,co2curr_denseffect,ecurr_denseffect,wcurr_denseffect)
      ! 解析：大功告成。
      ! 调用 initp2modelp 转换函数，计算时间跨度 tend-tbeg，把刚刚在临时缓冲区 initp 里算出来的、热乎的、
      ! 完全满足守恒律的全新状态变量和各路径通量数据，安全、完整地复制回 ED2 模型主网格系统（csite）的实际斑块中。
      !------------------------------------------------------------------------------------!


      return
   end subroutine integrate_patch_rk4
   !=======================================================================================!
   !=======================================================================================!
end module rk4_driver

!==========================================================================================!
!==========================================================================================!
