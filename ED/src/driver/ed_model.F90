!==========================================================================================!
!==========================================================================================!
!    This file contains the main driver for ED2.                                           !
!------------------------------------------------------------------------------------------!
! SUBROUTINE: ED_MODEL
!
!> \brief   Begins, updates, and outputs results from ecosystem simulation.
!> \details Coordinates meteorological forcing data, time-stepping, radiation, vegetation
!>          dynamics, hydraulogy, and the HDF5 output subsystem.
!> \author  Translated from ED1 by David Medvigy, Ryan Knox and Marcos Longo
!! 这段代码是 ED2 模型的核心逻辑纽带 —— ed_model() 子程序。
!! 它包含了陆地生态系统模型最核心的主时间步长循环（Main Time Loop）。
!! 在这段代码里，模型会反复经历：
!! 计算气象 ➔ 求解辐射 ➔ 生物物理积分（光合/呼吸） ➔ 推进时间 ➔ 计算植被动态（生长/死亡） ➔ 更新水文 ➔ 输出 HDF5 数据。
!------------------------------------------------------------------------------------------!
subroutine ed_model()
   use ed_misc_coms        , only : simtime                     & ! structure
                                  , ivegt_dynamics              & ! intent(in)
                                  , integration_scheme          & ! intent(in)
                                  , current_time                & ! intent(in)
                                  , frqfast                     & ! intent(in)
                                  , frqstate                    & ! intent(in)
                                  , out_time_fast               & ! intent(in)
                                  , dtlsm                       & ! intent(in)
                                  , ifoutput                    & ! intent(in)
                                  , isoutput                    & ! intent(in)
                                  , iqoutput                    & ! intent(in)
                                  , itoutput                    & ! intent(in)
                                  , iooutput                    & ! intent(in)
                                  , restore_file                & ! intent(in)
                                  , frqsum                      & ! intent(in)
                                  , unitfast                    & ! intent(in)
                                  , unitstate                   & ! intent(in)
                                  , imontha                     & ! intent(in)
                                  , iyeara                      & ! intent(in)
                                  , outstate                    & ! intent(in)
                                  , outfast                     & ! intent(in)
                                  , nrec_fast                   & ! intent(in)
                                  , nrec_state                  & ! intent(in)
                                  , runtype                     & ! intent(in)
                                  , month_yrstep                & ! intent(in)
                                  , writing_dail                & ! intent(in)
                                  , writing_mont                & ! intent(in)
                                  , writing_dcyc                & ! intent(in)
                                  , writing_eorq                & ! intent(in)
                                  , writing_long                & ! intent(in)
                                  , writing_year                ! ! intent(in)
   use ed_init             , only : remove_obstime              & ! sub-routine
                                  , is_obstime                  ! ! sub-routine
   use grid_coms           , only : ngrids                      & ! intent(in)
                                  , istp                        & ! intent(in)
                                  , time                        & ! intent(in)
                                  , timmax                      ! ! intent(in)
   use ed_state_vars       , only : edgrid_g                    & ! intent(in)
                                  , edtype                      & ! intent(in)
                                  , patchtype                   & ! intent(in)
                                  , filltab_alltypes            & ! intent(in)
                                  , filltables                  ! ! intent(in)
   use rk4_driver          , only : rk4_timestep                ! ! sub-routine
   use rk4_coms            , only : integ_err                   & ! intent(in)
                                  , integ_lab                   & ! intent(in)
                                  , record_err                  & ! intent(inout)
                                  , print_detailed              & ! intent(inout)
                                  , nerr                        & ! intent(in)
                                  , errmax_fout                 & ! intent(in)
                                  , sanity_fout                 & ! intent(in)
                                  , alloc_integ_err             & ! subroutine
                                  , assign_err_label            & ! subroutine
                                  , reset_integ_err             ! ! subroutine
   use ed_node_coms        , only : mynum                       & ! intent(in)
                                  , nnodetot                    ! ! intent(in)
   use mem_polygons        , only : n_ed_region                 & ! intent(in)
                                  , n_poi                       ! ! intent(in)
   use consts_coms         , only : day_sec                     ! ! intent(in)
   use average_utils       , only : update_ed_yearly_vars       & ! sub-routine
                                  , zero_ed_dmean_vars          & ! sub-routine
                                  , zero_ed_mmean_vars          & ! sub-routine
                                  , zero_ed_qmean_vars          & ! sub-routine
                                  , zero_ed_fmean_vars          & ! sub-routine
                                  , integrate_ed_fmean_met_vars & ! sub-routine
                                  , zero_ed_yearly_vars         ! ! sub-routine
   use edio                , only : ed_output                   ! ! sub-routine
   use ed_met_driver       , only : read_met_drivers            & ! sub-routine
                                  , update_met_drivers          ! ! sub-routine
   use euler_driver        , only : euler_timestep              ! ! sub-routine
   use heun_driver         , only : heun_timestep               ! ! sub-routine
   use hybrid_driver       , only : hybrid_timestep             ! ! sub-routine
   use lsm_hyd             , only : updateHydroParms            & ! sub-routine
                                  , calcHydroSubsurface         & ! sub-routine
                                  , calcHydroSurface            & ! sub-routine
                                  , writeHydro                  ! ! sub-routine
   use radiate_driver      , only : canopy_radiation            ! ! sub-routine
   use rk4_integ_utils     , only : initialize_rk4patches       & ! sub-routine
                                  , initialize_misc_stepvars    ! ! sub-routine
   use stable_cohorts      , only : flag_stable_cohorts         ! ! sub-routine
   use update_derived_utils, only : update_model_time_dm        ! ! sub-routine
   use budget_utils        , only : ed_init_budget              ! ! intent(in)
   use vegetation_dynamics , only : veg_dynamics_driver         ! ! sub-routine
   use ed_type_init        , only : ed_init_viable              ! ! sub-routine
   use soil_respiration    , only : zero_litter_inputs          ! ! sub-routine
   implicit none
   !----- Common blocks. ------------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !----- Local variables. ----------------------------------------------------------------!
   type(simtime)      :: daybefore
   character(len=28)  :: fmthead
   character(len=32)  :: fmtcntr
   integer            :: ifm
   integer            :: nn
   integer            :: ndays
   integer            :: dbndays
   integer            :: obstime_idx
   logical            :: analysis_time
   logical            :: observation_time
   logical            :: new_day
   logical            :: new_month
   logical            :: new_year
   logical            :: history_time
   logical            :: dcycle_time
   logical            :: annual_time
   logical            :: mont_analy_time
   logical            :: dail_analy_time
   logical            :: dcyc_analy_time
   logical            :: reset_time
   logical            :: past_one_day
   logical            :: past_one_month
   logical            :: printbanner
   logical            :: veget_dyn_on
   real               :: wtime_start
   real               :: t1
   real               :: wtime1
   real               :: wtime2
   real               :: t2
   real               :: wtime_tot
   real               :: dbndaysi
   real               :: gr_tfact0
   !----- Local variables (MPI only). -----------------------------------------------------!
#if defined(RAMS_MPI)
   integer            :: ierr
#endif
   !----- Local constants. ----------------------------------------------------------------!
   logical         , parameter :: whos_slow=.false. ! This will print out node numbers
                                                    !    during synchronization, so you
                                                    !    can find out which node is the
                                                    !    slow one
   !----- String for output format. -------------------------------------------------------!
   character(len=26), parameter :: fmtrest = '(i4.4,2(1x,i2.2),1x,2i2.2)'
   !----- External functions. -------------------------------------------------------------!
   real    , external :: walltime ! Wall time
   integer , external :: num_days ! Number of days in the current month
   ! 声明一系列局部变量。包括用来记录昨天日期信息的自定义结构体 daybefore；
   ! 用来动态拼接文件输出格式的字符串 fmthead 和 fmtcntr；
   ! 大量控制数据输出时机的布尔逻辑标记（如是否为新的一天/月/年）；
   ! 以及用于监控运行效率的 CPU/墙钟耗时浮点数。
   ! whos_slow 是一项调试开关，设为 .true. 时能监控并行状态下哪个计算节点遭遇了性能瓶颈。

   !---------------------------------------------------------------------------------------!

   past_one_day   = .false.
   past_one_month = .false.
   filltables     = .false.

   ! 重置这几个关键的时间跨越标记与数据表填充标记，为接下来的时间循环做准备。

   !----- Print the hour banner only for regional runs that aren't massively parallel. ----!
   printbanner = n_ed_region > 0 .and. edgrid_g(1)%npolygons > 50 .and. mynum == 1
   ! 设定“打印运行时状态横幅”的触发条件：必须是区域运行（n_ed_region > 0），
   ! 多边形网格数大于 50，且只有 1 号进程（主节点）有权打印，防止数百个并行节点同时向终端写数据引发画面崩溃。

   !----- Run with vegetation dynamics turned on?  ----------------------------------------!
   veget_dyn_on = ivegt_dynamics == 1
   ! 读取全局配置。若 ivegt_dynamics 的值为 1，则布尔变量 veget_dyn_on 设为真，激活生态学上的植物生死演替计算。

   wtime_start=walltime(0.)
   istp = 0
   ! 调用外部函数 walltime(0.) 抓取绝对起点时间，重置当前总模拟步数计数器 istp。

   !---------------------------------------------------------------------------------------!
   !     If we are going to record the integrator errors, here is the time to open it for  !
   ! the first time and write the header.  But just before we do it, we check whether this !
   ! is a single POI run, the only case where we will allow this recording.                !
   !---------------------------------------------------------------------------------------!
   record_err     = record_err     .and. n_ed_region == 0 .and. n_poi == 1
   print_detailed = print_detailed .and. n_ed_region == 0 .and. n_poi == 1
   ! 这是一层安全锁。规定只有在非区域模拟、且仅有一个单点（POI）测试运行的前提下，
   ! 才允许开启数值积分器（如龙格库塔 RK4）的误差记录和详细日志。多点并行时若写这个会产生巨大的磁盘 I/O 负担。
   if(record_err) then
      ! 若开启误差记录，分配内存，初始化误差数组，并为不同的方程误差打上文本标签（如温度、水分、碳通量）。
      !----- Initialise the error structures. ---------------------------------------------!
      call alloc_integ_err()
      call reset_integ_err()
      call assign_err_label()

      !----- Define the formats for both the header and the actual output. ----------------!
      write(fmthead,fmt='(a,i3.3,a)')  '(a4,1x,2(a3,1x),',nerr,'(a13,1x))'
      write(fmtcntr,fmt='(a,i3.3,a)')  '(i4.4,1x,2(i3.2,1x),',nerr,'(i13,1x))'
      ! 动态构建 Fortran 的格式化字符串。将总误差变量数 nerr 动态拼接进格式中，形成类似 (a4,1x,2(a3,1x), 025(a13,1x)) 的格式定义。

      open  (unit=77,file=trim(errmax_fout),form='formatted',status='replace')
      write (unit=77,fmt=fmthead) 'YEAR','MON','DAY',(integ_lab(nn),nn=1,nerr)
      close (unit=77,status='keep')

      open  (unit=78,file=trim(sanity_fout),form='formatted',status='replace')
      write (unit=78,fmt=fmthead) 'YEAR','MON','DAY',(integ_lab(nn),nn=1,nerr)
      close (unit=78,status='keep')
      ! 打开指定的误差输出文件和合理性检查文件（如果已存在则直接覆盖 status='replace'），写入由年、月、日和各变量标签组成的表格表头，随后暂时关闭文件。
   end if

   out_time_fast     = current_time
   out_time_fast%month = -1
   !  同步高频输出时间结构体，并人为将月份设为 -1 作为初始未激活状态。

   !---------------------------------------------------------------------------------------!
   !      If this is not a history restart, then zero out the long term diagnostics.       !
   !---------------------------------------------------------------------------------------!
   select case (trim(runtype))
   case ('HISTORY')
      continue
   case default

      do ifm=1,ngrids
         if (writing_long) call zero_ed_dmean_vars(edgrid_g(ifm))
         if (writing_eorq) call zero_ed_mmean_vars(edgrid_g(ifm))
         if (writing_dcyc) call zero_ed_qmean_vars(edgrid_g(ifm))
      end do

      !----- Long-term dynamics structure. ------------------------------------------------!
      do ifm=1,ngrids
         call update_ed_yearly_vars(edgrid_g(ifm))
      end do
      !------------------------------------------------------------------------------------!
   end select
   !! 检查运行类型。如果是 HISTORY（断点恢复运行），什么都不做，直接继承上一次的均值积分；
   !! 如果是全新运行（default），则遍历网格，将其日均值、月均值、日变化均值统计数组全部清零，并初始化年尺度变量。
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      The fast analysis is always reset, including history runs.                       !
   !---------------------------------------------------------------------------------------!
   do ifm=1,ngrids
      call zero_ed_fmean_vars(edgrid_g(ifm))
   end do
   ! 不论是否为历史续跑，高频/快速（Fast）分析输出的统计数组在每轮开始前必须强制清零。
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Allocate memory to the integration patch, Euler now utilises the RK4 buffers too.  !
   !---------------------------------------------------------------------------------------!
   ! 为积分算法（现在 Euler 算法也会借用 RK4 的高速缓存空间）分配斑块级的数组缓冲区，并初始化步长控制变量。
   call initialize_rk4patches(.true.)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Initialize some stepping variables.                                                !
   !---------------------------------------------------------------------------------------!
   call initialize_misc_stepvars()
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Here we must initialise or reset a group of variables.                           !
   ! 1.  Variable is_viable must be set to .true..  This variable is not saved in          !
   !     ed_init_history, and the default is .false., which would eliminate all cohorts.   !
   ! 2.  In the case of a initial simulation, we must reset all budget fluxes and set all  !
   !     budget stocks.  This should not be done in HISTORY initialisation, all variables  !
   !     should be read from history.                                                      !
   ! 3.  Litter inputs must be reset in the HISTORY initialisation, in case the history    !
   !     file is at midnight UTC (daily time step).  These variables are normally reset    !
   !     after writing the output so they are meaningful in the output.  Because history   !
   !     files are written before the inputs are reset, the inputs would be double-counted !
   !     in the second day of simulation.                                                  !
   !---------------------------------------------------------------------------------------!
   select case (trim(runtype))
   case ('INITIAL')
      do ifm=1,ngrids
         call ed_init_budget(edgrid_g(ifm),.true.)
         call ed_init_viable(edgrid_g(ifm))
      end do
   ! 若为全新启动：遍历网格，初始化系统内部的水分与能量收支平衡账本（ed_init_budget），
   ! 同时执行 ed_init_viable 将系统内植物同龄群（Cohorts）的生命力标记 is_viable 强行激活为 .true.。
   ! （注：如果不做此操作，这些群落默认是 .false.，会被系统判定为死群并从内存中直接抹除）。
   case ('HISTORY')
      new_day         = current_time%time < dtlsm
      do ifm=1,ngrids
         call flag_stable_cohorts(edgrid_g(ifm),.true.)
         call ed_init_viable(edgrid_g(ifm))      
         if (new_day) then
            call zero_litter_inputs(edgrid_g(ifm))
         end if
      end do
   ! 若为断点续跑：如果恢复时的系统时间刚好处于一天的第一个步长内（current_time%time < dtlsm），
   ! 则执行 zero_litter_inputs 清空枯落物输入流。
   ! 这行设计非常巧妙：历史恢复文件是在每天结束、数据清零前写入的。
   ! 如果新模拟的第一天不清零，之前的凋落物就会在第二天被重复计算。
   end select
   !---------------------------------------------------------------------------------------!


   if (ifoutput /= 0) call h5_output('INST')

   if (isoutput /= 0) then
      select case (trim(runtype))
      case ('INITIAL')
         call h5_output('HIST')
      case ('HISTORY')
         call h5_output('CONT')
      end select
   end if

   if (writing_year ) call h5_output('YEAR')
   ! 在正式跨入时间积分循环前，将当前的初始状态作为“第 0 步”数据写进 HDF5 文件中。
   ! 包括瞬时输出（INST）、历史状态（HIST/CONT）和年输出（YEAR）。


   ! 进入模型主时间积分循环 (timestep)。在这个循环里，模型会不断推进时间，更新状态，并根据设定的频率输出数据。
   !----- Start the timesteps. ------------------------------------------------------------!
   if (mynum == 1) write(unit=*,fmt='(a)') ' === Time integration starts (model) ==='
   ! 主节点打印模拟启动信息。开启有标签的 do while 循环，只要当前模拟累计时间 time 小于设定的最大上限 timmax，
   ! 循环就会不断向前推进。每进一次循环，总步数 istp 自增 1。

   !! 每次进循环，步数就加 1
   timestep: do while (time < timmax)

      istp = istp + 1

      !------------------------------------------------------------------------------------!
      !   CPU timing information & model timing information.                               !
      !------------------------------------------------------------------------------------!
      call timing(1,t1)
      wtime1=walltime(wtime_start)
      ! 在每个步长开始前，立刻抓取当前的 CPU 时间 t1 和墙钟时间 wtime1，用于后续计算单步耗时。
      !------------------------------------------------------------------------------------!

      if (current_time%time < dtlsm .and. mynum == 1) then
           write (unit=*,fmt='(a,3x,2(i2.2,a),i4.4,a,3(i2.2,a))')                          &
              ' - Simulating:',current_time%month,'/',current_time%date,'/'                &
                              ,current_time%year,' ',current_time%hour,':'                 &
                              ,current_time%min,':',current_time%sec,' UTC'
      end if
      !! 如果当前正处于一天的第一个时间步，主节点会在终端优雅地打印当前模拟推进到了现实世界中的哪一年、哪一月、
      !! 哪一日以及具体时间（UTC时间）。

      !----- Define which cohorts are to be solved prognostically. ------------------------!
      do ifm=1,ngrids
         call flag_stable_cohorts(edgrid_g(ifm),.false.)
         ! 遍历网格，筛选并标记哪些植物群落处于稳定状态（可以进行预测性解算）。
      end do
      !------------------------------------------------------------------------------------!


      !----- Solve the radiation profile. -------------------------------------------------!
      do ifm=1,ngrids
          call canopy_radiation(edgrid_g(ifm))
          !! 核心物理步骤 1：遍历网格，调用冠层辐射传输模型。计算太阳辐射在各层植物叶片、树冠之间的吸收、反射和透射。
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     At this point, all meteorologic driver data for the land surface model has     !
      ! been updated for the current timestep.  Perform the time average for the output    !
      ! diagnostic.                                                                        !
      !------------------------------------------------------------------------------------!
      do ifm=1,ngrids
         call integrate_ed_fmean_met_vars(edgrid_g(ifm))
         ! 当前步长的气象驱动数据已就绪，调用此子程序在时间层面上对高频气象诊断通量进行累加积分。
      end do
      !------------------------------------------------------------------------------------!


      !----- Solve the photosynthesis and biophysics. -------------------------------------!
      select case (integration_scheme)
      !! 核心物理步骤 2：这是 ED2 处理能量和气体交换（光合作用、气孔导度、蒸腾作用）的核心积分器选择逻辑。
      ! 根据用户的配置 integration_scheme（默认是1），选择不同的数值微分方程解算器来向前推进一个主步长（dtlsm）：
      case (0)
         do ifm=1,ngrids
            call euler_timestep(edgrid_g(ifm)) ! 一阶显式欧拉法（euler，速度快精度低）。
         end do
      case (1)
         do ifm=1,ngrids
            call rk4_timestep(edgrid_g(ifm)) ! 四阶龙格库塔法（rk4，速度慢精度高）。
            !! rk4_driver.F90 中的 rk4_timestep() 是 ED2 的默认积分器，采用经典的四阶龙格库塔方法，能在保持数值稳定性的同时提供较高的精度，但计算成本较高。
         end do
      case (2)
         do ifm=1,ngrids
            call heun_timestep(edgrid_g(ifm)) ! 二阶显式 Heun 法（heun，速度适中精度适中）。
         end do
      case (3)
         do ifm=1,ngrids
            call hybrid_timestep(edgrid_g(ifm)) ! 混合方法（hybrid，针对不同变量选择不同的积分器以优化效率和稳定性）。
         end do
      end select
      !------------------------------------------------------------------------------------!


      ! 8. 更新模拟时间与时间节点检测
      !------------------------------------------------------------------------------------!
      !     Update the model time.                                                         !
      !------------------------------------------------------------------------------------!
      time=time+dble(dtlsm)
      call update_model_time_dm(current_time, dtlsm)
      ! 物理方程解算完毕，将累计时间 time 加上当前步长。
      ! 同时更新人类可读的时间结构体 current_time（自动处理进位，如秒变成微、日变成月等）。
      !------------------------------------------------------------------------------------!


      !----- Check whether it is some special time... -------------------------------------!
      new_day         = current_time%time < dtlsm
      if (.not. past_one_day .and. new_day) past_one_day=.true.

      new_month       = current_time%date == 1  .and. new_day
      if (.not. past_one_month .and. new_month) past_one_month=.true.

      new_year        = current_time%month == month_yrstep .and. new_month
      mont_analy_time = new_month .and. writing_mont
      dail_analy_time = new_day   .and. writing_dail
      dcyc_analy_time = new_month .and. writing_dcyc
      reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
      annual_time     = new_year .and. writing_year
      ! 这一大段是一套极其紧密的时钟边沿触发逻辑。
      ! 通过检测当前秒数是否小于一个步长，来判断是否跨入了新的一天（new_day）；
      ! 同理检测是否为新一月（new_month）或新一年（new_year），并据此举起各类输出和重置的布尔逻辑“信号旗”。
      !------------------------------------------------------------------------------------!



      !----- Check whether this is time to write fast analysis output or not. -------------!
      ! 根据高频输出的单位类型（秒、月或年），通过取模运算（mod）判断当前时刻是否正好踩在了用户设定的高频数据输出周期（frqfast）
      ! 或日变化周期（iqoutput）的点上。如果是，则将 analysis_time 或 dcycle_time 激活为真。
      select case (unitfast)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         analysis_time   = mod(current_time%time, frqfast) < dtlsm .and.                   &
                           (ifoutput /= 0 .or. itoutput /= 0 .or. iooutput /= 0)
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      case (2)   !----- Months, analysis time is at the new month -------------------------!
         analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /=0) .and.         &
                           mod(real(12+current_time%month-imontha),frqfast) == 0.
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      case (3) !----- Year, analysis time is at the same month as initial time ------------!
         analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /= 0) .and.        &
                           current_time%month == imontha .and.                             &
                           mod(real(current_time%year-iyeara),frqfast) == 0.
         dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
      end select
      !------------------------------------------------------------------------------------!



      !----- Check whether it is an observation time --------------------------------------!
      ! 观测时间（Observation Time）匹配。如果用户提供了一份零散的、不规则的历史观测事件时间表，程序在此处调用 is_obstime 检查当前时间是否与表中的某个记录吻合。如果吻合（observation_time = .true.），数据输出后，将该时刻从待检查列表中移除（remove_obstime）。
      if (iooutput == 0 .or. unitfast /= 0) then
         !------ Observation_time is not used when unitfast /= 0 or iooutput is 0. --------!
         observation_time = .false. 
         !---------------------------------------------------------------------------------!
      else
         !------ check whether it is the observation time. --------------------------------!
         call is_obstime(current_time%year,current_time%month,current_time%date            &
                        ,current_time%time,observation_time,obstime_idx)
         !---------------------------------------------------------------------------------!

         !------ Get rid of the obstime record if observation_time is true. ---------------!
         if (observation_time) then
            call remove_obstime(obstime_idx)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !----- Check whether this is time to write restart output or not. -------------------!
      ! 同理，通过取模运算，判断当前步长是否达到了输出全状态恢复文件/重启历史文件（HIST）的指定周期（frqstate）。
      select case(unitstate)
      case (0,1) !----- Now both are in seconds -------------------------------------------!
         history_time   = mod(current_time%time, frqstate) < dtlsm .and. isoutput /= 0
      case (2)   !----- Months, history time is at the new month --------------------------!
         history_time   = new_month .and. isoutput /= 0 .and.                              &
                          mod(real(12+current_time%month-imontha),frqstate) == 0.
      case (3) !----- Year, history time is at the same month as initial time -------------!
         history_time   = new_month .and. isoutput /= 0 .and.                              &
                          current_time%month == imontha .and.                              &
                          mod(real(current_time%year-iyeara),frqstate) == 0.
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Update nrec_fast and nrec_state if it is a new month and outfast/outstate are   !
      ! monthly and frqfast/frqstate are daily or by seconds.                              !
      !------------------------------------------------------------------------------------!
      ! 如果进入了新的一月，调用外部函数 num_days 计算这个月实际有多少天（自动处理闰年）。如果用户配置中指定了按月切分文件（格式标记为 -2.），这里会实时动态计算并更新本月内高频和状态输出的数据记录总总条数（nrec_fast/nrec_state），用于初始化 HDF5 的数据集维度。
      if (new_month) then
         ndays = num_days(current_time%month,current_time%year)
         if (outfast  == -2.) nrec_fast  = ndays*ceiling(day_sec/frqfast)
         if (outstate == -2.) nrec_state = ndays*ceiling(day_sec/frqstate)
      end if
      !------------------------------------------------------------------------------------!


      !! 10. 日尺度触发：误差写入与群落动态学（Demographics）
      !----- Check if this is the beginning of a new simulated day. -----------------------!
      if (new_day) then
         if (record_err) then
         ! 如果今天是个新日期的起点，并且开启了 Debug 机制，则以追加模式（access='append'）重新打开第 77、78 号文件，
         ! 把昨天一整天数值积分器累计的最大误差值、合理性指标写入磁盘，随后清空误差缓存，重新开始计算。
            open (unit=77,file=trim(errmax_fout),form='formatted',access='append'          &
                 ,status='old')
            write (unit=77,fmt=fmtcntr) current_time%year,current_time%month               &
                                       ,current_time%date,(integ_err(nn,1),nn=1,nerr)
            close(unit=77,status='keep')

            open (unit=78,file=trim(sanity_fout),form='formatted',access='append'          &
                 ,status='old')
            write (unit=78,fmt=fmtcntr) current_time%year,current_time%month               &
                                       ,current_time%date,(integ_err(nn,2),nn=1,nerr)
            close(unit=78,status='keep')

            call reset_integ_err()
         end if
         


         !----- Find the number of days in this month and the previous month. -------------!
         call yesterday_info(current_time,daybefore,dbndays,dbndaysi)
         ! 解析：调用函数获取昨天的时钟快照，计算昨天所在月份的总天数 dbndays。
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     This cap limits the growth rate depending on the day of the month.  This is !
         ! to account for the different time scales between heartwood growth (monthly) and !
         ! growth of the other tissues (daily).  This factor ensures that growth           !
         ! respiration is evenly distributed during the month, as opposed to have a spike  !
         ! on the second day of every month, when live tissues are growing to catch up the !
         ! allometry after heartwood biomass had increased.  This factor grows as the time !
         ! step approaches the end of the month, but the biomass increment will be the     !
         ! same every day and trees will be back on allometry be the time of the following !
         ! month in case storage is not limiting.
         !---------------------------------------------------------------------------------!
         gr_tfact0 = 1.0 / (dbndays - daybefore%date + 1)
         ! 极为优雅的生态学设计。
         ! 计算一个生物生长权重调节因子。
         ! 在 ED2 模型中，植物的边材/心材生长、碳分配是在月尺度上算的，而活体组织（叶片、细根）的生长是日尺度算的。
         ! 为了防止每个月的第二天因为全分配对齐导致植物组织出现爆发性的“畸形暴长”（Growth Spike），
         ! 这个因子会随着日期接近月底而变大，从而确保生长的呼吸消耗能被均匀地分摊到这个月的每一天中。
         !------------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Compute phenology, growth, mortality, recruitment, disturbance, and check   !
         ! whether we will apply them to the ecosystem or not.                             !
         !---------------------------------------------------------------------------------!
         !! 解析：核心植物学步骤：每天触发一次。调用植被动态总驱动。在此处计算植物群落的物候（如季节性展叶/落叶）、
         !! 生物量分配与生长、自然死亡率、竞争干扰以及幼苗更新。
         call veg_dynamics_driver(new_month,new_year,gr_tfact0,veget_dyn_on)
         ! vegetation_dynamics.f90 里包含了 ED2 模型中所有与植物生长、死亡、更新相关的核心生态学算法。
         !---------------------------------------------------------------------------------!

         !----- First day of a month. -----------------------------------------------------!
         !! 月尺度触发：高性能集群同步与大重置
         if (new_month) then

            !------------------------------------------------------------------------------!
            !      On the monthly timestep we have performed various fusion/fission calls. !
            ! Therefore the var-table's pointer vectors must be updated, and the global    !
            ! definitions of the total numbers must be exported to all nodes.              !
            !      Also, if we do not need to fill the tables until we do I/O, so instead  !
            ! of running this routine every time the demographics change, we set this flag !
            ! and run the routine when the next IO occurs.                                 !
            !------------------------------------------------------------------------------!
            if (nnodetot > 1) then
               if (mynum == 1) write(unit=*,fmt='(a)')                                     &
                                               '-- Monthly node synchronization - waiting'
               if (whos_slow ) then
                  write(unit=*,fmt='(a,1x,i5,1x,a,1x,f7.1)') 'Node',mynum                  &
                                                            ,'time', walltime(wtime_start)
               end if
#if defined(RAMS_MPI)
               call MPI_Barrier(MPI_COMM_WORLD,ierr)
            ! 在经历了一个月的生态学计算后，由于不同斑块上的植物死伤不一，导致各并行节点（Node）手里的计算量产生了严重不均（负载失衡）。
            ! 为了防止跑得快的节点把跑得慢的节点越甩越远，在每月第一天，强制执行 MPI_Barrier 屏障。
            ! 所有人必须在此集合驻留，直到最慢的节点赶到并对齐步调。
#endif
               if (mynum == 1) write(unit=*,fmt='(a)') '-- Synchronized.'
            end if
            

            !! 月尺度更新完成后，升起 filltables 旗帜，提示下一次 I/O 前必须重新映射高维内存表；
            filltables=.true.   ! call filltab_alltypes

            !----- Read new met driver files only if this is the first timestep. ----------!
            !! 同时调用 read_met_drivers() 擦除旧数据，从磁盘里读入下一个月的全新气象强迫场文件；并重新初始化积分器缓冲区。
            call read_met_drivers()
            !----- Re-allocate integration buffer. ----------------------------------------!
            call initialize_rk4patches(.false.)
            
         end if
      end if
      !------------------------------------------------------------------------------------!


      ! 12. 数据输出、状态重置与水文解算
      !------------------------------------------------------------------------------------!
      !      Update the yearly variables.                                                  !
      !------------------------------------------------------------------------------------!
      !! 如果来到了新一年的首日，遍历网格更新并锁死年尺度诊断变量。
      if (analysis_time .and. new_year .and. new_day) then
         do ifm = 1,ngrids
            call update_ed_yearly_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Call the model output driver.                                                  !
      !------------------------------------------------------------------------------------!
      ! 解析：核心输出步骤：调用输出驱动程序。它在内部会集中检查上述传进去的所有布尔标志。
      ! 任何一面“信号旗”为真，它就会立刻把对应的瞬时、日均、月均、或断点状态写成 HDF5 文件输出到标准结果目录。
      call ed_output(observation_time,analysis_time,new_day,new_year,dail_analy_time       &
                    ,mont_analy_time,dcyc_analy_time,annual_time,history_time,dcycle_time)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Write a file with the current history time.                                    !
      !------------------------------------------------------------------------------------!
      ! 解析：如果是输出断点重启文件的时刻，主节点会单独创建或重写一个非常小的纯文本文件（restore_file），
      ! 里面只记录当前精准的年、月、日、时、分。这样下次用户重启模型时，只要读取这个文本，
      ! 就能知道该从哪个 HDF5 恢复文件开始续跑。
      if (history_time .and. mynum == 1) then
         open (unit=18,file=trim(restore_file),form='formatted',status='replace'           &
              ,action='write')
         write(unit=18,fmt=fmtrest) current_time%year,current_time%month,current_time%date &
                                   ,current_time%hour,current_time%min
         close(unit=18,status='keep')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Reset time happens every frqsum.  This is to avoid variables to build up when !
      ! history and analysis are off.  This should be done outside ed_output so I have a   !
      ! chance to copy some of these to BRAMS structures.                                  !
      !------------------------------------------------------------------------------------!
      ! 解析：如果达到了大重置周期（reset_time），清空网格上的高频临时平均统计变量，
      ! 防止在历史记录和快速分析关闭时变量无限堆叠。
      if (reset_time) then
         do ifm=1,ngrids
            call zero_ed_fmean_vars(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Reset inputs to soil carbon.                                                  !
      !------------------------------------------------------------------------------------!
      ! 解析：如果是新的日期，清空土壤碳的输入变量。
      if (new_day) then
         do ifm=1,ngrids
            call zero_litter_inputs(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update the meteorological driver, and the hydrology parameters.               !
      !------------------------------------------------------------------------------------!
      ! 解析：为下一个时间步长做准备：更新下一阶段的气象强迫场数据。如果是新的一月且新的一天，
      ! 则调用 updateHydroParms 依据最新的土壤、植被状况重新校准水文特征参数。
      do ifm=1,ngrids
         call update_met_drivers(edgrid_g(ifm))
      end do
      if (new_day .and. new_month) then
         do ifm = 1,ngrids
            call updateHydroParms(edgrid_g(ifm))
         end do
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update the yearly variables.                                                  !
      !------------------------------------------------------------------------------------!
      !if (analysis_time .and. new_month .and. new_day .and. current_time%month == 6) then
      !   do ifm = 1,ngrids
      !      call update_ed_yearly_vars(edgrid_g(ifm))
      !   end do
      !end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update lateral hydrology.                                                     !
      !------------------------------------------------------------------------------------!
      !! 解析：核心物理步骤 3：调用水文学地表/地下解算器。
      ! 计算多边形网格之间的地下侧向水流径流、地表径流交换以及水分入渗、土壤层间水动力学，并将最终的水文通量和状态结果写出。
      call calcHydroSubsurface()
      call calcHydroSurface()
      call writeHydro()
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update wall time.                                                             !
      !------------------------------------------------------------------------------------!
      ! 解析：当前步长彻底收尾。抓取当前的绝对流逝时间 wtime2 与系统内 CPU 消耗累计 t2。如果开启了横幅打印，主节点会在屏幕上吐出一行极为关键的进度报告。格式如：
      ! Timestep 1440; Sim time 06-15-2026 43200s; Wall 0.045s; CPU 0.038s
      ! 随后代码遇到 end do timestep，立刻带着最新的状态调头回到循环顶部（Step 6），直到将时间推向终点。
      wtime2=walltime(wtime_start)
      call timing(2,t2)
      if (printbanner) then
         write (unit=*,fmt='(a,i10,a,i2.2,a,i2.2,a,i4.4,a,f6.0,2(a,f7.3),a)')              &
             ' Timestep ',istp,'; Sim time  '                                              &
            ,current_time%month,'-',current_time%date,'-',current_time%year,' '            &
            ,current_time%time,'s; Wall',wtime2-wtime1,'s; CPU',t2-t1,'s'
      end if
   end do timestep
   !---------------------------------------------------------------------------------------!

   wtime_tot=walltime(wtime_start)
   write(unit=*,fmt='(a,1x,f10.1,1x,a)') ' === Time integration ends; Total elapsed time=' &
                                        ,wtime_tot," ==="
   ! 解析：当 time >= timmax，跳出整个主循环。计算并向终端报告由于这次模拟所耗费的真实总墙钟时间（秒）。
   ! 子程序安全返回，整个 ED2 模拟主任务顺利完工！
   return
end subroutine ed_model
!==========================================================================================!
!==========================================================================================!
